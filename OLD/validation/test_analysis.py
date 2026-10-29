import sys
import math
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import MagicMock, patch

from analysis import ElectronAnalysis
from analysis import AnalysisConfig
from analysis import inspect_event
from analysis import match_clusters
from analysis import angular_residuals, format_analysis_summary, match_clusters_angular, momentum_eta, position_phi
from histograms import plot_analysis
from OLD.validation.fixtures import base_event, cluster, selection_cases
from OLD.validation.reference import run_reference
from OLD.validation.root_stub import ROOT


class AnalysisTests(unittest.TestCase):
    def test_legacy_histograms_counts_and_records(self):
        events = selection_cases()
        with patch.dict(sys.modules, {"ROOT": ROOT}):
            reference = run_reference(events, ROOT)
            analysis = ElectronAnalysis(AnalysisConfig(distance_cut_mm=87))
            records = [analysis.process_event(event, i) for i, event in enumerate(events)]
            analysis.finalize()
        histograms = {name: value for name, value in reference.items() if name.startswith("h_")}
        self.assertTrue(set(histograms).issubset(analysis.histograms))
        self.assertEqual(set(analysis.histograms) - set(histograms), {
            "h_dphi_noEoP", "h_dphi_withEoP", "h_dtheta_zoom_noEoP", "h_dtheta_zoom_withEoP",
            "h_eop_per_event", "h_lost_track_mc_eta",
        })
        for name, histogram in histograms.items():
            if name in {"h_fail_recoE", "h_fail_mcE"}:
                # Failure-energy histograms now follow the angular selection.
                continue
            with self.subTest(histogram=name):
                result = analysis.histograms[name]
                self.assertEqual(histogram.contents, result.contents)
                self.assertEqual(histogram.entries, result.entries)
                self.assertEqual(histogram.fills, result.fills)
        counts = {
            "n_total_events": "events_processed", "n_all_status1_electrons": "status1_electrons",
            "n_events_with_scattered_electron": "selected_electrons", "n_scattered_with_reco": "with_reco",
            "n_scattered_with_track": "with_track", "n_scattered_with_projection": "with_projection",
            "n_scattered_with_ecal_point": "with_ecal_point",
            "n_scattered_with_cluster_noEoP": "with_cluster_no_eop",
            "n_scattered_with_cluster_withEoP": "with_cluster_with_eop",
            "n_pass_cut_noEoP": "pass_distance_no_eop", "n_pass_cut_withEoP": "pass_distance_with_eop",
        }
        for old, new in counts.items():
            self.assertEqual(reference[old], analysis.counts[new], old)
        for record in reference["scattered_records"]:
            for key, value in record.items():
                self.assertEqual(value, records[record["event"]][key], (record["event"], key))
        for record in reference["failed_records"]:
            for key, value in record.items():
                # Distance-only failures are no longer angular-selection failures.
                if key == "reason" and value == "cluster_too_far_after_EoP":
                    value = None
                if key == "reason" and value == "no_cluster_after_EoP":
                    value = ("clusters_fail_EoP" if records[record["event"]]["has_cluster_noEoP"]
                             else "no_cluster_available")
                self.assertEqual(value, records[record["event"]][key], (record["event"], key))
        summary = analysis.summary()
        for old_suffix, new_suffix in (("noEoP", "no_eop"), ("withEoP", "with_eop")):
            self.assertEqual(reference[f"fraction_at_cut_{old_suffix}"], summary["at_distance_cut"][new_suffix]["fraction"])
            self.assertEqual(reference[f"n_fail_cut_{old_suffix}"], summary["at_distance_cut"][new_suffix]["legacy_failed_count"])

    def test_eop_is_applied_before_nearest_selection(self):
        nearest, selected, values = match_clusters([cluster(1, 0.5), cluster(30, 1.0)], (1000, 0, 100), 10, 0.7, 1.3)
        self.assertEqual((nearest.index, selected.index), (0, 1))
        self.assertEqual(values, [0.5, 1.0])

    def test_summary_separates_missing_clusters_from_eop_failures(self):
        no_cluster = base_event()
        no_cluster.EcalBarrelClusters = []
        fails_eop = base_event()
        fails_eop.EcalBarrelClusters = [cluster(1, 0.5)]
        passes_eop = base_event()
        with patch.dict(sys.modules, {"ROOT": ROOT}):
            analysis = ElectronAnalysis()
            records = [analysis.process_event(event, index) for index, event in
                       enumerate((no_cluster, fails_eop, passes_eop))]
            summary = analysis.summary()
        self.assertEqual([record["reason"] for record in records],
                         ["no_cluster_available", "clusters_fail_EoP", None])
        self.assertEqual(summary["eop_diagnostics"], {
            "no_cluster": 1, "clusters_fail_eop": 1, "clusters_pass_eop": 1,
        })
        report = format_analysis_summary(summary)
        self.assertIn("No reconstructed ECAL cluster                 : 1", report)
        self.assertIn("Cluster(s) present, all fail E/p              : 1", report)

    def test_lost_track_eta_includes_missing_reco_and_empty_track_range(self):
        missing_reco = base_event()
        missing_reco._ReconstructedParticleAssociations_sim = []
        missing_track = base_event()
        missing_track.ReconstructedParticles[0].tracks_end = 0
        good = base_event()
        with patch.dict(sys.modules, {"ROOT": ROOT}):
            analysis = ElectronAnalysis()
            records = [analysis.process_event(event, index) for index, event in
                       enumerate((missing_reco, missing_track, good))]
            summary = analysis.summary()
        self.assertEqual(analysis.counts["lost_tracks"], 2)
        self.assertEqual(analysis.histograms["h_lost_track_mc_eta"].fills, [0.0, 0.0])
        self.assertEqual([record["lost_track"] for record in records], [True, True, False])
        self.assertEqual([record["mc_eta"] for record in records], [0.0, 0.0, 0.0])
        self.assertEqual(summary["lost_track_eta_entries"], 2)

    def test_eop_per_event_uses_one_nearest_cluster_value_before_cut(self):
        event = base_event()
        event.EcalBarrelClusters = [cluster(1, 0.5), cluster(30, 1.0)]
        with patch.dict(sys.modules, {"ROOT": ROOT}):
            analysis = ElectronAnalysis()
            record = analysis.process_event(event)
            summary = analysis.summary()
        self.assertEqual(record["eop_nearest_cluster"], 0.5)
        self.assertEqual(record["best_eop"], 1.0)
        self.assertEqual(analysis.histograms["h_eop_per_event"].fills, [0.5])
        self.assertEqual(summary["eop_per_event_entries"], 1)

    def test_momentum_eta_is_truth_momentum_pseudorapidity(self):
        electron = base_event().MCParticles[0]
        electron.momentum = type(electron.momentum)(x=3.0, y=4.0, z=5.0)
        self.assertAlmostEqual(momentum_eta(electron), math.asinh(1.0))
        electron.momentum.x = electron.momentum.y = 0.0
        self.assertIsNone(momentum_eta(electron))

    def test_diagnostic_figures_are_separate_and_show_eop_cut_limits(self):
        with patch.dict(sys.modules, {"ROOT": ROOT}):
            analysis = ElectronAnalysis()
            analysis.finalize()
            histogram_names = analysis.histograms
        histograms = {name: MagicMock() for name in histogram_names}
        histograms["h_eop_per_event"].GetMaximum.return_value = 4.0
        fake_root = MagicMock(kBlack=1, kRed=2, kBlue=3)
        canvases = {}

        def make_canvas(name, *_):
            canvases[name] = MagicMock()
            return canvases[name]

        fake_root.TCanvas.side_effect = make_canvas
        config = AnalysisConfig()
        with TemporaryDirectory() as directory, patch("histograms.load_root", return_value=fake_root):
            result = plot_analysis(histograms, config, directory)
        self.assertEqual(set(result), {
            "distance_no_eop", "distance_with_eop", "fraction_vs_distance",
            "failed_energies", "parent_status", "parent_pdg", "e_over_p",
            "position_differences", "e_over_p_per_event", "lost_track_eta",
        })
        self.assertEqual(
            [call.args for call in fake_root.TLine.call_args_list],
            [(config.eop_min, 0.0, config.eop_min, 4.2),
             (config.eop_max, 0.0, config.eop_max, 4.2)],
        )
        saved = {Path(canvas.SaveAs.call_args.args[0]).name for canvas in canvases.values()}
        self.assertIn("e_over_p_per_event.png", saved)
        self.assertIn("lost_track_eta.png", saved)

    def test_distance_failure_is_not_an_angular_failure(self):
        with patch.dict(sys.modules, {"ROOT": ROOT}):
            analysis = ElectronAnalysis(AnalysisConfig(distance_cut_mm=5))
            record = analysis.process_event(base_event())
        self.assertFalse(record["passed_withEoP"])
        self.assertTrue(record["passed_angular_withEoP"])
        self.assertIsNone(record["reason"])

    def test_angular_failure_reports_event_reason_and_cluster_count(self):
        event = base_event()
        event.EcalBarrelClusters = [cluster(1), cluster(2)]
        for candidate in event.EcalBarrelClusters:
            candidate.position.y = 100
        with patch.dict(sys.modules, {"ROOT": ROOT}):
            analysis = ElectronAnalysis()
            record = analysis.process_event(event)
            summary = analysis.summary()
        self.assertFalse(record["passed_angular_withEoP"])
        self.assertEqual(record["reason"], "all_clusters_outside_theta_phi_cuts")
        self.assertEqual(record["angular_clusters_after_eop"], 2)
        self.assertEqual(record["angular_clusters_outside_cuts"], 2)
        self.assertEqual(summary["failure_reasons"]["all_clusters_outside_theta_phi_cuts"], 1)
        self.assertEqual(summary["cutflow"]["angular_clusters_outside_cuts"], 2)
        report = format_analysis_summary(summary)
        self.assertIn("Clusters outside theta/phi cuts               : 2", report)

    def test_empty_sample_is_defined(self):
        with patch.dict(sys.modules, {"ROOT": ROOT}):
            analysis = ElectronAnalysis()
            analysis.finalize()
        self.assertEqual(analysis.summary()["at_distance_cut"]["with_eop"]["fraction"], 0)

    def test_event_display_records_are_copies(self):
        event = base_event()
        record = inspect_event(event, 12, AnalysisConfig())
        event._CalorimeterTrackProjections_points[0].position.x = -1
        self.assertEqual(record["tracks"][0]["points"][0]["x"], 1000)
        self.assertEqual(record["tracks"][0]["truth_pdg"], 11)

    def test_angular_selection_ignores_radius_and_compares_efficiencies(self):
        event = base_event()
        # The nearby cluster fails the angular window. The distant cluster has
        # identical angles, so it must be found independently of 3D matching.
        distant = cluster(1000)
        distant.position.z = 200
        nearby = cluster(0)
        nearby.position.y = 10
        event.EcalBarrelClusters = [nearby, distant]
        missing_reco = base_event()
        missing_reco._ReconstructedParticleAssociations_sim = []
        with patch.dict(sys.modules, {"ROOT": ROOT}):
            analysis = ElectronAnalysis(AnalysisConfig(distance_cut_mm=5, theta_cut_rad=0.002, phi_min_rad=-0.0012, phi_max_rad=0.0012))
            record = analysis.process_event(event)
            analysis.process_event(missing_reco)
            self.assertFalse(record["passed_withEoP"])
            self.assertTrue(record["passed_angular_withEoP"])
            self.assertEqual(record["angular_match"]["cluster_index"], 1)
            self.assertEqual(analysis.electron_efficiency(), 0)
            self.assertEqual(analysis.angular_electron_efficiency(), 0.5)
            self.assertEqual(analysis.counts["pass_angular_only"], 1)

    def test_phi_wrap_and_strict_angular_and_eop_limits(self):
        residuals = angular_residuals((-1, -0.0001, 0), (-1, 0.0001, 0))
        self.assertAlmostEqual(residuals[1], 0.0002, places=10)
        self.assertIsNone(position_phi(0, 0, 1))
        candidate = cluster(0)
        candidate.position.x, candidate.position.y, candidate.position.z = math.cos(0.002), math.sin(0.002), 0
        limit = angular_residuals((candidate.position.x, candidate.position.y, 0), (1, 0, 0))[1]
        self.assertIsNone(match_clusters_angular([candidate], (1, 0, 0), 10, 0.7, 1.3, 0.002, -limit, limit))
        self.assertIsNotNone(match_clusters_angular([candidate], (1, 0, 0), 10, 0.7, 1.3, 0.002, -limit, limit * 2))
        candidate.energy = 7
        self.assertIsNone(match_clusters_angular([candidate], (1, 0, 0), 10, 0.7, 1.3, 0.002, -limit, limit * 2))

    def test_angular_config_and_empty_efficiency(self):
        for value in (0, -1, math.inf, math.nan):
            with self.assertRaises(ValueError):
                AnalysisConfig(theta_cut_rad=value)
            with self.assertRaises(ValueError):
                AnalysisConfig(phi_min_rad=value, phi_max_rad=value)
        with patch.dict(sys.modules, {"ROOT": ROOT}):
            self.assertIsNone(ElectronAnalysis().angular_electron_efficiency())

    def test_asymmetric_phi_window_and_theta_boundaries(self):
        config = AnalysisConfig()
        self.assertEqual((config.theta_cut_rad, config.phi_min_rad, config.phi_max_rad),
                         (0.02, -0.012, 0.005))
        for theta, phi, passes in (
            (0.019, -0.011, True), (-0.019, 0.004, True),
            (0, 0.006, False), (0, -0.013, False),
            (0, -0.012, False), (0, 0.005, False),
            (0.02, 0, False), (-0.02, 0, False),
        ):
            with self.subTest(theta=theta, phi=phi), patch(
                "analysis.angular_residuals", return_value=(theta, phi)
            ):
                result = match_clusters_angular([cluster(0)], (1000, 0, 100), 10,
                                                0.7, 1.3, config.theta_cut_rad,
                                                config.phi_min_rad, config.phi_max_rad)
                self.assertEqual(result is not None, passes)


if __name__ == "__main__":
    unittest.main()
