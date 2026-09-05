import sys
import math
import unittest
from unittest.mock import patch

from analysis import ElectronAnalysis
from analysis import AnalysisConfig
from analysis import inspect_event
from analysis import match_clusters
from analysis import angular_residuals, match_clusters_angular, position_phi
from OLD.validation.fixtures import base_event, cluster, selection_cases
from OLD.validation.reference import run_reference
from OLD.validation.root_stub import ROOT


class AnalysisTests(unittest.TestCase):
    def test_legacy_histograms_counts_and_records(self):
        events = selection_cases()
        with patch.dict(sys.modules, {"ROOT": ROOT}):
            reference = run_reference(events, ROOT)
            analysis = ElectronAnalysis()
            records = [analysis.process_event(event, i) for i, event in enumerate(events)]
            analysis.finalize()
        histograms = {name: value for name, value in reference.items() if name.startswith("h_")}
        self.assertTrue(set(histograms).issubset(analysis.histograms))
        self.assertEqual(set(analysis.histograms) - set(histograms), {
            "h_dphi_noEoP", "h_dphi_withEoP", "h_dtheta_zoom_noEoP", "h_dtheta_zoom_withEoP",
        })
        for name, histogram in histograms.items():
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
                # The distance-failure label was clarified; selection is unchanged.
                if key == "reason" and value == "cluster_too_far_after_EoP":
                    value = "cluster_too_far_after_distance_cut"
                self.assertEqual(value, records[record["event"]][key], (record["event"], key))
        summary = analysis.summary()
        for old_suffix, new_suffix in (("noEoP", "no_eop"), ("withEoP", "with_eop")):
            self.assertEqual(reference[f"fraction_at_cut_{old_suffix}"], summary["at_distance_cut"][new_suffix]["fraction"])
            self.assertEqual(reference[f"n_fail_cut_{old_suffix}"], summary["at_distance_cut"][new_suffix]["legacy_failed_count"])

    def test_eop_is_applied_before_nearest_selection(self):
        nearest, selected, values = match_clusters([cluster(1, 0.5), cluster(30, 1.0)], (1000, 0, 100), 10, 0.7, 1.3)
        self.assertEqual((nearest.index, selected.index), (0, 1))
        self.assertEqual(values, [0.5, 1.0])

    def test_config_changes_selection(self):
        with patch.dict(sys.modules, {"ROOT": ROOT}):
            analysis = ElectronAnalysis(AnalysisConfig(distance_cut_mm=5))
            record = analysis.process_event(base_event())
        self.assertFalse(record["passed_withEoP"])
        self.assertEqual(record["reason"], "cluster_too_far_after_distance_cut")

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
