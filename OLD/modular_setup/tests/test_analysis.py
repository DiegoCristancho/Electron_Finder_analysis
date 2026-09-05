import sys
import unittest
from unittest.mock import patch

from electron_finder.analysis import ElectronAnalysis
from electron_finder.config import AnalysisConfig
from electron_finder.diagnostics import inspect_event
from electron_finder.matching import match_clusters
from tests.fixtures import base_event, cluster, selection_cases
from tests.reference import run_reference
from tests.root_stub import ROOT


class AnalysisTests(unittest.TestCase):
    def test_legacy_histograms_counts_and_records(self):
        events = selection_cases()
        with patch.dict(sys.modules, {"ROOT": ROOT}):
            reference = run_reference(events, ROOT)
            analysis = ElectronAnalysis()
            records = [analysis.process_event(event, i) for i, event in enumerate(events)]
            analysis.finalize()
        histograms = {name: value for name, value in reference.items() if name.startswith("h_")}
        self.assertEqual(set(histograms), set(analysis.histograms))
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
        self.assertEqual(record["reason"], "cluster_too_far_after_EoP")

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


if __name__ == "__main__":
    unittest.main()
