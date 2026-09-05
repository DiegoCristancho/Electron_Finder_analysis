"""Compare a small real ROOT sample against the original matching calculation.

Run from the project directory inside eic-shell:
    python3 -m tests.compare_root_sample --file-list inputs/local_example.txt --max-events 100

This does not execute the archived script's file writes or plotting. It checks every
histogram bin (including under/overflow), errors, entries, counts and records.
"""

import argparse
import math

from electron_finder.analysis import ElectronAnalysis, required_branches
from electron_finder.config import AnalysisConfig
from electron_finder.io import build_chain, iter_events, load_root, resolve_inputs
from run_analysis import nonnegative_int
from tests.reference import run_reference


def compare_values(expected, actual, description):
    if isinstance(expected, float):
        equal = isinstance(actual, (int, float)) and math.isclose(expected, actual, rel_tol=1e-7, abs_tol=1e-10)
    else:
        equal = expected == actual
    if not equal:
        raise AssertionError(f"{description}: original={expected!r}, refactored={actual!r}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--input", nargs="+")
    inputs.add_argument("--file-list")
    parser.add_argument("--max-events", type=nonnegative_int, default=100)
    args = parser.parse_args()
    root = load_root()
    config = AnalysisConfig()
    filenames = resolve_inputs(args.input, args.file_list)
    chain = build_chain(filenames, config.event_tree, required_branches(config))
    root.gROOT.cd()
    reference = run_reference((event for _, event in iter_events(chain, args.max_events)), root)
    original_histograms = {name: value for name, value in reference.items() if name.startswith("h_")}
    for histogram in original_histograms.values():
        histogram.SetDirectory(0)

    analysis = ElectronAnalysis(config)
    records = {}
    for index, event in iter_events(chain, args.max_events):
        records[index] = analysis.process_event(event, index)
    analysis.finalize()
    if set(original_histograms) != set(analysis.histograms):
        raise AssertionError("The original and refactored histogram names differ")
    for name, original in original_histograms.items():
        histogram = analysis.histograms[name]
        compare_values(original.GetNbinsX(), histogram.GetNbinsX(), f"{name} bin count")
        compare_values(original.GetEntries(), histogram.GetEntries(), f"{name} entries")
        compare_values(original.GetMean(), histogram.GetMean(), f"{name} mean")
        for i in range(original.GetNbinsX() + 2):
            compare_values(original.GetBinContent(i), histogram.GetBinContent(i), f"{name} bin {i}")
            compare_values(original.GetBinError(i), histogram.GetBinError(i), f"{name} bin {i} error")
    count_names = {
        "n_total_events": "events_processed", "n_all_status1_electrons": "status1_electrons",
        "n_events_with_scattered_electron": "selected_electrons", "n_scattered_with_reco": "with_reco",
        "n_scattered_with_track": "with_track", "n_scattered_with_projection": "with_projection",
        "n_scattered_with_ecal_point": "with_ecal_point", "n_scattered_with_cluster_noEoP": "with_cluster_no_eop",
        "n_scattered_with_cluster_withEoP": "with_cluster_with_eop",
        "n_pass_cut_noEoP": "pass_distance_no_eop", "n_pass_cut_withEoP": "pass_distance_with_eop",
    }
    for old_name, new_name in count_names.items():
        compare_values(reference[old_name], analysis.counts[new_name], old_name)
    for group in ("scattered_records", "failed_records"):
        for old_record in reference[group]:
            for key, value in old_record.items():
                compare_values(value, records[old_record["event"]][key], f"event {old_record['event']} {key}")
    for old_suffix, new_suffix in (("noEoP", "no_eop"), ("withEoP", "with_eop")):
        result = analysis.summary()["at_distance_cut"][new_suffix]
        compare_values(reference[f"fraction_at_cut_{old_suffix}"], result["fraction"], f"{old_suffix} fraction")
        compare_values(reference[f"n_fail_cut_{old_suffix}"], result["legacy_failed_count"], f"{old_suffix} legacy failures")
    print(f"PASS: {len(records)} events, all {len(original_histograms)} histograms, counters and legacy event records agree.")


if __name__ == "__main__":
    main()
