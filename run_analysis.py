#!/usr/bin/env python3
"""Analyze an entire ROOT file or file list.

    python3 run_analysis.py sample.root
    python3 run_analysis.py files.txt

Edit cuts/physics in analysis.py and histograms/plots in histograms.py.
"""

import sys

from analysis import AnalysisConfig, ElectronAnalysis, required_branches
from analysis_io import open_run, run_cli


def execute(args):
    config = AnalysisConfig.load(args.config, distance_cut_mm=args.distance_cut_mm,
                                 eop_min=args.eop_min, eop_max=args.eop_max)
    with open_run(args, config, required_branches(config)) as (chain, output):
        analysis = ElectronAnalysis(config)
        analysis.run(chain, output, args.max_events, args.progress_every, event_index=args.event)
        analysis.save_results(output, make_plots=not args.no_plots)
    analysis.print_summary()
    print(f"Completed. Details: {output / 'summary.txt'}")
    return output


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "display":  # Compatibility with the previous command.
        from run_event_display import main as display
        return display(argv[1:])
    return run_cli(execute, "analyze", argv)


if __name__ == "__main__":
    main()
