#!/usr/bin/env python3
"""Inspect one event, including tracks/projections and ECal clusters.

    python3 run_event_display.py sample.root --event 1
    python3 run_event_display.py files.txt --event 1

Edit event data in analysis.py and view styling in histograms.py.
"""

from analysis import AnalysisConfig, display_branches, format_event_details, inspect_event
from analysis_io import event_source, open_run, run_cli, write_json
from histograms import draw_event_display


def execute(args):
    config = AnalysisConfig.load(args.config)
    with open_run(args, config, display_branches(config, args.tracks_only)) as (chain, output):
        if chain.GetEntry(args.event) <= 0:
            raise RuntimeError(f"Could not read event {args.event}")
        record = inspect_event(chain, args.event, config, args.tracks_only)
        record.update(event_source(chain, args.event))
        write_json(output / "event.json", record)
        details = format_event_details(record)
        (output / "event.txt").write_text(details)
        print(details, end="", flush=True)
        files = draw_event_display(record, output)
        write_json(output / "summary.json", {
            "event": args.event, "tracks": len(record["tracks"]),
            "clusters": len(record["clusters"]), **files,
        })
    print("Saved event views:")
    for path in files["plots"]:
        print(f"  {path}")
    print(f"  {files['root_file']}")
    print(f"Event details: {output / 'event.txt'}")
    return output


def main(argv=None):
    return run_cli(execute, "display", argv)


if __name__ == "__main__":
    main()
