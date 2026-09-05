#!/usr/bin/env python3
"""One entry point for matching analysis, event displays, and parent checks."""

import argparse
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import sys

from electron_finder.config import AnalysisConfig
from electron_finder.io import build_chain, load_root, prepare_output, resolve_inputs, write_json


PROJECT_DIR = Path(__file__).resolve().parent


def nonnegative_int(value):
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError("must be zero or greater")
    return number


def make_parser(default_command=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    command_parsers = {}
    for name, help_text in (
        ("analyze", "Run the MC electron -> track -> ECAL cluster matching study"),
        ("display", "Inspect one event: tracks/projection points and ECal clusters in XY and RZ"),
        ("parents", "Check the first status-1 electron's first parent"),
    ):
        command = commands.add_parser(name, help=help_text, description=help_text)
        command_parsers[name] = command
        inputs = command.add_mutually_exclusive_group(required=True)
        inputs.add_argument("--input", nargs="+", metavar="ROOT_FILE", help="One or more ROOT paths or URLs")
        inputs.add_argument("--file-list", metavar="TXT", help="ROOT list; relative entries resolve beside this text file")
        command.add_argument("--config", default=str(PROJECT_DIR / "configs" / "baseline.json"), help="JSON settings file")
        command.add_argument("--output", help="New/empty output directory (default: results/<mode>_<timestamp>)")
        if name != "display":
            selection = command.add_mutually_exclusive_group() if name == "analyze" else command
            selection.add_argument("--max-events", type=nonnegative_int, help="Limit the first N entries; 0 creates empty results")
            if name == "analyze":
                selection.add_argument("--event", type=nonnegative_int, help="Run matching histograms for only this zero-based entry; use display for XY/RZ event views")
            command.add_argument("--no-plots", action="store_true", help="Save histograms and records without drawing figures")
        if name == "analyze":
            command.add_argument("--distance-cut", dest="distance_cut_mm", type=float, help="Override distance cut in mm")
            command.add_argument("--eop-min", type=float, help="Override strict lower E/p bound")
            command.add_argument("--eop-max", type=float, help="Override strict upper E/p bound")
            command.add_argument("--progress-every", type=nonnegative_int, default=1000, help="Progress interval; 0 is quiet")
        if name == "display":
            command.add_argument("--event", type=nonnegative_int, default=1, help="Zero-based entry in the whole input chain (default: 1)")
            command.add_argument("--tracks-only", action="store_true", help="Only require projection branches; draw 3D and XY")
    plot = commands.add_parser("plot", help="Redraw matching histograms from a saved results directory")
    plot.add_argument("--results", required=True, help="Completed analysis directory containing analysis.root and config.json")
    plot.add_argument("--output", help="New/empty output directory (default: results/plot_<timestamp>)")
    command_parsers["plot"] = plot
    if default_command is not None:
        parser = command_parsers[default_command]
        parser.prog = Path(sys.argv[0]).name
        parser.set_defaults(command=default_command)
    return parser


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def code_hashes():
    paths = [*sorted(PROJECT_DIR.glob("run_*.py")), *sorted((PROJECT_DIR / "electron_finder").glob("*.py"))]
    return {str(path.relative_to(PROJECT_DIR)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def execute(args):
    overrides = {name: getattr(args, name, None) for name in ("distance_cut_mm", "eop_min", "eop_max")}
    inputs, chain = [], None
    if args.command == "plot":
        results = Path(args.results).expanduser().resolve()
        config = AnalysisConfig.load(results / "config.json")
        if not (results / "analysis.root").is_file():
            raise FileNotFoundError(f"No saved matching analysis in {results}")
    else:
        config = AnalysisConfig.load(args.config, **overrides)
        inputs = resolve_inputs(args.input, args.file_list)
        if args.command == "analyze":
            from electron_finder.analysis import required_branches
            required = required_branches(config)
        else:
            from electron_finder.diagnostics import PARENT_BRANCHES, display_branches
            required = PARENT_BRANCHES if args.command == "parents" else display_branches(config, args.tracks_only)
        chain = build_chain(inputs, config.event_tree, required)
        selected_event = getattr(args, "event", None)
        if selected_event is not None and selected_event >= chain.GetEntries():
            raise ValueError(f"Event {selected_event} is out of range; input has {chain.GetEntries()} entries")
    root = load_root()
    default_name = args.command + "_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    output = prepare_output(args.output or PROJECT_DIR / "results" / default_name)
    manifest = {
        "schema_version": 1, "status": "running", "started_at_utc": timestamp(),
        "command": args.command, "arguments": vars(args), "inputs": inputs,
        "output_dir": str(output), "python_version": sys.version,
        "root_version": str(root.gROOT.GetVersion()), "code_sha256": code_hashes(),
    }
    write_json(output / "config.json", config.to_dict())
    write_json(output / "run.json", manifest)
    (output / "inputs.txt").write_text("\n".join(inputs) + ("\n" if inputs else ""))
    print(f"Results: {output}", flush=True)
    try:
        if args.command == "analyze":
            from electron_finder.analysis import ElectronAnalysis
            analysis = ElectronAnalysis(config)
            summary = analysis.run(chain, output, args.max_events, args.progress_every, event_index=args.event)
            analysis.save_results(output, make_plots=not args.no_plots)
        elif args.command == "parents":
            from electron_finder.diagnostics import run_parent_check
            summary = run_parent_check(chain, output, args.max_events, make_plots=not args.no_plots)
        elif args.command == "display":
            from electron_finder.diagnostics import draw_event_display, inspect_event
            from electron_finder.io import event_source
            from electron_finder.reporting import format_event_details
            if chain.GetEntry(args.event) <= 0:
                raise RuntimeError(f"Could not read event {args.event}")
            record = inspect_event(chain, args.event, config, args.tracks_only)
            record.update(event_source(chain, args.event))
            write_json(output / "event.json", record)
            details = format_event_details(record)
            (output / "event.txt").write_text(details)
            print(details, end="", flush=True)
            files = draw_event_display(record, output)
            summary = {"event": args.event, "tracks": len(record["tracks"]), "clusters": len(record["clusters"]), **files}
            write_json(output / "summary.json", summary)
        else:
            from electron_finder.plotting import plot_saved_results
            summary = plot_saved_results(results, output)
            write_json(output / "summary.json", summary)
        manifest.update(status="complete", finished_at_utc=timestamp())
        write_json(output / "run.json", manifest)
    except BaseException as exc:
        manifest.update(status="interrupted" if isinstance(exc, KeyboardInterrupt) else "failed",
                        finished_at_utc=timestamp(), error=f"{type(exc).__name__}: {exc}")
        write_json(output / "run.json", manifest)
        raise
    if args.command == "analyze":
        analysis.print_summary()
    elif args.command == "display":
        print("Saved event views:")
        for path in summary["plots"]:
            print(f"  {path}")
        print(f"  {summary['root_file']}")
        print(f"Event details: {output / 'event.txt'}")
    else:
        for name, value in summary.get("cutflow", {}).items():
            print(f"  {name}: {value}")
    print(f"Completed. Details: {output / 'summary.json'}")
    return output


def main(argv=None, *, default_command=None):
    parser = make_parser(default_command)
    args = parser.parse_args(argv)
    try:
        execute(args)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(2, f"Error: {exc}\n")
    except KeyboardInterrupt:
        parser.exit(130, "Interrupted. The output run.json records the incomplete run.\n")


if __name__ == "__main__":
    main()
