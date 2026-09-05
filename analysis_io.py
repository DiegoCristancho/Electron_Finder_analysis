"""Command-line inputs, ROOT reading, and run/event bookkeeping.

Both runners accept a ROOT path or a .txt list as their first argument.
This file is named analysis_io.py because Python already has a built-in io module.
Physics and cuts are in analysis.py; ROOT output is in histograms.py.
"""

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

PROJECT_DIR = Path(__file__).resolve().parent

def load_root():
    try:
        import ROOT
    except ImportError as exc:
        raise RuntimeError(
            "PyROOT is unavailable in this Python environment. Run this command "
            "inside your EIC software environment (eic-shell), which also provides "
            "the EDM4hep/EDM4eic dictionaries. --help works without ROOT."
        ) from exc
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(1110)
    return ROOT


def resolve_inputs(inputs=None, file_list=None):
    """Direct paths use the working directory; list entries use the list directory.

    Blank lines and lines beginning with # are ignored. Remote URLs are preserved.
    """
    if bool(inputs) == bool(file_list):
        raise ValueError("Supply either input ROOT file(s) or one text file list")
    if isinstance(inputs, (str, Path)):
        inputs = [str(inputs)]
    if inputs and len(inputs) == 1 and str(inputs[0]).lower().endswith(".txt"):
        file_list, inputs = inputs[0], None
    elif inputs and any(str(value).lower().endswith(".txt") for value in inputs):
        raise ValueError("Supply a text file list by itself, or supply ROOT files")
    base = Path.cwd()
    if file_list:
        file_list = Path(file_list).expanduser().resolve()
        base = file_list.parent
        inputs = [line.strip() for line in file_list.read_text().splitlines()
                  if line.strip() and not line.lstrip().startswith("#")]
    resolved = []
    for value in inputs:
        value = str(value)
        if "://" not in value:
            path = Path(value).expanduser()
            path = (base / path).resolve()
            if not path.is_file():
                raise FileNotFoundError(f"Input file does not exist: {path}")
            value = str(path)
        if value in resolved:
            raise ValueError(f"Duplicate input would double count events: {value}")
        resolved.append(value)
    if not resolved:
        raise ValueError("The input file list is empty")
    return resolved


def build_chain(inputs, tree_name, required=()):
    """Check every file/tree/branch before processing, including chained files."""
    root = load_root()
    chain = root.TChain(tree_name)
    for filename in inputs:
        source = root.TFile.Open(filename, "READ")
        if not source or source.IsZombie():
            raise RuntimeError(f"Could not open ROOT input: {filename}")
        try:
            tree = source.Get(tree_name)
            if not tree or not tree.InheritsFrom("TTree"):
                raise RuntimeError(f"Missing TTree {tree_name!r} in {filename}")
            missing = [name for name in required if not tree.GetBranch(name)]
            if missing:
                raise RuntimeError(f"Missing branches in {filename}: {', '.join(missing)}")
        finally:
            source.Close()
        if chain.Add(filename) == 0:
            raise RuntimeError(f"Could not add input to TChain: {filename}")
    return chain


def iter_events(chain, max_events=None, *, event_index=None):
    """Yield the first N entries or one selected zero-based global chain entry."""
    if event_index is not None and max_events is not None:
        raise ValueError("Choose either event_index or max_events")
    if max_events is not None and (type(max_events) is not int or max_events < 0):
        raise ValueError("max_events must be a nonnegative integer")
    total = int(chain.GetEntries())
    if event_index is not None:
        if type(event_index) is not int or not 0 <= event_index < total:
            raise ValueError(f"Event {event_index} is out of range; input has {total} entries")
        indices = (event_index,)
    else:
        limit = total if max_events is None else min(total, max_events)
        indices = range(limit)
    for index in indices:
        if chain.GetEntry(index) <= 0:
            raise RuntimeError(f"Could not read event entry {index}")
        yield index, chain


def event_source(chain, index):
    current = chain.GetCurrentFile()
    return {
        "event": index,
        "source_file": current.GetName() if current else None,
        "source_entry": int(chain.GetTree().GetReadEntry()),
    }


def prepare_output(path):
    path = Path(path).expanduser().resolve()
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise FileExistsError(f"Output directory is not empty; choose a new run directory: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


class EventWriter:
    """Write records immediately, so memory does not grow with the event count."""

    def __init__(self, directory):
        self.directory = Path(directory)

    def __enter__(self):
        self.events = (self.directory / "events.jsonl").open("w")
        self.failures = (self.directory / "failures.jsonl").open("w")
        return self

    def write(self, record):
        line = json.dumps(record, allow_nan=False) + "\n"
        self.events.write(line)
        if record.get("reason") is not None:
            self.failures.write(line)

    def __exit__(self, *exc):
        self.events.close()
        self.failures.close()


# Command-line input shared by the two runners.
def nonnegative_int(value):
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError("must be zero or greater")
    return number


def make_parser(mode="analyze"):
    description = ("Analyze every event and print the electron-matching summary."
                   if mode == "analyze" else "Draw one event's tracks and ECal clusters.")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("files", nargs="*", metavar="FILE", help="A ROOT file, several ROOT files, or one .txt list")
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument("--input", nargs="+", metavar="FILE", help="Alternative to the positional input")
    inputs.add_argument("--file-list", metavar="TXT", help="Alternative for a .txt input list")
    parser.add_argument("--output", help="New/empty result directory; otherwise results/<mode>_<timestamp>")
    parser.add_argument("--config", help="Optional saved JSON settings; defaults are in analysis.py")
    if mode == "analyze":
        selection = parser.add_mutually_exclusive_group()
        selection.add_argument("--max-events", type=nonnegative_int, help="Process only the first N events")
        selection.add_argument("--event", type=nonnegative_int, help="Run matching for one zero-based entry")
        parser.add_argument("--distance-cut", dest="distance_cut_mm", type=float, help="Override the distance cut in mm")
        parser.add_argument("--eop-min", type=float, help="Override the strict lower E/p bound")
        parser.add_argument("--eop-max", type=float, help="Override the strict upper E/p bound")
        parser.add_argument("--no-plots", action="store_true", help="Save histograms/records without drawing PNG figures")
        parser.add_argument("--progress-every", type=nonnegative_int, default=1000, help="Print progress every N events; 0 is quiet")
    else:
        parser.add_argument("--event", type=nonnegative_int, default=1, help="Zero-based entry to display (default: 1)")
        parser.add_argument("--tracks-only", action="store_true", help="Draw projection points without requiring clusters or MC")
    parser.set_defaults(command=mode)
    return parser


def run_cli(execute, mode, argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    # Keep previous commands such as `run_analysis.py analyze --input ...` working.
    if argv and argv[0] == mode:
        argv.pop(0)
    parser = make_parser(mode)
    args = parser.parse_args(argv)
    if args.files and (args.input or args.file_list):
        parser.error("Positional input is not allowed with --input or --file-list")
    if not (args.files or args.input or args.file_list):
        parser.error("Supply a ROOT file or a .txt list")
    if args.files:
        args.input = args.files
    del args.files
    try:
        return execute(args)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(2, f"Error: {exc}\n")
    except KeyboardInterrupt:
        parser.exit(130, "Interrupted. See the output run.json for this run's status.\n")


@contextmanager
def open_run(args, config, required):
    """Validate inputs, prepare the run folder, and record completion/failure."""
    inputs = resolve_inputs(args.input, args.file_list)
    chain = build_chain(inputs, config.event_tree, required)
    selected = getattr(args, "event", None)
    if selected is not None and selected >= chain.GetEntries():
        raise ValueError(f"Event {selected} is out of range; input has {chain.GetEntries()} entries")
    name = args.command + "_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    output = prepare_output(args.output or PROJECT_DIR / "results" / name)
    sources = ("run_analysis.py", "run_event_display.py", "analysis_io.py", "analysis.py", "histograms.py")
    manifest = {
        "schema_version": 1, "status": "running",
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": args.command, "arguments": vars(args), "inputs": inputs,
        "output_dir": str(output), "python_version": sys.version,
        "root_version": str(load_root().gROOT.GetVersion()),
        "code_sha256": {name: hashlib.sha256((PROJECT_DIR / name).read_bytes()).hexdigest() for name in sources},
    }
    write_json(output / "config.json", config.to_dict())
    write_json(output / "run.json", manifest)
    (output / "inputs.txt").write_text("\n".join(inputs) + "\n")
    print(f"Results: {output}", flush=True)
    try:
        yield chain, output
    except BaseException as exc:
        manifest.update(status="interrupted" if isinstance(exc, KeyboardInterrupt) else "failed",
                        error=f"{type(exc).__name__}: {exc}")
        raise
    else:
        manifest["status"] = "complete"
    finally:
        manifest["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
        write_json(output / "run.json", manifest)
