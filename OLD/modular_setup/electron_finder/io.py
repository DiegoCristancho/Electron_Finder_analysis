"""ROOT input validation, output files, and streaming event records."""

from contextlib import contextmanager
import json
from pathlib import Path


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


@contextmanager
def root_output(path):
    root = load_root()
    output = root.TFile(str(path), "RECREATE")
    if not output or output.IsZombie():
        raise RuntimeError(f"Could not create ROOT output: {path}")
    try:
        output.cd()
        yield output
    finally:
        output.Close()
