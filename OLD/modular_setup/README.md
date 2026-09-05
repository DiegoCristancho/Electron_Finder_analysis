# Electron finder analysis

Run the full-file matching analysis with `run_analysis.py analyze`, then inspect
individual events with `run_event_display.py`. Reusable calculations live in `electron_finder/`; inputs and
cuts are supplied separately. The original scripts in `OLD/other/` remain the
reference for the refactor.

## Your two-step workflow

From this directory inside your EIC environment:

```bash
# 1. Analyze every event and print the detailed Trial.py-style summary.
python3 run_analysis.py analyze --file-list inputs/local_example.txt

# 2. Inspect event 1 with the projection/track and ECal cluster views from Analysis.py.
python3 run_event_display.py --file-list inputs/local_example.txt --event 1
```

The two entry scripts share the same package. `run_event_display.py` selects the
existing `run_analysis.py display` command directly, so both display commands
produce the same plots and event details. Use `--input filename.root` instead of
`--file-list` when selecting one ROOT file beside the runners.

## Start with a small sample

Use Python 3.9 or newer **inside your EIC software environment**, with PyROOT and
the EDM4hep/EDM4eic dictionaries available. No additional Python packages are
required. The system Python can run the help commands and unit tests, but cannot
read the analysis inputs without that environment.

From a terminal, enter your existing environment:

```bash
/home/diego/eic/eic-shell
cd /home/diego/eic/Analysis_epic/Electron_finder
```

Then try 100 events from the local file already referenced in
`inputs/local_example.txt`:

```bash
python3 run_analysis.py analyze \
  --file-list inputs/local_example.txt \
  --max-events 100 \
  --output results/first_100
```

If you are already in eic-shell, run the analysis command there directly. Use a
new output directory for each run, or omit `--output` to generate a timestamped
directory under this project's `results/`. Existing nonempty directories are
rejected so a later run cannot replace an earlier result.

Remove `--max-events 100` to process the whole input. Add `--no-plots` to save
histograms and records now and draw figures later.

To draw event 1 of a ROOT file beside `run_analysis.py`, run from this
directory inside your EIC environment:

```bash
python3 run_event_display.py \
  --input 'pythia8NCDIS_26_718x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_1.0008.eicrecon.edm4eic.root' \
  --event 1
```

Event numbering is zero-based: `--event 1` selects the second entry;
`--event 0` selects the first. The selected global entry index is retained in
the output records. The separate `run_analysis.py analyze --event 1` option
runs the matching histograms and selection summary for that single entry;
use `run_event_display.py` for the geometric track/cluster views. In `analyze`,
use either `--event` or `--max-events`, not both. Omitting both processes every
event in the inputs.

Each matching run prints the detailed `Trial.py`-style summary: counts through
the matching stages, mean distances, position differences, parent statuses,
plateau information, passing/failing counts, cut fractions, and failure reasons.
The same report is saved as `summary.txt` alongside the JSON summary. This report
is printed even with `--no-plots` or `--progress-every 0`.

## Change files or cuts

One file, several files, and a text file list are supported:

```bash
python3 run_analysis.py analyze --input /path/to/sample.root
python3 run_analysis.py analyze --input /path/to/first.root /path/to/second.root
python3 run_analysis.py analyze --file-list inputs/dis_18x275_q2_100_26_02.txt
```

Direct relative input paths are relative to your working directory. Relative
paths **inside a file list are relative to that list's directory**. Blank lines
and full-line `#` comments are ignored. `root://` URLs are passed through to ROOT;
remote inputs need the corresponding network access. Duplicate inputs are
rejected to prevent double counting.

Edit `configs/baseline.json` or copy it to a named configuration for a study:

```bash
python3 run_analysis.py analyze --file-list inputs/local_example.txt \
  --config configs/baseline.json --distance-cut 60 --eop-min 0.8 --eop-max 1.2
```

Command-line cuts override the JSON file. The effective settings are saved with
each run. Unknown JSON keys and invalid cuts produce an error.

| Setting | Baseline | Meaning |
| --- | --- | --- |
| `event_tree` | `events` | Input TTree |
| `cluster_collection` | `EcalBarrelClusters` | Cluster branch |
| `ecal_system` | `101` | Projection point system ID |
| `distance_cut_mm` | `87.0` | Strict 3D distance cut |
| `eop_min`, `eop_max` | `0.7`, `1.3` | Strict E/p bounds |
| `cut_scan_bins` | `100` | Bins for the distance-cut scan |
| `cut_scan_max_mm` | `200.0` | Upper edge of the cut scan |

Other histogram ranges remain those in the original analysis and can be edited
in `electron_finder/histograms.py`.

## Inspect events, parents, and plots

```bash
# Event numbering is zero-based across the complete input chain.
python3 run_event_display.py --file-list inputs/local_example.txt \
  --event 12 --output results/event12

# Projection-only files: no MC, reconstructed particles, or clusters required.
python3 run_event_display.py --file-list inputs/local_example.txt \
  --event 12 --tracks-only --output results/event12_tracks

# First PDG=11/status=1 electron's parent, without the analysis parent veto.
python3 run_analysis.py parents --file-list inputs/local_example.txt \
  --max-events 100 --output results/parents_100

# Redraw matching histograms using the settings saved with that run.
python3 run_analysis.py plot --results results/first_100 --output results/replotted
```

Display mode prints the original `Analysis.py` event details: collection counts,
track identities, reconstructed-particle and truth matches, each projection
point's system/surface and XYZ/R coordinates, and ECal cluster coordinates.
It displays every projected track and cluster in the chosen event, independently
of the scattered-electron selection used by the full-file matching study.

Its output directory contains:

```text
event.json                    # Copied event coordinates and associations
event.txt                     # The same detailed printout shown in the terminal
event_display.root            # All four canvases, with their drawing objects
plots/event_3d.png             # 3D track projection points
plots/event_xy.png             # XY track projection points
plots/event_xy_clusters.png    # XY tracks/projection points + ECal clusters
plots/event_rz_clusters.png    # RZ tracks/projection points + ECal clusters
```

In RZ, Z is the horizontal axis and R = sqrt(x*x+y*y) is the vertical axis, as
in the original script. The projection-only option saves 3D and XY. Empty events
have defined plotting ranges. Saved plot paths are printed at the end of the run
and recorded in `summary.json`. Config, input paths, and run status are saved too.

The plot command reads `analysis.root` and `config.json` from the result directory.
It does not reopen the original event files. Change figure styling in
`electron_finder/plotting.py` and rerun `plot` to a new directory.

## What a matching run saves

```text
results/first_100/
├── analysis.root       # Original histogram names, plus canvases when requested
├── plots/              # Eight PNG figures, unless --no-plots
├── events.jsonl        # One record for every processed event
├── failures.jsonl      # Events that did not pass the full selection
├── summary.json        # Counts by selection step, failures, fractions, plateaus
├── summary.txt         # The same detailed report printed in the terminal
├── config.json         # Effective settings, including command-line overrides
├── inputs.txt          # Resolved input paths in processing order
└── run.json            # Run status, arguments, versions, timestamps, code hashes
```

Records include the global event index, source file and local entry, selected MC
and reconstructed particle indices, track identity, projection and cluster
indices, energies and momentum in GeV, distances in mm, E/p, and a failure reason.
They contain copied Python values and are streamed to disk rather than retained
for every event in memory. `failures.jsonl` also includes events with no selected
electron, which the original `failed_records` list did not include.

`run.json` records `complete`, `failed`, or `interrupted`; a forcibly killed
process can leave it at `running`. Only complete runs should be used as final
results. Files, trees, and required branches are validated before processing;
missing branches produce an error instead of silently skipping events.

## Analysis definitions retained from the reference

This is currently a **truth-seeded electron matching study**. It starts with the
first MC particle with `PDG == 11`, `generatorStatus == 1`, and first-parent
status different from 2. The original convention lets a missing or invalid
parent index pass this veto because its status is `None`. Parent checking uses
the first status-1 electron without that veto, intentionally.

The matching then uses the first MC/reco association, the first track on that
reconstructed particle, the first matching projection, and the first projection
point with the chosen system ID. Track identity includes both collection ID and
index. Cluster matching minimizes 3D distance. For the E/p branch, the strict
E/p window is applied **before** selecting the nearest cluster. Equal distances
retain the first cluster. Both E/p and distance bounds are strict.

The histogram fractions and reported matching fraction are conditional on having
a reconstructed ECAL projection point. They are not fractions of all generated
electrons. `summary.json` names that denominator explicitly. The original printed
failure count included missing ECAL points despite excluding them from the
fraction denominator; this value is retained as `legacy_failed_count`, alongside
the consistent `failed_among_ecal_projected` count.

Position-based delta-r uses the spherical radius `sqrt(x*x+y*y+z*z)`. The RZ event
display uses the cylindrical radius `sqrt(x*x+y*y)`. These distinct original
definitions are preserved. Histogram division and its ROOT error treatment are
also preserved; no new statistical prescription is introduced in this refactor.

## Where to work on the code

| Module | Responsibility |
| --- | --- |
| `run_analysis.py` | Command-line entry point and run bookkeeping |
| `run_event_display.py` | Dedicated entry point for one event's XY/RZ and 3D views |
| `electron_finder/config.py` | Validated settings |
| `electron_finder/io.py` | Input lists, ROOT chain, output files and records |
| `electron_finder/physics.py` | Momentum, energy and position calculations |
| `electron_finder/truth.py` | MC electron selection and parent information |
| `electron_finder/matching.py` | Associations, projections and nearest clusters |
| `electron_finder/analysis.py` | One `ElectronAnalysis` instance per run |
| `electron_finder/histograms.py` | Histogram names, titles and binning |
| `electron_finder/plotting.py` | Matching figures and replotting saved results |
| `electron_finder/reporting.py` | Trial.py summary and Analysis.py event details |
| `electron_finder/diagnostics.py` | Event displays and parent checks |

For use from another Python script:

```python
from electron_finder.config import AnalysisConfig
from electron_finder.analysis import ElectronAnalysis, required_branches
from electron_finder.io import build_chain, prepare_output, resolve_inputs, write_json

config = AnalysisConfig.load("configs/baseline.json")
inputs = resolve_inputs(file_list="inputs/local_example.txt")
chain = build_chain(inputs, config.event_tree, required_branches(config))
output = prepare_output("results/python_run")
write_json(output / "config.json", config.to_dict())

analysis = ElectronAnalysis(config)
summary = analysis.run(chain, output, max_events=100)
analysis.save_results(output)
```

Use `process_event(event, index)` to inspect a single returned record in Python.
`summary()` finalizes the run, so create another instance to process more events
after finalization. Call `analysis.print_summary()` for the detailed terminal
report when using the class directly. To process one entry through this API,
use `analysis.run(chain, output, event_index=1)` instead of `max_events`.
The CLI additionally saves provenance and run status.

## Validation

Without ROOT:

```bash
python3 -m unittest discover -v
```

The tests execute the original calculation from
`OLD/other/Distances_Between_cluster_track.py` on synthetic events and compare every
histogram fill/bin, selection counter, legacy event/failure record, and chosen-cut
fraction. Cases cover missing associations, missing ECAL points, no clusters,
zero momentum, strict boundaries, E/p ordering, ties, invalid/missing parents,
and histogram overflow. The histogram recorder verifies the Python calculation;
it does not validate ROOT drawing, serialization, or EDM dictionaries.
Additional checks compare every original `Trial.py` report value and verify that
single-event selection reads only the requested entry.
The event-display tests also check full runner dispatch, track/cluster coordinates
in XY and cylindrical RZ, and projection-only inputs. Graphics calls are recorded
for these tests; they do not render real ROOT images.

For a real-file comparison **inside eic-shell**:

```bash
python3 -m tests.compare_root_sample --file-list inputs/local_example.txt --max-events 100
```

This runs the original calculation and the refactored calculation on the same
entries, comparing all ROOT histogram bins, errors, entries, means, counters,
and legacy record fields. It does not write over the original ROOT outputs.
Then run the `analyze`, `display`, `parents`, and `plot` examples above to check
file output and figures in that environment.

The 19 ROOT-independent tests passed. The real-file
check and visual verification have not been run: system Python lacks PyROOT,
the local EIC container failed to start in the sandbox, and permission to run it
outside the sandbox was declined.
