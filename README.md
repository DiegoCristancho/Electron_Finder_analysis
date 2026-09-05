# Electron finder: two runners and three helper files

All five files are beside your ROOT inputs. Run the commands inside your existing
EIC environment, where Python can import ROOT and the EDM dictionaries.

| File | What you do here |
| --- | --- |
| `run_analysis.py` | Launch the full-file matching analysis. |
| `run_event_display.py` | Launch a single-event track/cluster display. |
| `analysis_io.py` | Command arguments, ROOT paths or `.txt` lists, event reading, and run records. |
| `analysis.py` | **All physics, electron selection, matching, cuts, event inspection, and printed summaries.** |
| `histograms.py` | **Histogram definitions, plot styling, and ROOT output for both workflows.** |

The input helper is called `analysis_io.py` because Python already has a built-in
module named `io`. There is no analysis package or required configuration folder.

## Run an analysis

From this directory, inside your EIC environment:

```bash
# Process every event in one ROOT file.
python3 run_analysis.py sample.root

# Or process the ROOT files named in a text list.
python3 run_analysis.py files.txt

# Try the first 100 events of the local sample already in this directory.
python3 run_analysis.py local_example.txt --max-events 100
```

The local example list points to
`pythia8NCDIS_26_718x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_1.0008.eicrecon.edm4eic.root`
beside the runners. Remove `--max-events` to analyze the whole file.

The detailed **Trial.py-style summary** is printed at completion and saved as
`summary.txt`: counts, mean distances, position differences, parent information,
plateaus, passing/failing counts, matching fractions, and failure reasons.

## Display one event

```bash
python3 run_event_display.py sample.root --event 1
python3 run_event_display.py files.txt --event 1

# Your local sample:
python3 run_event_display.py local_example.txt --event 1
```

This prints the detailed projection/track/cluster information and makes the four
views: **3D tracks, XY tracks, XY tracks + ECal clusters, and RZ tracks + ECal
clusters**. It displays every projected track in the selected event.

Event indices start at zero: `--event 1` is the second entry. RZ uses Z horizontally
and cylindrical R = sqrt(x*x + y*y) vertically. Add `--tracks-only` for projection
files without clusters or truth associations; that saves the 3D and XY views.

## Change cuts or calculations

Edit `AnalysisConfig` near the top of **analysis.py**. Current defaults include:

```python
ecal_system: int = 101
distance_cut_mm: float = 87.0
eop_min: float = 0.7
eop_max: float = 1.3
```

For a temporary override:

```bash
python3 run_analysis.py local_example.txt --distance-cut 60 --eop-min 0.8 --eop-max 1.2
```

The same file contains `ElectronAnalysis.process_event()` (the sequence for one
event), the physics calculations, `find_first_scattered_electron_mc_index()`
(electron selection), `match_clusters()` (cluster selection), and `inspect_event()`
(single-event coordinates and associations). The summary formatting is there too.

In **histograms.py**, edit `book_histograms()` for bins/titles,
`plot_analysis()` for matching figures, or `draw_event_display()` for event views.
`save_analysis_results()` writes the matching ROOT file and summary outputs.

The matching definitions are preserved: the first status-1 MC electron with PDG 11
and first-parent status different from 2; a missing/invalid parent index passes the
original veto. Matching follows its reconstructed particle, first track and ECAL
projection point. It finds the nearest cluster in 3D both without an E/p cut and
among clusters passing the strict E/p window. Distance cuts are strict too.
The matching-fraction denominator is electrons with a reconstructed ECAL projection
point. The printed legacy failing counts still include missing ECAL points.

## Inputs and outputs

You can supply several ROOT paths after the script name. A `.txt` list must be
supplied by itself and contains one ROOT path/URL per line. Blank lines and lines
starting with `#` are ignored. Relative paths inside a list are resolved beside
that list; direct relative paths use the current working directory. `root://` URLs
are supported. Duplicate input paths are rejected.

Each invocation writes a new folder in `results/`. Use `--output results/my_run`
to choose a name; an existing nonempty folder is protected from overwriting.

| Full-file analysis | Single-event display |
| --- | --- |
| `analysis.root` | `event_display.root` |
| Eight figures in `plots/` | Four event figures in `plots/` |
| `events.jsonl`, `failures.jsonl` | `event.json`, `event.txt` |
| `summary.json`, `summary.txt` | `summary.json` |

Both save the effective settings in `config.json`, input paths in `inputs.txt`,
and completion/failure information in `run.json`. JSONL records are written during
processing, so all event records do not need to stay in memory.

`--no-plots` skips drawing the full-file figures while keeping the histograms and
summary. `--progress-every 0` hides progress messages but keeps the final summary.
`run_analysis.py --event 1` applies the matching study to one entry; use the
dedicated display runner for the geometric views.

The previous `analyze --input ...`, `analyze --file-list ...`, and
`display --input ...` command forms still work. Example lists formerly in `inputs/`
are now beside the runners, so use `local_example.txt` instead of
`inputs/local_example.txt`. An optional `--config path.json` can load saved settings;
ordinary runs use the defaults in `analysis.py`.

## Previous setup and validation

`OLD/modular_setup/` preserves the previous package, runners, configuration,
documentation, and tests. Its Word guide describes that archived layout. The extra
`parents` and `plot` subcommands belong to that archived runner. The active runners
handle the two workflows above. Earlier studies remain in `OLD/other/` and
`OLD/single_electron/`. Existing results are kept in `results/`.

The independent `4_E-pz_histograms.py` study is unchanged and is separate from this
five-file pipeline.

The current regression checks are in `OLD/validation/`, outside the everyday
workflow. From the project directory:

```bash
python3 -m unittest discover -s OLD/validation -v

# Optional real-file comparison, inside the EIC environment:
python3 -m OLD.validation.compare_root_sample --file-list local_example.txt --max-events 100
```

All 22 ROOT-independent checks pass. They cover the original histogram contents,
selection counts and event records, Trial.py report values, event-view coordinates,
direct ROOT/list commands, saved summaries, and failure reporting. Real ROOT
serialization/rendering has not been rerun for this simplification; this terminal's
Python does not have PyROOT.
