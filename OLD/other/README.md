# Original analysis references

The organized workflow starts at `../../run_analysis.py`; see the project README
for environment setup and commands. These original scripts and existing ROOT
files are retained as references while the new workflow is validated on data.

| Original | New entry point |
| --- | --- |
| `Distances_Between_cluster_track.py` | `analyze` |
| `Trial.py` | Overlapping matching study; use `analyze` for new work |
| `Analysis.py` | `run_event_display.py` (or `run_analysis.py display`) |
| `Tracks_Vs_Clusters.py` | `run_event_display.py --tracks-only` |
| `checking_parent.py` | `parents` |

Run new commands from the project root. The local ROOT sample is already listed
in `../../inputs/local_example.txt`. The reference scripts are still standalone
scripts with relative paths and fixed output names; their functions should not
be imported into new work. Import from `electron_finder` instead.
