#!/usr/bin/env python3
"""Inspect one event and draw its tracks/projection points and ECal clusters.

Example (inside the EIC environment):
    python3 run_event_display.py --input sample.root --event 1

This is the single-event counterpart to run_analysis.py analyze. Both runners
use the same input handling, configuration, and reusable analysis package.
"""

from run_analysis import main as run_main


def main(argv=None):
    run_main(argv, default_command="display")


if __name__ == "__main__":
    main()
