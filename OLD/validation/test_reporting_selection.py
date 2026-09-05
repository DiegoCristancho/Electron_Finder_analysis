import ast
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from analysis import ElectronAnalysis
from analysis_io import iter_events
from analysis import format_analysis_summary
from analysis_io import make_parser
from OLD.validation.fixtures import EventChain, selection_cases
from OLD.validation.reference import PROJECT, assigns, run_reference
from OLD.validation.root_stub import ROOT


class ReportingAndSelectionTests(unittest.TestCase):
    def test_report_preserves_every_trial_summary_value(self):
        events = selection_cases()
        reference = run_reference(events, ROOT)
        trial = PROJECT / "OLD" / "other" / "Trial.py"
        tree = ast.parse(trial.read_text())
        start = next(i for i, node in enumerate(tree.body)
                     if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
                     and node.value.args and isinstance(node.value.args[0], ast.Constant)
                     and node.value.args[0].value == "\n========== Summary ==========")
        end = next(i for i, node in enumerate(tree.body) if i > start and assigns(node, "c1"))
        original_stdout = io.StringIO()
        with redirect_stdout(original_stdout):
            exec(compile(ast.Module(body=tree.body[start:end], type_ignores=[]), str(trial), "exec"), reference)
        with patch.dict(sys.modules, {"ROOT": ROOT}):
            analysis = ElectronAnalysis()
            for index, event in enumerate(events):
                analysis.process_event(event, index)
            summary = analysis.summary()
            new_stdout = io.StringIO()
            with redirect_stdout(new_stdout):
                analysis.print_summary()

        def report_values(report):
            return {" ".join(label.split()): value.strip()
                    for line in report.splitlines() if ":" in line
                    for label, value in [line.split(":", 1)]}

        old_values, new_values = report_values(original_stdout.getvalue()), report_values(new_stdout.getvalue())
        self.assertGreater(len(old_values), 20)
        for label, value in old_values.items():
            self.assertEqual(value, new_values[label], label)
        # Saved JSON data can reproduce the exact same terminal summary.
        self.assertEqual(new_stdout.getvalue(), format_analysis_summary(json.loads(json.dumps(summary))))

    def test_single_event_run_loads_only_requested_entry_and_keeps_its_index(self):
        chain = EventChain(selection_cases())
        with TemporaryDirectory() as directory, patch.dict(sys.modules, {"ROOT": ROOT}):
            analysis = ElectronAnalysis()
            summary = analysis.run(chain, directory, event_index=1)
            records = [json.loads(line) for line in (Path(directory) / "events.jsonl").read_text().splitlines()]
        self.assertEqual(chain.read_entries, [1])
        self.assertEqual(summary["cutflow"]["events_processed"], 1)
        self.assertEqual(summary["cutflow"]["selected_electrons"], 0)
        self.assertEqual(len(records), 1)
        self.assertEqual((records[0]["event"], records[0]["source_entry"]), (1, 1))
        self.assertEqual(records[0]["reason"], "no_scattered_electron")

    def test_invalid_event_selection_fails_before_reading(self):
        for selection in ({"event_index": -1}, {"event_index": 16}, {"event_index": 1, "max_events": 3}):
            chain = EventChain(selection_cases())
            with self.subTest(selection=selection), self.assertRaises(ValueError):
                list(iter_events(chain, **selection))
            self.assertEqual(chain.read_entries, [])

    def test_cli_supports_event_one_and_rejects_max_events_together(self):
        parser = make_parser()
        args = parser.parse_args(["--input", "sample.root", "--event", "1"])
        self.assertEqual(args.event, 1)
        self.assertIsNone(args.max_events)
        with patch("sys.stderr", new=io.StringIO()), self.assertRaises(SystemExit) as error:
            parser.parse_args(["--input", "sample.root", "--event", "1", "--max-events", "10"])
        self.assertEqual(error.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
