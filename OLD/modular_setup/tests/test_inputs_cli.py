import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

from electron_finder.config import AnalysisConfig
from electron_finder.io import EventWriter, prepare_output, resolve_inputs
from tests.reference import PROJECT


class InputsAndCliTests(unittest.TestCase):
    def test_file_list_resolves_beside_list_and_preserves_urls(self):
        with TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "data.root").touch()
            listing = folder / "files.txt"
            listing.write_text("# comment\n\ndata.root\nroot://host//data.root\n")
            self.assertEqual(resolve_inputs(file_list=listing), [str(folder / "data.root"), "root://host//data.root"])

    def test_duplicates_and_empty_inputs_are_errors(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            resolve_inputs(["root://host/file.root", "root://host/file.root"])
        with self.assertRaises(ValueError):
            resolve_inputs([])

    def test_output_does_not_overwrite_previous_run(self):
        with TemporaryDirectory() as tmp:
            output = prepare_output(Path(tmp) / "run")
            (output / "config.json").write_text("{}")
            with self.assertRaises(FileExistsError):
                prepare_output(output)

    def test_config_overrides_and_validation(self):
        baseline = PROJECT / "configs" / "baseline.json"
        self.assertEqual(AnalysisConfig.load(baseline, distance_cut_mm=30).distance_cut_mm, 30)
        for changes in ({"eop_min": 2}, {"distance_cut_mm": float("nan")}, {"cut_scan_bins": 2.5}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                AnalysisConfig(**changes)
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.json"
            path.write_text('{"distance_typo": 30}')
            with self.assertRaisesRegex(ValueError, "Unknown"):
                AnalysisConfig.load(path)

    def test_records_include_only_failures_in_failure_file(self):
        with TemporaryDirectory() as tmp:
            with EventWriter(tmp) as writer:
                writer.write({"event": 0, "reason": None})
                writer.write({"event": 1, "reason": "no_reco_match"})
            events = (Path(tmp) / "events.jsonl").read_text().splitlines()
            failures = (Path(tmp) / "failures.jsonl").read_text().splitlines()
            self.assertEqual(len(events), 2)
            self.assertEqual([json.loads(line)["event"] for line in failures], [1])

    def test_cli_help_does_not_require_root(self):
        for command in ([], ["analyze"], ["display"], ["parents"], ["plot"]):
            result = subprocess.run([sys.executable, str(PROJECT / "run_analysis.py"), *command, "--help"], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("usage:", result.stdout)

    def test_cli_rejects_incompatible_inputs(self):
        result = subprocess.run([sys.executable, str(PROJECT / "run_analysis.py"), "analyze", "--input", "sample.root", "--file-list", "files.txt"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("not allowed", result.stderr)


if __name__ == "__main__":
    unittest.main()
