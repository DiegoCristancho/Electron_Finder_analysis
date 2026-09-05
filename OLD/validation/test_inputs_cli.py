import json
from contextlib import nullcontext, redirect_stdout
import io
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from analysis import AnalysisConfig
from analysis_io import EventWriter, prepare_output, resolve_inputs
import run_analysis
from OLD.validation.fixtures import EventChain, selection_cases
from OLD.validation.reference import PROJECT
from OLD.validation.root_stub import Histogram, ROOT


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
        self.assertEqual(AnalysisConfig.load(distance_cut_mm=30).distance_cut_mm, 30)
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
        for script, command in (("run_analysis.py", []), ("run_analysis.py", ["analyze"]),
                                ("run_analysis.py", ["display"]), ("run_event_display.py", [])):
            result = subprocess.run([sys.executable, str(PROJECT / script), *command, "--help"], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("usage:", result.stdout)

    def test_cli_rejects_incompatible_inputs(self):
        result = subprocess.run([sys.executable, str(PROJECT / "run_analysis.py"), "analyze", "--input", "sample.root", "--file-list", "files.txt"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("not allowed", result.stderr)

    def test_direct_root_and_text_list_run_the_same_analysis_and_save_summary(self):
        """Exercise parsing, processing and output together through both input forms."""
        with TemporaryDirectory() as tmp:
            folder = Path(tmp)
            sample = folder / "sample.root"
            sample.touch()
            listing = folder / "files.txt"
            listing.write_text("# one file\nsample.root\n")
            summaries = []
            for value in (sample, listing):
                chain = EventChain(selection_cases())
                output = folder / value.suffix[1:]
                stdout = io.StringIO()
                with patch.dict(sys.modules, {"ROOT": ROOT}), \
                        patch.object(ROOT.gROOT, "GetVersion", return_value="test", create=True), \
                        patch("analysis_io.build_chain", return_value=chain) as build, \
                        patch("histograms.root_output", return_value=nullcontext()) as root_output, \
                        patch.object(Histogram, "Write", create=True), \
                        redirect_stdout(stdout):
                    run_analysis.main([str(value), "--no-plots", "--progress-every", "0",
                                       "--output", str(output)])
                self.assertEqual(build.call_args.args[0], [str(sample)])
                self.assertEqual(root_output.call_args.args[0], output / "analysis.root")
                self.assertEqual(chain.read_entries, list(range(16)))
                self.assertEqual(len((output / "events.jsonl").read_text().splitlines()), 16)
                summaries.append(json.loads((output / "summary.json").read_text()))
                self.assertIn((output / "summary.txt").read_text(), stdout.getvalue())
                self.assertEqual(json.loads((output / "run.json").read_text())["status"], "complete")
                self.assertEqual(json.loads((output / "config.json").read_text()), AnalysisConfig().to_dict())
            self.assertEqual(summaries[0], summaries[1])

    def test_new_cli_rejects_missing_or_conflicting_positional_inputs(self):
        for argv in ([], ["sample.root", "--file-list", "files.txt"],
                     ["sample.root", "--input", "second.root"]):
            with self.subTest(argv=argv), patch("sys.stderr", new=io.StringIO()), \
                    self.assertRaises(SystemExit) as error:
                run_analysis.main(argv)
            self.assertEqual(error.exception.code, 2)

    def test_failed_run_is_recorded_after_processing_error(self):
        with TemporaryDirectory() as tmp:
            folder = Path(tmp)
            sample, output = folder / "sample.root", folder / "failed"
            sample.touch()
            with patch("analysis_io.build_chain", return_value=EventChain(selection_cases())), \
                    patch("analysis_io.load_root", return_value=SimpleNamespace(gROOT=SimpleNamespace(GetVersion=lambda: "test"))), \
                    patch("run_analysis.ElectronAnalysis", side_effect=RuntimeError("deliberate test failure")), \
                    redirect_stdout(io.StringIO()), patch("sys.stderr", new=io.StringIO()), \
                    self.assertRaises(SystemExit) as error:
                run_analysis.main([str(sample), "--output", str(output)])
            self.assertEqual(error.exception.code, 2)
            manifest = json.loads((output / "run.json").read_text())
            self.assertEqual(manifest["status"], "failed")
            self.assertIn("deliberate test failure", manifest["error"])


if __name__ == "__main__":
    unittest.main()
