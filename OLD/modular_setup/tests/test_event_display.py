"""Check event-view coordinates and dispatch without requiring ROOT graphics."""

from contextlib import nullcontext, redirect_stdout
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace as NS
import unittest
from unittest.mock import MagicMock, patch

from electron_finder.config import AnalysisConfig
from electron_finder.diagnostics import draw_event_display, inspect_event
from electron_finder.reporting import format_event_details
import run_event_display
from tests.fixtures import EventChain, base_event, vector


def event_with_two_tracks():
    event = base_event()
    event.CalorimeterTrackProjections = [NS(points_begin=0, points_end=2), NS(points_begin=2, points_end=3)]
    event._CalorimeterTrackProjections_track.append(NS(collectionID=7, index=9))
    event._CalorimeterTrackProjections_points = [
        NS(position=vector(3, 4, 12), system=100, surface=1),
        NS(position=vector(6, 8, 24), system=101, surface=2),
        NS(position=vector(5, 12, 0), system=101, surface=2),
    ]
    event.EcalBarrelClusters = [NS(position=vector(-8, 15, 30), energy=10),
                                NS(position=vector(0, -7, -20), energy=5)]
    return event


class GraphicsRecorder:
    """Record ROOT graph coordinates and canvas saves, without rendering images."""

    def __init__(self):
        self.root = MagicMock()
        for index, color in enumerate(("kRed", "kBlue", "kGreen", "kMagenta", "kOrange", "kCyan", "kViolet", "kBlack")):
            setattr(self.root, color, index)
        self.canvases = {}
        self.graphs = {}
        self.root.TCanvas.side_effect = self.canvas
        self.root.TGraph.side_effect = self.graph

    def canvas(self, name, *_):
        self.current = name
        self.canvases[name] = MagicMock()
        self.graphs[name] = []
        return self.canvases[name]

    def graph(self, *_):
        graph = MagicMock()
        self.graphs[self.current].append(graph)
        return graph

    def coordinates(self, view):
        return [[call.args for call in graph.SetPoint.call_args_list] for graph in self.graphs[view]]


class EventDisplayTests(unittest.TestCase):
    def test_xy_and_rz_include_every_track_and_cluster(self):
        record = inspect_event(event_with_two_tracks(), 1, AnalysisConfig())
        graphics = GraphicsRecorder()
        with TemporaryDirectory() as directory, \
                patch("electron_finder.diagnostics.load_root", return_value=graphics.root), \
                patch("electron_finder.diagnostics.root_output", return_value=nullcontext()):
            files = draw_event_display(record, directory)
        self.assertEqual(graphics.coordinates("event_xy_clusters"), [
            [(0, 3, 4), (1, 6, 8)], [(0, 5, 12)], [(0, -8, 15), (1, 0, -7)],
        ])
        # R is cylindrical sqrt(x*x+y*y); nonzero z must not enter the radius.
        self.assertEqual(graphics.coordinates("event_rz_clusters"), [
            [(0, 12, 5), (1, 24, 10)], [(0, 0, 13)], [(0, 30, 17), (1, -20, 7)],
        ])
        self.assertEqual(len(graphics.coordinates("event_xy")), 2)
        self.assertEqual(len(files["plots"]), 4)
        for name, canvas in graphics.canvases.items():
            self.assertEqual(Path(canvas.SaveAs.call_args.args[0]).name, name + ".png")
            canvas.Write.assert_called_once()
        details = format_event_details(record)
        self.assertIn("Projection 1", details)
        self.assertIn("track.index        = 9", details)
        self.assertIn("x=3.000, y=4.000, z=12.000, r=5.000", details)
        self.assertIn("cluster 0: x=-8.000, y=15.000, z=30.000, r=17.000", details)

    def test_projection_only_mode_does_not_require_mc_or_clusters(self):
        full = event_with_two_tracks()
        event = NS(**{name: getattr(full, name) for name in (
            "CalorimeterTrackProjections", "_CalorimeterTrackProjections_track", "_CalorimeterTrackProjections_points",
        )})
        record = inspect_event(event, 1, AnalysisConfig(), tracks_only=True)
        graphics = GraphicsRecorder()
        with TemporaryDirectory() as directory, \
                patch("electron_finder.diagnostics.load_root", return_value=graphics.root), \
                patch("electron_finder.diagnostics.root_output", return_value=nullcontext()):
            files = draw_event_display(record, directory)
        self.assertEqual(set(graphics.canvases), {"event_3d", "event_xy"})
        self.assertEqual(len(files["plots"]), 2)
        self.assertNotIn("truth PDG", format_event_details(record))

    def test_dedicated_runner_selects_display_and_saves_event_details(self):
        chain = EventChain([base_event(), event_with_two_tracks()])
        stdout = io.StringIO()
        with TemporaryDirectory() as directory:
            input_path, output = Path(directory) / "sample.root", Path(directory) / "event1"
            input_path.touch()
            files = {"root_file": str(output / "event_display.root"),
                     "plots": [str(output / "plots" / "event_xy_clusters.png"),
                               str(output / "plots" / "event_rz_clusters.png")]}
            with patch("run_analysis.build_chain", return_value=chain), \
                    patch("run_analysis.load_root", return_value=NS(gROOT=NS(GetVersion=lambda: "test"))), \
                    patch("electron_finder.diagnostics.draw_event_display", return_value=files) as draw, \
                    patch("electron_finder.analysis.ElectronAnalysis", side_effect=AssertionError("Wrong workflow")), \
                    redirect_stdout(stdout):
                run_event_display.main(["--input", str(input_path), "--event", "1", "--output", str(output)])
            record = json.loads((output / "event.json").read_text())
            manifest = json.loads((output / "run.json").read_text())
            self.assertEqual(record["event"], 1)
            self.assertEqual(record["counts"]["projections"], 2)
            self.assertEqual(manifest["command"], "display")
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(chain.read_entries, [1])
            self.assertIn((output / "event.txt").read_text(), stdout.getvalue())
            self.assertIn("Saved event views:", stdout.getvalue())
            draw.assert_called_once()


if __name__ == "__main__":
    unittest.main()
