"""Single-event inspection and first-status-1-electron parent diagnostics."""

from collections import Counter
import math
from pathlib import Path

from .io import EventWriter, event_source, iter_events, load_root, root_output, write_json
from .matching import find_mc_pdg_for_reco_index, find_reco_particle_for_track


PROJECTION_BRANCHES = (
    "CalorimeterTrackProjections", "_CalorimeterTrackProjections_track",
    "_CalorimeterTrackProjections_points",
)
PARENT_BRANCHES = ("MCParticles", "_MCParticles_parents")


def display_branches(config, tracks_only=False):
    if tracks_only:
        return PROJECTION_BRANCHES
    return PROJECTION_BRANCHES + (
        config.cluster_collection, "ReconstructedParticles", "_ReconstructedParticles_tracks",
        "_ReconstructedParticleAssociations_rec", "_ReconstructedParticleAssociations_sim", "MCParticles",
    )


def inspect_event(event, index, config, tracks_only=False):
    """Copy positions and associations out of ROOT's current event view."""
    tracks = []
    points = event._CalorimeterTrackProjections_points
    for projection_index, projection in enumerate(event.CalorimeterTrackProjections):
        track = event._CalorimeterTrackProjections_track[projection_index]
        reco_index = truth_pdg = None
        if not tracks_only:
            reco_index = find_reco_particle_for_track(
                track, event.ReconstructedParticles, event._ReconstructedParticles_tracks,
            )
            truth_pdg = find_mc_pdg_for_reco_index(
                reco_index, event._ReconstructedParticleAssociations_rec,
                event._ReconstructedParticleAssociations_sim, event.MCParticles,
            )
        positions = []
        for point_index in range(int(projection.points_begin), int(projection.points_end)):
            point = points[point_index]
            x, y, z = float(point.position.x), float(point.position.y), float(point.position.z)
            positions.append({
                "point_index": point_index, "x": x, "y": y, "z": z,
                "rho": math.sqrt(x*x + y*y), "system": int(point.system), "surface": int(point.surface),
            })
        tracks.append({
            "projection_index": projection_index, "track_index": int(track.index),
            "collection_id": int(track.collectionID), "reco_index": reco_index,
            "truth_pdg": truth_pdg, "points_begin": int(projection.points_begin),
            "points_end": int(projection.points_end), "points": positions,
        })
    clusters = []
    if not tracks_only:
        for cluster_index, cluster in enumerate(getattr(event, config.cluster_collection)):
            x, y, z = float(cluster.position.x), float(cluster.position.y), float(cluster.position.z)
            clusters.append({"index": cluster_index, "x": x, "y": y, "z": z,
                             "rho": math.sqrt(x*x + y*y), "energy": float(cluster.energy)})
    counts = {
        "projections": len(event.CalorimeterTrackProjections),
        "track_refs": len(event._CalorimeterTrackProjections_track),
        "projection_points": len(points),
    }
    if not tracks_only:
        counts.update(
            clusters=len(clusters), reconstructed_particles=len(event.ReconstructedParticles),
            reconstructed_particle_tracks=len(event._ReconstructedParticles_tracks),
            reco_mc_associations=len(event._ReconstructedParticleAssociations_rec),
            mc_particles=len(event.MCParticles),
        )
    return {"event": index, "tracks_only": tracks_only, "tracks": tracks,
            "cluster_collection": config.cluster_collection, "clusters": clusters, "counts": counts}


def draw_event_display(record, directory):
    """Draw the original 3D, XY, XY+clusters and cylindrical RZ views."""
    root = load_root()
    root.gStyle.SetOptStat(0)
    directory = Path(directory)
    plots = directory / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    positions = [point for track in record["tracks"] for point in track["points"]] + record["clusters"]
    ranges = {}
    for axis in ("x", "y", "z", "rho"):
        values = [point[axis] for point in positions] or [0.0]
        lower = min(values) - 20.0
        ranges[axis] = (max(0.0, lower) if axis == "rho" else lower, max(values) + 20.0)
    colors = [root.kRed + 1, root.kBlue + 1, root.kGreen + 2, root.kMagenta + 1,
              root.kOrange + 7, root.kCyan + 1, root.kViolet + 1, root.kBlack]
    views = ["3d", "xy"] if record["tracks_only"] else ["3d", "xy", "xy_clusters", "rz_clusters"]
    keepalive = []
    canvases = []
    with root_output(directory / "event_display.root"):
        for view in views:
            canvas = root.TCanvas(f"event_{view}", f"Event {record['event']} {view}", 1100, 900)
            canvases.append(canvas)
            legend = root.TLegend(0.60, 0.60, 0.93, 0.90)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.SetTextSize(0.025)
            is3d = view == "3d"
            axis_x, axis_y = ("z", "rho") if view == "rz_clusters" else ("x", "y")
            title = f"Event {record['event']} track projection points"
            if view.endswith("clusters"):
                title += " and ECal clusters"
            if is3d:
                frame = root.TH3F(f"frame_{view}", title + ";x [mm];y [mm];z [mm]",
                                 10, *ranges["x"], 10, *ranges["y"], 10, *ranges["z"])
            else:
                ylabel = "r = #sqrt{x^{2}+y^{2}}" if axis_y == "rho" else axis_y
                frame = root.TH2F(f"frame_{view}", title + f";{axis_x} [mm];{ylabel} [mm]",
                                 200, *ranges[axis_x], 200, *ranges[axis_y])
            frame.SetDirectory(0)
            frame.SetStats(0)
            frame.Draw()
            keepalive.extend([frame, legend])
            for i, track in enumerate(record["tracks"]):
                points = track["points"]
                if not points:
                    continue
                color = colors[i % len(colors)]
                if is3d:
                    graph = root.TPolyLine3D(len(points))
                    markers = root.TPolyMarker3D(len(points))
                    for j, point in enumerate(points):
                        graph.SetPoint(j, point["x"], point["y"], point["z"])
                        markers.SetPoint(j, point["x"], point["y"], point["z"])
                    graph.SetLineColor(color)
                    graph.SetLineWidth(2)
                    graph.Draw("SAME")
                    markers.SetMarkerColor(color)
                    markers.SetMarkerStyle(20)
                    markers.SetMarkerSize(1.0)
                    markers.Draw("SAME")
                    keepalive.extend([graph, markers])
                    entry = markers
                else:
                    graph = root.TGraph(len(points))
                    for j, point in enumerate(points):
                        graph.SetPoint(j, point[axis_x], point[axis_y])
                    graph.SetLineColor(color)
                    graph.SetMarkerColor(color)
                    graph.SetLineWidth(2)
                    graph.SetMarkerStyle(20)
                    graph.SetMarkerSize(1.1)
                    graph.Draw("LP SAME")
                    keepalive.append(graph)
                    entry = graph
                pdg = track["truth_pdg"] if track["truth_pdg"] is not None else "?"
                legend.AddEntry(entry, f"track {track['track_index']} (PDG {pdg})", "p" if is3d else "lp")
            if view.endswith("clusters") and record["clusters"]:
                graph = root.TGraph(len(record["clusters"]))
                for j, point in enumerate(record["clusters"]):
                    graph.SetPoint(j, point[axis_x], point[axis_y])
                graph.SetMarkerStyle(20)
                graph.SetMarkerColor(root.kBlack)
                graph.SetMarkerSize(1.2)
                graph.Draw("P SAME")
                legend.AddEntry(graph, record["cluster_collection"], "p")
                keepalive.append(graph)
            legend.Draw()
            canvas.SaveAs(str(plots / f"event_{view}.png"))
            canvas.Write()
    # Keep ROOT-owned drawing primitives alive until canvases have been written.
    for canvas in canvases:
        canvas.Close()
    return {"root_file": str(directory / "event_display.root"),
            "plots": [str(plots / f"event_{view}.png") for view in views]}


def run_parent_check(chain, directory, max_events=None, make_plots=True):
    """Check the first PDG=11/status=1 electron, without the parent-status veto."""
    root = load_root()
    status_hist = root.TH1F("h_parent_status", "Parent status of first status-1 electron;Parent generator status;Events", 60, -0.5, 59.5)
    pdg_hist = root.TH1F("h_parent_pdg", "Parent PDG of first status-1 electron;Parent PDG;Events", 100, -50, 50)
    for histogram in (status_hist, pdg_hist):
        histogram.SetDirectory(0)
    counts = Counter({key: 0 for key in ("events_processed", "with_electron", "with_parent", "no_parent", "bad_parent_index", "parent_status_23")})
    directory = Path(directory)
    with EventWriter(directory) as writer:
        for index, event in iter_events(chain, max_events):
            counts["events_processed"] += 1
            mc, refs = event.MCParticles, event._MCParticles_parents
            electron_index = next((i for i, part in enumerate(mc)
                                   if int(part.PDG) == 11 and int(part.generatorStatus) == 1), None)
            record = dict(event_source(chain, index), electron_index=electron_index,
                          parent_index=None, parent_status=None, parent_pdg=None, reason=None)
            if electron_index is None:
                record["reason"] = "no_status1_electron"
            else:
                counts["with_electron"] += 1
                particle = mc[electron_index]
                begin, end = int(particle.parents_begin), int(particle.parents_end)
                if end <= begin:
                    counts["no_parent"] += 1
                    record["reason"] = "no_parent"
                else:
                    parent_index = int(refs[begin].index)
                    record["parent_index"] = parent_index
                    if not 0 <= parent_index < len(mc):
                        counts["bad_parent_index"] += 1
                        record["reason"] = "bad_parent_index"
                    else:
                        parent = mc[parent_index]
                        status, pdg = int(parent.generatorStatus), int(parent.PDG)
                        counts["with_parent"] += 1
                        counts["parent_status_23"] += status == 23
                        record.update(parent_status=status, parent_pdg=pdg)
                        status_hist.Fill(status)
                        pdg_hist.Fill(pdg)
            writer.write(record)
    canvases = []
    with root_output(directory / "parents.root"):
        for histogram in (status_hist, pdg_hist):
            histogram.Write()
            if make_plots:
                plots = directory / "plots"
                plots.mkdir(exist_ok=True)
                canvas = root.TCanvas(histogram.GetName() + "_canvas", histogram.GetTitle(), 900, 700)
                histogram.SetLineWidth(2)
                histogram.Draw("HIST")
                canvas.SaveAs(str(plots / (histogram.GetName() + ".png")))
                canvas.Write()
                canvases.append(canvas)
    for canvas in canvases:
        canvas.Close()
    summary = {"cutflow": dict(counts), "selection": "First MC electron with PDG=11 and generatorStatus=1; first parent only",
               "fraction_parent_status_23": counts["parent_status_23"] / counts["with_parent"] if counts["with_parent"] else 0.0}
    write_json(directory / "summary.json", summary)
    return summary
