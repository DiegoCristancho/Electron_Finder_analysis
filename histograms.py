"""Histogram definitions, figure styling, and ROOT output for both runners.

Edit book_histograms for bins/titles, plot_analysis for full-file figures,
or draw_event_display for 3D, XY and RZ event views.
"""

from contextlib import contextmanager
from pathlib import Path

from analysis_io import load_root, write_json



# ========================================================================
# HISTOGRAM DEFINITIONS: names and binning preserved from the original
# ========================================================================

def book_histograms(config):
    ROOT = load_root()
    histograms = {}
    histograms['h_d3D_noEoP'] = ROOT.TH1F('h_d3D_noEoP', 'Nearest ECal cluster distance for scattered electron (no E/p cut);d_{3D} [mm];Events', 100, 0, 200)
    histograms['h_d3D_withEoP'] = ROOT.TH1F('h_d3D_withEoP', f'Nearest ECal cluster distance for scattered electron ({config.eop_min:g} < E/p < {config.eop_max:g});d_{{3D}} [mm];Events', 100, 0, 200)
    histograms['h_dr_noEoP'] = ROOT.TH1F('h_dr_noEoP', 'Cluster - projection #Delta r (no E/p cut);r_{cluster} - r_{proj} [mm];Events', 120, -150, 150)
    histograms['h_dtheta_noEoP'] = ROOT.TH1F('h_dtheta_noEoP', 'Cluster - projection #Delta#theta (no E/p cut);#theta_{cluster} - #theta_{proj} [rad];Events', 120, -0.2, 0.2)
    histograms['h_deta_noEoP'] = ROOT.TH1F('h_deta_noEoP', 'Cluster - projection #Delta#eta (no E/p cut);#eta_{cluster} - #eta_{proj};Events', 120, -0.5, 0.5)
    histograms['h_dr_withEoP'] = ROOT.TH1F('h_dr_withEoP', f'Cluster - projection #Delta r ({config.eop_min:g} < E/p < {config.eop_max:g});r_{{cluster}} - r_{{proj}} [mm];Events', 120, -150, 150)
    histograms['h_dtheta_withEoP'] = ROOT.TH1F('h_dtheta_withEoP', f'Cluster - projection #Delta#theta ({config.eop_min:g} < E/p < {config.eop_max:g});#theta_{{cluster}} - #theta_{{proj}} [rad];Events', 120, -0.2, 0.2)
    histograms['h_deta_withEoP'] = ROOT.TH1F('h_deta_withEoP', f'Cluster - projection #Delta#eta ({config.eop_min:g} < E/p < {config.eop_max:g});#eta_{{cluster}} - #eta_{{proj}};Events', 120, -0.5, 0.5)
    histograms['h_eop_all'] = ROOT.TH1F('h_eop_all', 'E/p of all ECal clusters considered for scattered electron;E/p;Clusters', 120, 0, 3)
    histograms['h_eop_selected'] = ROOT.TH1F('h_eop_selected', f'E/p of selected clusters passing {config.eop_min:g} < E/p < {config.eop_max:g};E/p;Clusters', 120, 0, 3)
    histograms['h_total'] = ROOT.TH1F('h_total', 'Total scattered electrons with ECAL projection point;distance cut [mm];Events', config.cut_scan_bins, 0, config.cut_scan_max_mm)
    histograms['h_pass_d3D_noEoP'] = ROOT.TH1F('h_pass_d3D_noEoP', 'Scattered electrons with nearest cluster within d_{3D} cut (no E/p);distance cut [mm];Events', config.cut_scan_bins, 0, config.cut_scan_max_mm)
    histograms['h_pass_d3D_withEoP'] = ROOT.TH1F('h_pass_d3D_withEoP', f'Scattered electrons with nearest cluster within d_{{3D}} cut ({config.eop_min:g} < E/p < {config.eop_max:g});distance cut [mm];Events', config.cut_scan_bins, 0, config.cut_scan_max_mm)
    histograms['h_fail_recoE'] = ROOT.TH1F('h_fail_recoE', f'Reco energy of scattered electrons failing d_{{3D}} < {config.distance_cut_mm:g} mm;Reco electron energy [GeV];Events', 100, 0, 30)
    histograms['h_fail_mcE'] = ROOT.TH1F('h_fail_mcE', f'MC energy of scattered electrons failing d_{{3D}} < {config.distance_cut_mm:g} mm;MC electron energy [GeV];Events', 100, 0, 30)
    histograms['h_parent_status_all'] = ROOT.TH1F('h_parent_status_all', 'Parent generator status of all electrons with status=1 and PDG=11;Parent generator status;Electrons', 60, -0.5, 59.5)
    histograms['h_parent_status_selected'] = ROOT.TH1F('h_parent_status_selected', 'Parent generator status of selected scattered electrons (parent status != 2);Parent generator status;Events', 60, -0.5, 59.5)
    histograms['h_parent_pdg_all'] = ROOT.TH1F('h_parent_pdg_all', 'Parent PDG of all electrons with status=1 and PDG=11;Parent PDG;Electrons', 200, -100, 100)
    for histogram in histograms.values():
        histogram.SetDirectory(0)
    return histograms


# ========================================================================
# FULL-FILE FIGURES
# ========================================================================

def plot_analysis(histograms, config, directory):
    ROOT = load_root()
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    c1 = ROOT.TCanvas('c1', 'd3D without E/p', 900, 700)
    histograms['h_d3D_noEoP'].SetLineWidth(2)
    histograms['h_d3D_noEoP'].Draw('HIST')
    c2 = ROOT.TCanvas('c2', 'd3D with E/p', 900, 700)
    histograms['h_d3D_withEoP'].SetLineWidth(2)
    histograms['h_d3D_withEoP'].SetLineColor(ROOT.kBlue + 1)
    histograms['h_d3D_withEoP'].Draw('HIST')
    c3 = ROOT.TCanvas('c3', 'fraction vs distance cut', 900, 700)
    histograms['h_frac_d3D_noEoP'].SetLineWidth(2)
    histograms['h_frac_d3D_noEoP'].SetLineColor(ROOT.kBlack)
    histograms['h_frac_d3D_noEoP'].SetMinimum(0.0)
    histograms['h_frac_d3D_noEoP'].SetMaximum(1.05)
    histograms['h_frac_d3D_noEoP'].Draw('HIST')
    histograms['h_frac_d3D_withEoP'].SetLineWidth(2)
    histograms['h_frac_d3D_withEoP'].SetLineColor(ROOT.kRed)
    histograms['h_frac_d3D_withEoP'].Draw('HIST SAME')
    leg = ROOT.TLegend(0.52, 0.73, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(histograms['h_frac_d3D_noEoP'], 'No E/p cut', 'l')
    leg.AddEntry(histograms['h_frac_d3D_withEoP'], f'{config.eop_min:g} < E/p < {config.eop_max:g}', 'l')
    leg.Draw()
    c4 = ROOT.TCanvas('c4', 'failing scattered electron energies', 1200, 500)
    c4.Divide(2, 1)
    c4.cd(1)
    histograms['h_fail_recoE'].SetLineWidth(2)
    histograms['h_fail_recoE'].Draw('HIST')
    c4.cd(2)
    histograms['h_fail_mcE'].SetLineWidth(2)
    histograms['h_fail_mcE'].Draw('HIST')
    c5 = ROOT.TCanvas('c5', 'parent status histograms', 1200, 500)
    c5.Divide(2, 1)
    c5.cd(1)
    histograms['h_parent_status_all'].SetLineWidth(2)
    histograms['h_parent_status_all'].Draw('HIST')
    c5.cd(2)
    histograms['h_parent_status_selected'].SetLineWidth(2)
    histograms['h_parent_status_selected'].SetLineColor(ROOT.kBlue + 1)
    histograms['h_parent_status_selected'].Draw('HIST')
    c6 = ROOT.TCanvas('c6', 'parent PDG all electrons', 900, 700)
    histograms['h_parent_pdg_all'].SetLineWidth(2)
    histograms['h_parent_pdg_all'].Draw('HIST')
    c7 = ROOT.TCanvas('c7', 'E over p', 1200, 500)
    c7.Divide(2, 1)
    c7.cd(1)
    histograms['h_eop_all'].SetLineWidth(2)
    histograms['h_eop_all'].Draw('HIST')
    c7.cd(2)
    histograms['h_eop_selected'].SetLineWidth(2)
    histograms['h_eop_selected'].SetLineColor(ROOT.kBlue + 1)
    histograms['h_eop_selected'].Draw('HIST')
    c8 = ROOT.TCanvas('c8', 'position based cluster projection differences', 1500, 900)
    c8.Divide(3, 2)
    c8.cd(1)
    histograms['h_dr_noEoP'].SetLineWidth(2)
    histograms['h_dr_noEoP'].Draw('HIST')
    c8.cd(2)
    histograms['h_dtheta_noEoP'].SetLineWidth(2)
    histograms['h_dtheta_noEoP'].Draw('HIST')
    c8.cd(3)
    histograms['h_deta_noEoP'].SetLineWidth(2)
    histograms['h_deta_noEoP'].Draw('HIST')
    c8.cd(4)
    histograms['h_dr_withEoP'].SetLineWidth(2)
    histograms['h_dr_withEoP'].SetLineColor(ROOT.kBlue + 1)
    histograms['h_dr_withEoP'].Draw('HIST')
    c8.cd(5)
    histograms['h_dtheta_withEoP'].SetLineWidth(2)
    histograms['h_dtheta_withEoP'].SetLineColor(ROOT.kBlue + 1)
    histograms['h_dtheta_withEoP'].Draw('HIST')
    c8.cd(6)
    histograms['h_deta_withEoP'].SetLineWidth(2)
    histograms['h_deta_withEoP'].SetLineColor(ROOT.kBlue + 1)
    histograms['h_deta_withEoP'].Draw('HIST')
    canvases = {
        "distance_no_eop": c1, "distance_with_eop": c2,
        "fraction_vs_distance": c3, "failed_energies": c4,
        "parent_status": c5, "parent_pdg": c6,
        "e_over_p": c7, "position_differences": c8,
    }
    for name, canvas in canvases.items():
        canvas._keepalive = [leg, *histograms.values()]
        canvas.SaveAs(str(directory / (name + ".png")))
    return canvases


# ========================================================================
# SINGLE-EVENT FIGURES: 3D, XY, XY + clusters, RZ + clusters
# ========================================================================

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


# ========================================================================
# ROOT FILES AND SAVED ANALYSIS RESULTS
# ========================================================================

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


def save_analysis_results(histograms, config, summary, report, directory, make_plots=True):
    """Save ROOT histograms/canvases, figures, and the detailed summary."""
    directory = Path(directory)
    canvases = plot_analysis(histograms, config, directory / "plots") if make_plots else {}
    with root_output(directory / "analysis.root"):
        for histogram in histograms.values():
            histogram.Write()
        for canvas in canvases.values():
            canvas.Write()
    write_json(directory / "summary.json", summary)
    (directory / "summary.txt").write_text(report)
