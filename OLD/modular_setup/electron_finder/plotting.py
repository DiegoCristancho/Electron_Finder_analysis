"""Draw saved analysis histograms without rereading event files."""

from pathlib import Path
from .io import load_root, root_output


def plot_saved_results(results_dir, output_dir):
    """Use the settings saved with a completed run so E/p labels stay correct."""
    from .config import AnalysisConfig
    from .histograms import book_histograms
    results_dir, output_dir = Path(results_dir), Path(output_dir)
    config = AnalysisConfig.load(results_dir / "config.json")
    ROOT = load_root()
    source = ROOT.TFile.Open(str(results_dir / "analysis.root"), "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"Could not read {results_dir / 'analysis.root'}")
    histograms = {}
    try:
        names = list(book_histograms(config)) + ["h_frac_d3D_noEoP", "h_frac_d3D_withEoP"]
        for name in names:
            histogram = source.Get(name)
            if not histogram or not histogram.InheritsFrom("TH1"):
                raise RuntimeError(f"Missing histogram {name!r} in saved analysis")
            histogram.SetDirectory(0)
            histograms[name] = histogram
    finally:
        source.Close()
    canvases = plot_analysis(histograms, config, output_dir / "plots")
    with root_output(output_dir / "plots.root"):
        for canvas in canvases.values():
            canvas.Write()
    return {"figures": len(canvases), "source_results": str(results_dir.resolve())}


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
