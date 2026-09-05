"""Histogram definitions; baseline names and binning match the original script."""

from .io import load_root


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
