import ROOT
import math

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(1110)

# ---------------------------------
# Input / output
# ---------------------------------

# Set this to True if you want to use the txt file list again.
# Set this to False if you want to use only one ROOT file.
USE_LISTFILE = False

# Old option: read many ROOT files from this txt list
LISTFILE = "list_18x275_Q2_100_26_02.txt"

# New option: read one ROOT file directly
INPUTROOT = "pythia8NCDIS_26_718x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_1.0008.eicrecon.edm4eic.root"

EVENT_TREE = "events"
OUTROOT = "scattered_electron_cluster_fraction_with_EoP.root"

ECAL_SYSTEM = 101
CHOSEN_CUT_MM = 87.0

# E/p window for electron-like clusters
EOP_MIN = 0.7
EOP_MAX = 1.3

# ---------------------------------
# Helpers
# ---------------------------------
def particle_energy(obj):
    px = float(obj.momentum.x)
    py = float(obj.momentum.y)
    pz = float(obj.momentum.z)
    m  = float(obj.mass) if hasattr(obj, "mass") else 0.0
    return math.sqrt(px * px + py * py + pz * pz + m * m)

def particle_p(obj):
    px = float(obj.momentum.x)
    py = float(obj.momentum.y)
    pz = float(obj.momentum.z)
    return math.sqrt(px * px + py * py + pz * pz)

def get_first_parent_status(part, mc_particles, parent_refs):
    pb = int(part.parents_begin)
    pe = int(part.parents_end)

    if pe <= pb:
        return None

    pref = parent_refs[pb]
    parent_idx = int(pref.index)

    if parent_idx < 0 or parent_idx >= len(mc_particles):
        return None

    return int(mc_particles[parent_idx].generatorStatus)

def get_first_parent_pdg(part, mc_particles, parent_refs):
    pb = int(part.parents_begin)
    pe = int(part.parents_end)

    if pe <= pb:
        return None

    pref = parent_refs[pb]
    parent_idx = int(pref.index)

    if parent_idx < 0 or parent_idx >= len(mc_particles):
        return None

    return int(mc_particles[parent_idx].PDG)

def find_first_scattered_electron_mc_index(mc_particles, parent_refs):
    """
    Scattered electron definition:
      generatorStatus == 1
      PDG == 11
      parent_status != 2
    """
    for i, part in enumerate(mc_particles):
        if int(part.generatorStatus) != 1:
            continue
        if int(part.PDG) != 11:
            continue

        parent_status = get_first_parent_status(part, mc_particles, parent_refs)
        if parent_status != 2:
            return i

    return None

def find_reco_index_for_mc_index(mc_index, assoc_rec, assoc_sim):
    if mc_index is None:
        return None

    for ia in range(len(assoc_sim)):
        sim_ref = assoc_sim[ia]
        if int(sim_ref.index) == int(mc_index):
            return int(assoc_rec[ia].index)

    return None

def find_track_ref_for_reco_index(reco_index, reco_particles, reco_particle_tracks):
    if reco_index is None:
        return None

    if reco_index < 0 or reco_index >= len(reco_particles):
        return None

    reco = reco_particles[reco_index]
    b = int(reco.tracks_begin)
    e = int(reco.tracks_end)

    if e <= b:
        return None

    return reco_particle_tracks[b]

def same_track_ref(a, b):
    return (int(a.collectionID) == int(b.collectionID)) and (int(a.index) == int(b.index))

# ---------------------------------
# Build chain
# ---------------------------------
chain = ROOT.TChain(EVENT_TREE)

if USE_LISTFILE:
    with open(LISTFILE) as f:
        for line in f:
            fn = line.strip()
            if fn:
                chain.Add(fn)
    print(f"✅ Added {chain.GetNtrees()} files from txt list, total {chain.GetEntries()} events")
else:
    chain.Add(INPUTROOT)
    print(f"✅ Added {chain.GetNtrees()} file, total {chain.GetEntries()} events")

# ---------------------------------
# Histograms
# ---------------------------------
# d3D without E/p
h_d3D_noEoP = ROOT.TH1F(
    "h_d3D_noEoP",
    "Nearest ECal cluster distance for scattered electron (no E/p cut);d_{3D} [mm];Events",
    100, 0, 200
)

# d3D with E/p
h_d3D_withEoP = ROOT.TH1F(
    "h_d3D_withEoP",
    f"Nearest ECal cluster distance for scattered electron ({EOP_MIN:.1f} < E/p < {EOP_MAX:.1f});d_{{3D}} [mm];Events",
    100, 0, 200
)

# optional E/p diagnostic
h_eop_all = ROOT.TH1F(
    "h_eop_all",
    "E/p of all ECal clusters considered for scattered electron;E/p;Clusters",
    120, 0, 3
)

h_eop_selected = ROOT.TH1F(
    "h_eop_selected",
    f"E/p of selected clusters passing {EOP_MIN:.1f} < E/p < {EOP_MAX:.1f};E/p;Clusters",
    120, 0, 3
)

NBINS_CUT = 100
DCUT_MAX = 200.0

# denominator
h_total = ROOT.TH1F(
    "h_total",
    "Total scattered electrons with ECAL projection point;distance cut [mm];Events",
    NBINS_CUT, 0, DCUT_MAX
)

# pass counts without E/p
h_pass_d3D_noEoP = ROOT.TH1F(
    "h_pass_d3D_noEoP",
    "Scattered electrons with nearest cluster within d_{3D} cut (no E/p);distance cut [mm];Events",
    NBINS_CUT, 0, DCUT_MAX
)

# pass counts with E/p
h_pass_d3D_withEoP = ROOT.TH1F(
    "h_pass_d3D_withEoP",
    f"Scattered electrons with nearest cluster within d_{{3D}} cut ({EOP_MIN:.1f} < E/p < {EOP_MAX:.1f});distance cut [mm];Events",
    NBINS_CUT, 0, DCUT_MAX
)

h_fail_recoE = ROOT.TH1F(
    "h_fail_recoE",
    f"Reco energy of scattered electrons failing d_{{3D}} < {CHOSEN_CUT_MM:.0f} mm;Reco electron energy [GeV];Events",
    100, 0, 30
)

h_fail_mcE = ROOT.TH1F(
    "h_fail_mcE",
    f"MC energy of scattered electrons failing d_{{3D}} < {CHOSEN_CUT_MM:.0f} mm;MC electron energy [GeV];Events",
    100, 0, 30
)

# parent-status histograms
h_parent_status_all = ROOT.TH1F(
    "h_parent_status_all",
    "Parent generator status of all electrons with status=1 and PDG=11;Parent generator status;Electrons",
    60, -0.5, 59.5
)

h_parent_status_selected = ROOT.TH1F(
    "h_parent_status_selected",
    "Parent generator status of selected scattered electrons (parent status != 2);Parent generator status;Events",
    60, -0.5, 59.5
)

h_parent_pdg_all = ROOT.TH1F(
    "h_parent_pdg_all",
    "Parent PDG of all electrons with status=1 and PDG=11;Parent PDG;Electrons",
    200, -100, 100
)

# ---------------------------------
# Counters and bookkeeping
# ---------------------------------
n_total_events = 0
n_all_status1_electrons = 0
n_events_with_scattered_electron = 0
n_scattered_with_reco = 0
n_scattered_with_track = 0
n_scattered_with_projection = 0
n_scattered_with_ecal_point = 0

n_scattered_with_cluster_noEoP = 0
n_scattered_with_cluster_withEoP = 0

scattered_records = []
failed_records = []

# ---------------------------------
# Event loop
# ---------------------------------
for ievt, event in enumerate(chain):
    n_total_events += 1

    required = [
        "MCParticles",
        "_MCParticles_parents",
        "CalorimeterTrackProjections",
        "_CalorimeterTrackProjections_track",
        "_CalorimeterTrackProjections_points",
        "EcalBarrelClusters",
        "ReconstructedParticles",
        "_ReconstructedParticles_tracks",
        "_ReconstructedParticleAssociations_rec",
        "_ReconstructedParticleAssociations_sim",
    ]
    if any(not hasattr(event, name) for name in required):
        continue

    mc_particles         = event.MCParticles
    parent_refs          = event._MCParticles_parents
    projs                = event.CalorimeterTrackProjections
    proj_tracks          = event._CalorimeterTrackProjections_track
    points               = event._CalorimeterTrackProjections_points
    clusters             = event.EcalBarrelClusters
    reco_particles       = event.ReconstructedParticles
    reco_particle_tracks = event._ReconstructedParticles_tracks
    assoc_rec            = event._ReconstructedParticleAssociations_rec
    assoc_sim            = event._ReconstructedParticleAssociations_sim

    # ---------------------------------
    # Fill histogram WITHOUT parent cut:
    # all electrons with status=1 and PDG=11
    # ---------------------------------
    for part in mc_particles:
        if int(part.generatorStatus) == 1 and int(part.PDG) == 11:
            n_all_status1_electrons += 1

            parent_status = get_first_parent_status(part, mc_particles, parent_refs)
            parent_pdg = get_first_parent_pdg(part, mc_particles, parent_refs)

            if parent_status is not None:
                h_parent_status_all.Fill(parent_status)
            if parent_pdg is not None:
                h_parent_pdg_all.Fill(parent_pdg)

    # ---------------------------------
    # Step 1: find scattered electron in MC with parent_status != 2
    # ---------------------------------
    scattered_mc_index = find_first_scattered_electron_mc_index(mc_particles, parent_refs)
    if scattered_mc_index is None:
        continue

    n_events_with_scattered_electron += 1
    scattered_mc = mc_particles[scattered_mc_index]
    mc_energy = particle_energy(scattered_mc)

    parent_status_selected = get_first_parent_status(scattered_mc, mc_particles, parent_refs)
    if parent_status_selected is not None:
        h_parent_status_selected.Fill(parent_status_selected)

    # ---------------------------------
    # Step 2: find associated reconstructed particle
    # ---------------------------------
    reco_index = find_reco_index_for_mc_index(scattered_mc_index, assoc_rec, assoc_sim)
    if reco_index is None or reco_index < 0 or reco_index >= len(reco_particles):
        failed_records.append({
            "event": ievt,
            "reason": "no_reco_match",
            "nearest_d3_noEoP": None,
            "nearest_d3_withEoP": None,
            "reco_energy": None,
            "mc_energy": mc_energy,
            "parent_status": parent_status_selected,
        })
        h_fail_mcE.Fill(mc_energy)
        continue

    n_scattered_with_reco += 1
    reco_obj = reco_particles[reco_index]
    reco_energy = particle_energy(reco_obj)
    reco_p = particle_p(reco_obj)

    # ---------------------------------
    # Step 3: get the track ref for that reco particle
    # ---------------------------------
    target_track_ref = find_track_ref_for_reco_index(reco_index, reco_particles, reco_particle_tracks)
    if target_track_ref is None:
        failed_records.append({
            "event": ievt,
            "reason": "no_track_on_reco_particle",
            "nearest_d3_noEoP": None,
            "nearest_d3_withEoP": None,
            "reco_energy": reco_energy,
            "mc_energy": mc_energy,
            "parent_status": parent_status_selected,
        })
        h_fail_recoE.Fill(reco_energy)
        h_fail_mcE.Fill(mc_energy)
        continue

    n_scattered_with_track += 1

    # ---------------------------------
    # Step 4: find the matching projected track
    # ---------------------------------
    matched_proj_index = None
    for i in range(len(projs)):
        if same_track_ref(proj_tracks[i], target_track_ref):
            matched_proj_index = i
            break

    if matched_proj_index is None:
        failed_records.append({
            "event": ievt,
            "reason": "no_projection_for_scattered_electron_track",
            "nearest_d3_noEoP": None,
            "nearest_d3_withEoP": None,
            "reco_energy": reco_energy,
            "mc_energy": mc_energy,
            "parent_status": parent_status_selected,
        })
        h_fail_recoE.Fill(reco_energy)
        h_fail_mcE.Fill(mc_energy)
        continue

    n_scattered_with_projection += 1

    # ---------------------------------
    # Step 5: find ECAL projection point
    # ---------------------------------
    proj = projs[matched_proj_index]
    b = int(proj.points_begin)
    e = int(proj.points_end)

    ecal_proj_pt = None
    for j in range(b, e):
        pt = points[j]
        if int(pt.system) == ECAL_SYSTEM:
            ecal_proj_pt = (
                float(pt.position.x),
                float(pt.position.y),
                float(pt.position.z)
            )
            break

    if ecal_proj_pt is None:
        failed_records.append({
            "event": ievt,
            "reason": "no_ecal_projection_point",
            "nearest_d3_noEoP": None,
            "nearest_d3_withEoP": None,
            "reco_energy": reco_energy,
            "mc_energy": mc_energy,
            "parent_status": parent_status_selected,
        })
        h_fail_recoE.Fill(reco_energy)
        h_fail_mcE.Fill(mc_energy)

        scattered_records.append({
            "event": ievt,
            "nearest_d3_noEoP": None,
            "nearest_d3_withEoP": None,
            "has_cluster_noEoP": False,
            "has_cluster_withEoP": False,
            "has_ecal_point": False,
            "reco_energy": reco_energy,
            "mc_energy": mc_energy,
            "parent_status": parent_status_selected,
        })
        continue

    n_scattered_with_ecal_point += 1

    px, py, pz = ecal_proj_pt

    # ---------------------------------
    # Step 6a: nearest reconstructed ECAL cluster WITHOUT E/p cut
    # ---------------------------------
    best_d3_noEoP = None
    best_dxy_noEoP = None

    if len(clusters) > 0:
        best_d3_tmp = 1e30
        best_dxy_tmp = 1e30
        found_noEoP = False

        for cl in clusters:
            cx = float(cl.position.x)
            cy = float(cl.position.y)
            cz = float(cl.position.z)

            dxy = math.sqrt((cx - px) ** 2 + (cy - py) ** 2)
            d3  = math.sqrt((cx - px) ** 2 + (cy - py) ** 2 + (cz - pz) ** 2)

            if d3 < best_d3_tmp:
                best_d3_tmp = d3
                best_dxy_tmp = dxy
                found_noEoP = True

        if found_noEoP:
            best_d3_noEoP = best_d3_tmp
            best_dxy_noEoP = best_dxy_tmp
            n_scattered_with_cluster_noEoP += 1
            h_d3D_noEoP.Fill(best_d3_noEoP)

    # ---------------------------------
    # Step 6b: nearest reconstructed ECAL cluster WITH E/p cut
    # ---------------------------------
    best_d3_withEoP = None
    best_dxy_withEoP = None
    best_eop = None

    if len(clusters) > 0 and reco_p > 0:
        best_d3_tmp = 1e30
        best_dxy_tmp = 1e30
        best_eop_tmp = None
        found_withEoP = False

        for cl in clusters:
            cx = float(cl.position.x)
            cy = float(cl.position.y)
            cz = float(cl.position.z)
            cE = float(cl.energy)

            dxy = math.sqrt((cx - px) ** 2 + (cy - py) ** 2)
            d3  = math.sqrt((cx - px) ** 2 + (cy - py) ** 2 + (cz - pz) ** 2)

            eop = cE / reco_p
            h_eop_all.Fill(eop)

            if not (EOP_MIN < eop < EOP_MAX):
                continue

            if d3 < best_d3_tmp:
                best_d3_tmp = d3
                best_dxy_tmp = dxy
                best_eop_tmp = eop
                found_withEoP = True

        if found_withEoP:
            best_d3_withEoP = best_d3_tmp
            best_dxy_withEoP = best_dxy_tmp
            best_eop = best_eop_tmp
            n_scattered_with_cluster_withEoP += 1
            h_d3D_withEoP.Fill(best_d3_withEoP)
            h_eop_selected.Fill(best_eop)

    scattered_records.append({
        "event": ievt,
        "nearest_d3_noEoP": best_d3_noEoP,
        "nearest_d3_withEoP": best_d3_withEoP,
        "nearest_dxy_noEoP": best_dxy_noEoP,
        "nearest_dxy_withEoP": best_dxy_withEoP,
        "has_cluster_noEoP": (best_d3_noEoP is not None),
        "has_cluster_withEoP": (best_d3_withEoP is not None),
        "has_ecal_point": True,
        "reco_energy": reco_energy,
        "mc_energy": mc_energy,
        "parent_status": parent_status_selected,
        "best_eop": best_eop,
    })

    # denominator for fraction:
    # every scattered electron with an ECAL projection point
    for ibin in range(1, NBINS_CUT + 1):
        dcut = h_total.GetBinCenter(ibin)
        h_total.Fill(dcut)

        if best_d3_noEoP is not None and best_d3_noEoP < dcut:
            h_pass_d3D_noEoP.Fill(dcut)

        if best_d3_withEoP is not None and best_d3_withEoP < dcut:
            h_pass_d3D_withEoP.Fill(dcut)

    # failing chosen cut WITH E/p version
    passed_cut_withEoP = (best_d3_withEoP is not None and best_d3_withEoP < CHOSEN_CUT_MM)
    if not passed_cut_withEoP:
        reason = "no_cluster_after_EoP" if best_d3_withEoP is None else "cluster_too_far_after_EoP"
        failed_records.append({
            "event": ievt,
            "reason": reason,
            "nearest_d3_noEoP": best_d3_noEoP,
            "nearest_d3_withEoP": best_d3_withEoP,
            "reco_energy": reco_energy,
            "mc_energy": mc_energy,
            "parent_status": parent_status_selected,
            "best_eop": best_eop,
        })
        h_fail_recoE.Fill(reco_energy)
        h_fail_mcE.Fill(mc_energy)

# ---------------------------------
# Build fraction histograms
# ---------------------------------
h_frac_d3D_noEoP = h_pass_d3D_noEoP.Clone("h_frac_d3D_noEoP")
h_frac_d3D_noEoP.SetTitle("Fraction of scattered electrons with nearby reconstructed ECal cluster (no E/p cut);d_{3D} cut [mm];Fraction")
h_frac_d3D_noEoP.Divide(h_total)

h_frac_d3D_withEoP = h_pass_d3D_withEoP.Clone("h_frac_d3D_withEoP")
h_frac_d3D_withEoP.SetTitle(
    f"Fraction of scattered electrons with nearby reconstructed ECal cluster ({EOP_MIN:.1f} < E/p < {EOP_MAX:.1f});d_{{3D}} cut [mm];Fraction"
)
h_frac_d3D_withEoP.Divide(h_total)

# ---------------------------------
# Plateau information
# ---------------------------------
plateau_bin_noEoP = h_frac_d3D_noEoP.GetMaximumBin()
plateau_fraction_noEoP = h_frac_d3D_noEoP.GetBinContent(plateau_bin_noEoP)
plateau_cut_noEoP = h_frac_d3D_noEoP.GetBinCenter(plateau_bin_noEoP)

first_plateau_bin_noEoP = None
for ibin in range(1, h_frac_d3D_noEoP.GetNbinsX() + 1):
    if abs(h_frac_d3D_noEoP.GetBinContent(ibin) - plateau_fraction_noEoP) < 1e-12:
        first_plateau_bin_noEoP = ibin
        break
first_plateau_cut_noEoP = h_frac_d3D_noEoP.GetBinCenter(first_plateau_bin_noEoP) if first_plateau_bin_noEoP else None

plateau_bin_withEoP = h_frac_d3D_withEoP.GetMaximumBin()
plateau_fraction_withEoP = h_frac_d3D_withEoP.GetBinContent(plateau_bin_withEoP)
plateau_cut_withEoP = h_frac_d3D_withEoP.GetBinCenter(plateau_bin_withEoP)

first_plateau_bin_withEoP = None
for ibin in range(1, h_frac_d3D_withEoP.GetNbinsX() + 1):
    if abs(h_frac_d3D_withEoP.GetBinContent(ibin) - plateau_fraction_withEoP) < 1e-12:
        first_plateau_bin_withEoP = ibin
        break
first_plateau_cut_withEoP = h_frac_d3D_withEoP.GetBinCenter(first_plateau_bin_withEoP) if first_plateau_bin_withEoP else None

# ---------------------------------
# Failure list for chosen cut
# ---------------------------------
n_pass_cut_noEoP = 0
n_fail_cut_noEoP = 0
n_pass_cut_withEoP = 0
n_fail_cut_withEoP = 0

for rec in scattered_records:
    if not rec["has_ecal_point"]:
        n_fail_cut_noEoP += 1
        n_fail_cut_withEoP += 1
        continue

    d3_noEoP = rec["nearest_d3_noEoP"]
    d3_withEoP = rec["nearest_d3_withEoP"]

    if d3_noEoP is not None and d3_noEoP < CHOSEN_CUT_MM:
        n_pass_cut_noEoP += 1
    else:
        n_fail_cut_noEoP += 1

    if d3_withEoP is not None and d3_withEoP < CHOSEN_CUT_MM:
        n_pass_cut_withEoP += 1
    else:
        n_fail_cut_withEoP += 1

fraction_at_cut_noEoP = (
    n_pass_cut_noEoP / n_scattered_with_ecal_point
    if n_scattered_with_ecal_point > 0 else 0.0
)

fraction_at_cut_withEoP = (
    n_pass_cut_withEoP / n_scattered_with_ecal_point
    if n_scattered_with_ecal_point > 0 else 0.0
)

# ---------------------------------
# Print summary
# ---------------------------------
print("\n========== Summary ==========")
print(f"Total events                                 : {n_total_events}")
print(f"All electrons with status=1 and PDG=11       : {n_all_status1_electrons}")
print(f"Events with selected scattered electron      : {n_events_with_scattered_electron}")
print(f"Scattered electrons with reco match          : {n_scattered_with_reco}")
print(f"Scattered electrons with reco track          : {n_scattered_with_track}")
print(f"Scattered electrons with projection          : {n_scattered_with_projection}")
print(f"Scattered electrons with ECAL projection     : {n_scattered_with_ecal_point}")
print(f"Scattered electrons with cluster (no E/p)    : {n_scattered_with_cluster_noEoP}")
print(f"Scattered electrons with cluster (with E/p)  : {n_scattered_with_cluster_withEoP}")

if h_d3D_noEoP.GetEntries() > 0:
    print(f"Mean nearest d3D (no E/p)                    : {h_d3D_noEoP.GetMean():.3f} mm")
if h_d3D_withEoP.GetEntries() > 0:
    print(f"Mean nearest d3D (with E/p)                  : {h_d3D_withEoP.GetMean():.3f} mm")

print("\n========== Parent-status summary ==========")
status_counts = {}
for ibin in range(1, h_parent_status_all.GetNbinsX() + 1):
    count = int(h_parent_status_all.GetBinContent(ibin))
    if count > 0:
        status_val = int(h_parent_status_all.GetBinCenter(ibin))
        status_counts[status_val] = count

for status_val in sorted(status_counts):
    print(f"Parent generator status {status_val:2d} : {status_counts[status_val]}")

bin2_all = h_parent_status_all.FindBin(2)
n_parent2_before = int(h_parent_status_all.GetBinContent(bin2_all))
print(f"Electrons removed by parent_status != 2 cut : {n_parent2_before}")

print("\n========== Plateau info ==========")
print(f"No E/p cut  - plateau fraction             : {plateau_fraction_noEoP:.6f}")
print(f"No E/p cut  - first plateau cut            : {first_plateau_cut_noEoP:.3f} mm")
print(f"With E/p    - plateau fraction             : {plateau_fraction_withEoP:.6f}")
print(f"With E/p    - first plateau cut            : {first_plateau_cut_withEoP:.3f} mm")

print(f"\n========== At chosen cut = {CHOSEN_CUT_MM:.1f} mm ==========")
print(f"No E/p cut  - passing scattered electrons  : {n_pass_cut_noEoP}")
print(f"No E/p cut  - failing scattered electrons  : {n_fail_cut_noEoP}")
print(f"No E/p cut  - fraction at chosen cut       : {fraction_at_cut_noEoP:.6f}")
print(f"With E/p    - passing scattered electrons  : {n_pass_cut_withEoP}")
print(f"With E/p    - failing scattered electrons  : {n_fail_cut_withEoP}")
print(f"With E/p    - fraction at chosen cut       : {fraction_at_cut_withEoP:.6f}")

# ---------------------------------
# Canvases (ROOT only, no PNGs)
# ---------------------------------
c1 = ROOT.TCanvas("c1", "d3D without E/p", 900, 700)
h_d3D_noEoP.SetLineWidth(2)
h_d3D_noEoP.Draw("HIST")

c2 = ROOT.TCanvas("c2", "d3D with E/p", 900, 700)
h_d3D_withEoP.SetLineWidth(2)
h_d3D_withEoP.SetLineColor(ROOT.kBlue + 1)
h_d3D_withEoP.Draw("HIST")

c3 = ROOT.TCanvas("c3", "fraction vs distance cut", 900, 700)
h_frac_d3D_noEoP.SetLineWidth(2)
h_frac_d3D_noEoP.SetLineColor(ROOT.kBlack)
h_frac_d3D_noEoP.SetMinimum(0.0)
h_frac_d3D_noEoP.SetMaximum(1.05)
h_frac_d3D_noEoP.Draw("HIST")

h_frac_d3D_withEoP.SetLineWidth(2)
h_frac_d3D_withEoP.SetLineColor(ROOT.kRed)
h_frac_d3D_withEoP.Draw("HIST SAME")

leg = ROOT.TLegend(0.52, 0.73, 0.88, 0.88)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.AddEntry(h_frac_d3D_noEoP, "No E/p cut", "l")
leg.AddEntry(h_frac_d3D_withEoP, f"{EOP_MIN:.1f} < E/p < {EOP_MAX:.1f}", "l")
leg.Draw()

c4 = ROOT.TCanvas("c4", "failing scattered electron energies", 1200, 500)
c4.Divide(2, 1)
c4.cd(1)
h_fail_recoE.SetLineWidth(2)
h_fail_recoE.Draw("HIST")
c4.cd(2)
h_fail_mcE.SetLineWidth(2)
h_fail_mcE.Draw("HIST")

c5 = ROOT.TCanvas("c5", "parent status histograms", 1200, 500)
c5.Divide(2, 1)
c5.cd(1)
h_parent_status_all.SetLineWidth(2)
h_parent_status_all.Draw("HIST")
c5.cd(2)
h_parent_status_selected.SetLineWidth(2)
h_parent_status_selected.SetLineColor(ROOT.kBlue + 1)
h_parent_status_selected.Draw("HIST")

c6 = ROOT.TCanvas("c6", "parent PDG all electrons", 900, 700)
h_parent_pdg_all.SetLineWidth(2)
h_parent_pdg_all.Draw("HIST")

c7 = ROOT.TCanvas("c7", "E over p", 1200, 500)
c7.Divide(2, 1)
c7.cd(1)
h_eop_all.SetLineWidth(2)
h_eop_all.Draw("HIST")
c7.cd(2)
h_eop_selected.SetLineWidth(2)
h_eop_selected.SetLineColor(ROOT.kBlue + 1)
h_eop_selected.Draw("HIST")

# ---------------------------------
# Save ROOT file
# ---------------------------------
out = ROOT.TFile(OUTROOT, "RECREATE")
h_d3D_noEoP.Write()
h_d3D_withEoP.Write()
h_eop_all.Write()
h_eop_selected.Write()
h_total.Write()
h_pass_d3D_noEoP.Write()
h_pass_d3D_withEoP.Write()
h_frac_d3D_noEoP.Write()
h_frac_d3D_withEoP.Write()
h_fail_recoE.Write()
h_fail_mcE.Write()
h_parent_status_all.Write()
h_parent_status_selected.Write()
h_parent_pdg_all.Write()
c1.Write()
c2.Write()
c3.Write()
c4.Write()
c5.Write()
c6.Write()
c7.Write()
out.Close()

print(f"\n✅ Wrote ROOT file: {OUTROOT}")
