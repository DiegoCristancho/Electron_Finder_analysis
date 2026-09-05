import ROOT
import math

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(1110)

# ---------------------------------
# Input / output
# ---------------------------------
LISTFILE = "single_electron.txt"
EVENT_TREE = "events"

OUTROOT = "electron_track_cluster_fraction.root"
OUTPNG1 = "nearest_cluster_d3D_electron_tracks.png"
OUTPNG2 = "nearest_cluster_dxy_electron_tracks.png"
OUTPNG3 = "fraction_vs_dcut.png"

ECAL_SYSTEM = 101

# ---------------------------------
# Helpers
# ---------------------------------
def find_reco_particle_for_track(track_ref, reco_particles, reco_particle_tracks):
    trk_coll = int(track_ref.collectionID)
    trk_idx  = int(track_ref.index)

    for ireco, reco in enumerate(reco_particles):
        b = int(reco.tracks_begin)
        e = int(reco.tracks_end)

        for j in range(b, e):
            tref = reco_particle_tracks[j]
            if int(tref.collectionID) == trk_coll and int(tref.index) == trk_idx:
                return ireco

    return None


def find_mc_pdg_for_reco_index(reco_index, assoc_rec, assoc_sim, mc_particles):
    if reco_index is None:
        return None

    for ia in range(len(assoc_rec)):
        rec_ref = assoc_rec[ia]
        if int(rec_ref.index) == int(reco_index):
            sim_ref = assoc_sim[ia]
            sim_idx = int(sim_ref.index)
            if 0 <= sim_idx < len(mc_particles):
                return int(mc_particles[sim_idx].PDG)

    return None


# ---------------------------------
# Build chain
# ---------------------------------
chain = ROOT.TChain(EVENT_TREE)
with open(LISTFILE) as f:
    for line in f:
        fn = line.strip()
        if fn:
            chain.Add(fn)

print(f"✅ Added {chain.GetNtrees()} files, total {chain.GetEntries()} events")

# ---------------------------------
# Histograms
# ---------------------------------
h_d3D = ROOT.TH1F(
    "h_d3D",
    "Nearest ECal cluster distance for truth-matched electron tracks;d_{3D} [mm];Electron tracks",
    100, 0, 200
)

h_dxy = ROOT.TH1F(
    "h_dxy",
    "Nearest ECal cluster transverse distance for truth-matched electron tracks;d_{xy} [mm];Electron tracks",
    100, 0, 200
)

# Fraction vs distance cut
NBINS_CUT = 100
DCUT_MAX = 200.0

h_total = ROOT.TH1F(
    "h_total",
    "Total truth-matched electron tracks;distance cut [mm];Tracks",
    NBINS_CUT, 0, DCUT_MAX
)

h_pass_d3D = ROOT.TH1F(
    "h_pass_d3D",
    "Electron tracks with nearest cluster within d_{3D} cut;distance cut [mm];Tracks",
    NBINS_CUT, 0, DCUT_MAX
)

h_pass_dxy = ROOT.TH1F(
    "h_pass_dxy",
    "Electron tracks with nearest cluster within d_{xy} cut;distance cut [mm];Tracks",
    NBINS_CUT, 0, DCUT_MAX
)

# ---------------------------------
# Counters
# ---------------------------------
n_total_events = 0
n_events_with_projs = 0
n_events_with_clusters = 0
n_electron_tracks = 0
n_electron_tracks_with_ecal_point = 0
n_electron_tracks_with_cluster = 0

# ---------------------------------
# Event loop
# ---------------------------------
for ievt, event in enumerate(chain):
    n_total_events += 1

    required = [
        "CalorimeterTrackProjections",
        "_CalorimeterTrackProjections_track",
        "_CalorimeterTrackProjections_points",
        "EcalBarrelClusters",
        "ReconstructedParticles",
        "_ReconstructedParticles_tracks",
        "_ReconstructedParticleAssociations_rec",
        "_ReconstructedParticleAssociations_sim",
        "MCParticles",
    ]
    if any(not hasattr(event, name) for name in required):
        continue

    projs                = event.CalorimeterTrackProjections
    tracks               = event._CalorimeterTrackProjections_track
    points               = event._CalorimeterTrackProjections_points
    clusters             = event.EcalBarrelClusters
    reco_particles       = event.ReconstructedParticles
    reco_particle_tracks = event._ReconstructedParticles_tracks
    assoc_rec            = event._ReconstructedParticleAssociations_rec
    assoc_sim            = event._ReconstructedParticleAssociations_sim
    mc_particles         = event.MCParticles

    if len(projs) > 0:
        n_events_with_projs += 1
    if len(clusters) > 0:
        n_events_with_clusters += 1

    # Loop over all projected tracks in the event
    for i in range(len(projs)):
        proj = projs[i]
        trk  = tracks[i]

        # Truth-match the track only
        reco_index = find_reco_particle_for_track(trk, reco_particles, reco_particle_tracks)
        truth_pdg  = find_mc_pdg_for_reco_index(reco_index, assoc_rec, assoc_sim, mc_particles)

        # Require electron track
        if truth_pdg is None or abs(truth_pdg) != 11:
            continue

        n_electron_tracks += 1

        # Find ECAL projection point for this track
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
            continue

        n_electron_tracks_with_ecal_point += 1

        if len(clusters) == 0:
            continue

        px, py, pz = ecal_proj_pt

        # Nearest reconstructed cluster
        best_d3 = 1e30
        best_dxy = 1e30

        for cl in clusters:
            cx = float(cl.position.x)
            cy = float(cl.position.y)
            cz = float(cl.position.z)

            dxy = math.sqrt((cx - px)**2 + (cy - py)**2)
            d3  = math.sqrt((cx - px)**2 + (cy - py)**2 + (cz - pz)**2)

            if d3 < best_d3:
                best_d3 = d3
                best_dxy = dxy

        n_electron_tracks_with_cluster += 1

        h_d3D.Fill(best_d3)
        h_dxy.Fill(best_dxy)

        # Fill pass/fail cumulatively vs cut
        for ibin in range(1, NBINS_CUT + 1):
            dcut = h_total.GetBinCenter(ibin)
            h_total.Fill(dcut)

            if best_d3 < dcut:
                h_pass_d3D.Fill(dcut)
            if best_dxy < dcut:
                h_pass_dxy.Fill(dcut)

# ---------------------------------
# Build fraction histograms
# ---------------------------------
h_frac_d3D = h_pass_d3D.Clone("h_frac_d3D")
h_frac_d3D.SetTitle("Fraction of electron tracks with nearby reconstructed ECal cluster;d_{3D} cut [mm];Fraction")
h_frac_d3D.Divide(h_total)

h_frac_dxy = h_pass_dxy.Clone("h_frac_dxy")
h_frac_dxy.SetTitle("Fraction of electron tracks with nearby reconstructed ECal cluster;d_{xy} cut [mm];Fraction")
h_frac_dxy.Divide(h_total)

# ---------------------------------
# Print summary
# ---------------------------------
print("\n========== Summary ==========")
print(f"Total events                            : {n_total_events}")
print(f"Events with projections                 : {n_events_with_projs}")
print(f"Events with clusters                    : {n_events_with_clusters}")
print(f"Truth-matched electron tracks           : {n_electron_tracks}")
print(f"Electron tracks with ECAL point         : {n_electron_tracks_with_ecal_point}")
print(f"Electron tracks with >=1 cluster        : {n_electron_tracks_with_cluster}")

if h_d3D.GetEntries() > 0:
    print(f"Mean nearest d3D                        : {h_d3D.GetMean():.3f} mm")
    print(f"Mean nearest dxy                        : {h_dxy.GetMean():.3f} mm")

# Some example thresholds
for thr in [5, 10, 20, 40, 60, 100]:
    b = h_frac_d3D.FindBin(thr)
    print(f"Fraction with d3D < {thr:3d} mm         : {h_frac_d3D.GetBinContent(b):.4f}")

# ---------------------------------
# Draw histograms
# ---------------------------------
c1 = ROOT.TCanvas("c1", "nearest d3D", 900, 700)
h_d3D.SetLineWidth(2)
h_d3D.Draw("HIST")
c1.SaveAs(OUTPNG1)

c2 = ROOT.TCanvas("c2", "nearest dxy", 900, 700)
h_dxy.SetLineWidth(2)
h_dxy.Draw("HIST")
c2.SaveAs(OUTPNG2)

c3 = ROOT.TCanvas("c3", "fraction vs distance cut", 900, 700)
h_frac_d3D.SetLineWidth(2)
h_frac_d3D.SetLineColor(ROOT.kBlue + 1)
h_frac_d3D.SetMinimum(0.0)
h_frac_d3D.SetMaximum(1.05)
h_frac_d3D.Draw("HIST")

h_frac_dxy.SetLineWidth(2)
h_frac_dxy.SetLineColor(ROOT.kRed)
h_frac_dxy.Draw("HIST SAME")

leg = ROOT.TLegend(0.58, 0.75, 0.88, 0.88)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.AddEntry(h_frac_d3D, "Using d_{3D} cut", "l")
leg.AddEntry(h_frac_dxy, "Using d_{xy} cut", "l")
leg.Draw()

c3.SaveAs(OUTPNG3)

# ---------------------------------
# Save ROOT file
# ---------------------------------
out = ROOT.TFile(OUTROOT, "RECREATE")
h_d3D.Write()
h_dxy.Write()
h_total.Write()
h_pass_d3D.Write()
h_pass_dxy.Write()
h_frac_d3D.Write()
h_frac_dxy.Write()
c1.Write()
c2.Write()
c3.Write()
out.Close()

print(f"\n✅ Wrote ROOT file: {OUTROOT}")
print(f"✅ Wrote image:     {OUTPNG1}")
print(f"✅ Wrote image:     {OUTPNG2}")
print(f"✅ Wrote image:     {OUTPNG3}")
