import ROOT
import math

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# ---------------------------------
# Inputs
# ---------------------------------

# Set this to True if you want to use the txt file list.
# Set this to False if you want to use only one ROOT file.
USE_LISTFILE = False

# Old option: read many ROOT files from this txt list
LISTFILE = "list_18x275_Q2_100_26_02.txt"

# New option: read one ROOT file directly
#INPUTROOT = "pythia8NCDIS_18x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_1.1601.eicrecon.edm4eic.root"

INPUTROOT = "pythia8NCDIS_26_718x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_1.0008.eicrecon.edm4eic.root"

EVENT_TREE = "events"
TARGET_EVENT = 1
OUTROOT = "one_event_track_projection_points.root"

# ---------------------------------
# Helper functions
# ---------------------------------
def find_reco_particle_for_track(track_ref, reco_particles, reco_particle_tracks):
    """
    Given one track reference from _CalorimeterTrackProjections_track,
    find which ReconstructedParticle owns that track.
    """
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
    """
    Given the index of a reconstructed particle in ReconstructedParticles,
    find the associated MC particle PDG using ReconstructedParticleAssociations.
    """
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

if USE_LISTFILE:
    # Old version: read ROOT files from a txt list
    with open(LISTFILE) as f:
        for line in f:
            fn = line.strip()
            if fn:
                chain.Add(fn)

    print(f"Added {chain.GetNtrees()} files from txt list, total {chain.GetEntries()} events")

else:
    # New version: read one ROOT file directly
    chain.Add(INPUTROOT)

    print(f"Added {chain.GetNtrees()} file, total {chain.GetEntries()} events")

# ---------------------------------
# Grab one event
# ---------------------------------
target = None
for ievt, event in enumerate(chain):
    if ievt == TARGET_EVENT:
        target = event
        break

if target is None:
    raise RuntimeError(f"Could not find event {TARGET_EVENT}")

# Required branches
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
for name in required:
    if not hasattr(target, name):
        raise RuntimeError(f"Missing required branch: {name}")

projs                = target.CalorimeterTrackProjections
tracks               = target._CalorimeterTrackProjections_track
points               = target._CalorimeterTrackProjections_points
clusters             = target.EcalBarrelClusters
reco_particles       = target.ReconstructedParticles
reco_particle_tracks = target._ReconstructedParticles_tracks
assoc_rec            = target._ReconstructedParticleAssociations_rec
assoc_sim            = target._ReconstructedParticleAssociations_sim
mc_particles         = target.MCParticles

print(f"\n========== Event {TARGET_EVENT} ==========")
print(f"n projections                    = {len(projs)}")
print(f"n track refs                     = {len(tracks)}")
print(f"n total projection points        = {len(points)}")
print(f"n ECal clusters                  = {len(clusters)}")
print(f"n ReconstructedParticles         = {len(reco_particles)}")
print(f"n _ReconstructedParticles_tracks = {len(reco_particle_tracks)}")
print(f"n Reco->MC associations          = {len(assoc_rec)}")
print(f"n MCParticles                    = {len(mc_particles)}\n")

# ---------------------------------
# Colors
# ---------------------------------
colors = [
    ROOT.kRed + 1,
    ROOT.kBlue + 1,
    ROOT.kGreen + 2,
    ROOT.kMagenta + 1,
    ROOT.kOrange + 7,
    ROOT.kCyan + 1,
    ROOT.kViolet + 1,
    ROOT.kBlack,
]

# ---------------------------------
# Find global ranges
# ---------------------------------
xmin = ymin = zmin =  1e30
xmax = ymax = zmax = -1e30
rmin =  1e30
rmax = -1e30

track_data = []

for i in range(len(projs)):
    proj = projs[i]
    trk  = tracks[i]

    b = int(proj.points_begin)
    e = int(proj.points_end)

    # ---- find reco particle and truth PDG for this projected track
    reco_index = find_reco_particle_for_track(trk, reco_particles, reco_particle_tracks)
    truth_pdg  = find_mc_pdg_for_reco_index(reco_index, assoc_rec, assoc_sim, mc_particles)

    pts = []
    print(f"Projection {i}")
    print(f"  track.index        = {int(trk.index)}")
    print(f"  track.collectionID = {int(trk.collectionID)}")
    print(f"  points_begin       = {b}")
    print(f"  points_end         = {e}")
    print(f"  n points           = {e-b}")
    print(f"  matched reco index = {reco_index}")
    print(f"  truth PDG          = {truth_pdg}")

    for j in range(b, e):
        pt = points[j]
        x = float(pt.position.x)
        y = float(pt.position.y)
        z = float(pt.position.z)
        r = math.sqrt(x*x + y*y)
        sysid = int(pt.system)
        surf  = int(pt.surface)

        pts.append((x, y, z, r, sysid, surf))

        xmin = min(xmin, x)
        xmax = max(xmax, x)
        ymin = min(ymin, y)
        ymax = max(ymax, y)
        zmin = min(zmin, z)
        zmax = max(zmax, z)
        rmin = min(rmin, r)
        rmax = max(rmax, r)

        print(
            f"    point {j}: "
            f"system={sysid}, surface={surf}, "
            f"x={x:.3f}, y={y:.3f}, z={z:.3f}, r={r:.3f}"
        )

    print()
    track_data.append({
        "proj_index": i,
        "track_index": int(trk.index),
        "collectionID": int(trk.collectionID),
        "reco_index": reco_index,
        "truth_pdg": truth_pdg,
        "points_begin": b,
        "points_end": e,
        "pts": pts
    })

# ---------------------------------
# Also include cluster positions in the plotting range
# ---------------------------------
cluster_xyz = []
print("EcalBarrelClusters")
for icl, cl in enumerate(clusters):
    x = float(cl.position.x)
    y = float(cl.position.y)
    z = float(cl.position.z)
    r = math.sqrt(x*x + y*y)

    cluster_xyz.append((x, y, z, r))

    xmin = min(xmin, x)
    xmax = max(xmax, x)
    ymin = min(ymin, y)
    ymax = max(ymax, y)
    zmin = min(zmin, z)
    zmax = max(zmax, z)
    rmin = min(rmin, r)
    rmax = max(rmax, r)

    print(f"  cluster {icl}: x={x:.3f}, y={y:.3f}, z={z:.3f}, r={r:.3f}")

margin = 20.0
xmin -= margin
xmax += margin
ymin -= margin
ymax += margin
zmin -= margin
zmax += margin
rmin = max(0.0, rmin - margin)
rmax += margin

# ---------------------------------
# 3D canvas
# ---------------------------------
c3d = ROOT.TCanvas("c3d", f"Event {TARGET_EVENT} 3D", 1100, 900)

frame3d = ROOT.TH3F(
    "frame3d",
    f"Event {TARGET_EVENT} track projection points; x [mm]; y [mm]; z [mm]",
    10, xmin, xmax,
    10, ymin, ymax,
    10, zmin, zmax
)
frame3d.SetStats(0)
frame3d.Draw()

saved = [frame3d]

leg3d = ROOT.TLegend(0.62, 0.60, 0.92, 0.90)
leg3d.SetBorderSize(0)
leg3d.SetFillStyle(0)
leg3d.SetTextSize(0.025)

for i, td in enumerate(track_data):
    pts = td["pts"]
    color = colors[i % len(colors)]

    line = ROOT.TPolyLine3D(len(pts))
    mark = ROOT.TPolyMarker3D(len(pts))

    for k, (x, y, z, r, sysid, surf) in enumerate(pts):
        line.SetPoint(k, x, y, z)
        mark.SetPoint(k, x, y, z)

    line.SetLineColor(color)
    line.SetLineWidth(2)

    mark.SetMarkerColor(color)
    mark.SetMarkerStyle(20)
    mark.SetMarkerSize(1.0)

    line.Draw("same")
    mark.Draw("same")

    pdg_label = td["truth_pdg"] if td["truth_pdg"] is not None else "?"
    leg3d.AddEntry(mark, f"track {td['track_index']} (PDG {pdg_label})", "p")

    saved.extend([line, mark])

leg3d.Draw()
saved.append(leg3d)

c3d.SaveAs("one_event_track_projection_points_3D.png")

# ---------------------------------
# X-Y canvas: track projection points only
# ---------------------------------
cxy = ROOT.TCanvas("cxy", f"Event {TARGET_EVENT} XY", 1100, 900)

framexy = ROOT.TH2F(
    "framexy",
    f"Event {TARGET_EVENT} track projection points; x [mm]; y [mm]",
    200, xmin, xmax,
    200, ymin, ymax
)
framexy.SetStats(0)
framexy.Draw()

saved.append(framexy)

legxy = ROOT.TLegend(0.62, 0.60, 0.92, 0.90)
legxy.SetBorderSize(0)
legxy.SetFillStyle(0)
legxy.SetTextSize(0.025)

for i, td in enumerate(track_data):
    pts = td["pts"]
    color = colors[i % len(colors)]

    g = ROOT.TGraph(len(pts))
    for k, (x, y, z, r, sysid, surf) in enumerate(pts):
        g.SetPoint(k, x, y)

    g.SetLineColor(color)
    g.SetMarkerColor(color)
    g.SetLineWidth(2)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.1)

    g.Draw("LP SAME")

    pdg_label = td["truth_pdg"] if td["truth_pdg"] is not None else "?"
    legxy.AddEntry(g, f"track {td['track_index']} (PDG {pdg_label})", "lp")

    saved.append(g)

legxy.Draw()
saved.append(legxy)

cxy.SaveAs("one_event_track_projection_points_XY.png")

# ---------------------------------
# X-Y canvas: track projection points + ECal clusters
# ---------------------------------
cxyclus = ROOT.TCanvas("cxyclus", f"Event {TARGET_EVENT} XY with ECal clusters", 1100, 900)

framexyclus = ROOT.TH2F(
    "framexyclus",
    f"Event {TARGET_EVENT} track projection points and ECal clusters; x [mm]; y [mm]",
    200, xmin, xmax,
    200, ymin, ymax
)
framexyclus.SetStats(0)
framexyclus.Draw()

saved.append(framexyclus)

legxyclus = ROOT.TLegend(0.58, 0.56, 0.92, 0.90)
legxyclus.SetBorderSize(0)
legxyclus.SetFillStyle(0)
legxyclus.SetTextSize(0.025)

for i, td in enumerate(track_data):
    pts = td["pts"]
    color = colors[i % len(colors)]

    g = ROOT.TGraph(len(pts))
    for k, (x, y, z, r, sysid, surf) in enumerate(pts):
        g.SetPoint(k, x, y)

    g.SetLineColor(color)
    g.SetMarkerColor(color)
    g.SetLineWidth(2)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.1)

    g.Draw("LP SAME")

    pdg_label = td["truth_pdg"] if td["truth_pdg"] is not None else "?"
    legxyclus.AddEntry(g, f"track {td['track_index']} (PDG {pdg_label})", "lp")

    saved.append(g)

gclus = ROOT.TGraph(len(cluster_xyz))
for i, (x, y, z, r) in enumerate(cluster_xyz):
    gclus.SetPoint(i, x, y)

gclus.SetMarkerStyle(20)
gclus.SetMarkerSize(1.2)
gclus.SetMarkerColor(ROOT.kBlack)
gclus.SetLineColor(ROOT.kBlack)

gclus.Draw("P SAME")
saved.append(gclus)
legxyclus.AddEntry(gclus, "EcalBarrelClusters", "p")

legxyclus.Draw()
saved.append(legxyclus)

cxyclus.SaveAs("one_event_track_projection_points_XY_with_EcalClusters.png")

# ---------------------------------
# R-Z canvas: track projection points + ECal clusters
# ---------------------------------
crz = ROOT.TCanvas("crz", f"Event {TARGET_EVENT} RZ with ECal clusters", 1100, 900)

framerz = ROOT.TH2F(
    "framerz",
    f"Event {TARGET_EVENT} track projection points and ECal clusters; z [mm]; r = #sqrt{{x^2+y^2}} [mm]",
    200, zmin, zmax,
    200, rmin, rmax
)
framerz.SetStats(0)
framerz.Draw()

saved.append(framerz)

legrz = ROOT.TLegend(0.58, 0.56, 0.92, 0.90)
legrz.SetBorderSize(0)
legrz.SetFillStyle(0)
legrz.SetTextSize(0.025)

for i, td in enumerate(track_data):
    pts = td["pts"]
    color = colors[i % len(colors)]

    g = ROOT.TGraph(len(pts))
    for k, (x, y, z, r, sysid, surf) in enumerate(pts):
        g.SetPoint(k, z, r)

    g.SetLineColor(color)
    g.SetMarkerColor(color)
    g.SetLineWidth(2)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.1)

    g.Draw("LP SAME")

    pdg_label = td["truth_pdg"] if td["truth_pdg"] is not None else "?"
    legrz.AddEntry(g, f"track {td['track_index']} (PDG {pdg_label})", "lp")

    saved.append(g)

gclus_rz = ROOT.TGraph(len(cluster_xyz))
for i, (x, y, z, r) in enumerate(cluster_xyz):
    gclus_rz.SetPoint(i, z, r)

gclus_rz.SetMarkerStyle(20)
gclus_rz.SetMarkerSize(1.2)
gclus_rz.SetMarkerColor(ROOT.kBlack)
gclus_rz.SetLineColor(ROOT.kBlack)

gclus_rz.Draw("P SAME")
saved.append(gclus_rz)
legrz.AddEntry(gclus_rz, "EcalBarrelClusters", "p")

legrz.Draw()
saved.append(legrz)

crz.SaveAs("one_event_track_projection_points_RZ_with_EcalClusters.png")

# ---------------------------------
# Save ROOT file
# ---------------------------------
out = ROOT.TFile(OUTROOT, "RECREATE")
c3d.Write()
cxy.Write()
cxyclus.Write()
crz.Write()
for obj in saved:
    obj.Write()
out.Close()

print("Saved:")
print("  one_event_track_projection_points_3D.png")
print("  one_event_track_projection_points_XY.png")
print("  one_event_track_projection_points_XY_with_EcalClusters.png")
print("  one_event_track_projection_points_RZ_with_EcalClusters.png")
print(f"  {OUTROOT}")
