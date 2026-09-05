import ROOT
import math

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(1110)

# ---------------------------------
# Input / output
# ---------------------------------
LISTFILE = "single_electron.txt"
OUTROOT  = "track_cluster_distance_hist.root"
OUTPNG1  = "nearest_cluster_d3D.png"
OUTPNG2  = "nearest_cluster_dxy.png"
OUTPNG3  = "cluster_minus_projection_deltas.png"
OUTPNG4  = "cluster_minus_projection_spherical_deltas.png"

EVENT_TREE = "events"

ECAL_SYSTEM = 101

# ---------------------------------
# Helpers
# ---------------------------------
def wrap_phi(dphi):
    while dphi > math.pi:
        dphi -= 2.0 * math.pi
    while dphi <= -math.pi:
        dphi += 2.0 * math.pi
    return dphi

def safe_r3(x, y, z):
    return math.sqrt(x*x + y*y + z*z)

def safe_theta(x, y, z):
    r = safe_r3(x, y, z)
    if r <= 0:
        return 0.0
    c = z / r
    if c > 1.0:
        c = 1.0
    if c < -1.0:
        c = -1.0
    return math.acos(c)

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
h_d3D = ROOT.TH1F("h_d3D", "Nearest ECAL cluster distance;d_{3D} [mm];Events", 100, 0, 200)
h_dxy = ROOT.TH1F("h_dxy", "Nearest ECAL cluster transverse distance;d_{xy} [mm];Events", 100, 0, 200)

h_dx   = ROOT.TH1F("h_dx",   "x_{cluster} - x_{proj};x_{cluster} - x_{proj} [mm];Events", 100, -100, 100)
h_dy   = ROOT.TH1F("h_dy",   "y_{cluster} - y_{proj};y_{cluster} - y_{proj} [mm];Events", 100, -100, 100)
h_dz   = ROOT.TH1F("h_dz",   "z_{cluster} - z_{proj};z_{cluster} - z_{proj} [mm];Events", 100, -100, 100)
h_dphi = ROOT.TH1F("h_dphi", "#phi_{cluster} - #phi_{proj};#phi_{cluster} - #phi_{proj} [rad];Events", 100, -0.2, 0.2)

h_dr     = ROOT.TH1F("h_dr",     "r_{cluster} - r_{proj};r_{cluster} - r_{proj} [mm];Events", 100, -150, 150)
h_dtheta = ROOT.TH1F("h_dtheta", "#theta_{cluster} - #theta_{proj};#theta_{cluster} - #theta_{proj} [rad];Events", 100, -0.2, 0.2)

# Optional bookkeeping
h_nclus = ROOT.TH1F("h_nclus", "Number of ECal clusters per event;N clusters;Events", 10, 0, 10)
h_nproj = ROOT.TH1F("h_nproj", "Number of track projections per event;N projections;Events", 10, 0, 10)

# ---------------------------------
# Counters
# ---------------------------------
n_total = 0
n_with_cluster = 0
n_with_proj = 0
n_with_both = 0
n_filled = 0
n_no_ecal_proj_point = 0

# ---------------------------------
# Event loop
# ---------------------------------
for ievt, event in enumerate(chain):
    n_total += 1

    if not hasattr(event, "CalorimeterTrackProjections"):
        continue
    if not hasattr(event, "_CalorimeterTrackProjections_points"):
        continue
    if not hasattr(event, "EcalBarrelClusters"):
        continue

    projs    = event.CalorimeterTrackProjections
    points   = event._CalorimeterTrackProjections_points
    clusters = event.EcalBarrelClusters

    nproj = len(projs)
    nclus = len(clusters)

    h_nproj.Fill(nproj)
    h_nclus.Fill(nclus)

    if nclus > 0:
        n_with_cluster += 1
    if nproj > 0:
        n_with_proj += 1

    if nclus == 0 or nproj == 0:
        continue

    n_with_both += 1

    # ---------------------------------
    # Use the first projection in the event
    # ---------------------------------
    proj = projs[0]
    b = int(proj.points_begin)
    e = int(proj.points_end)

    ecal_proj_pt = None

    for j in range(b, e):
        pt = points[j]
        sysid = int(pt.system)

        if sysid == ECAL_SYSTEM:
            ecal_proj_pt = (
                float(pt.position.x),
                float(pt.position.y),
                float(pt.position.z)
            )
            break

    if ecal_proj_pt is None:
        n_no_ecal_proj_point += 1
        continue

    px, py, pz = ecal_proj_pt
    pphi = math.atan2(py, px)
    pr   = safe_r3(px, py, pz)
    ptheta = safe_theta(px, py, pz)

    # ---------------------------------
    # Find nearest ECAL cluster
    # ---------------------------------
    best_d3 = 1e30
    best_dxy = 1e30
    best_cluster = None

    for cl in clusters:
        cx = float(cl.position.x)
        cy = float(cl.position.y)
        cz = float(cl.position.z)

        dxy = math.sqrt((cx - px)**2 + (cy - py)**2)
        d3  = math.sqrt((cx - px)**2 + (cy - py)**2 + (cz - pz)**2)

        if d3 < best_d3:
            best_d3 = d3
            best_dxy = dxy
            best_cluster = (cx, cy, cz)

    h_d3D.Fill(best_d3)
    h_dxy.Fill(best_dxy)

    if best_cluster is not None:
        cx, cy, cz = best_cluster
        cphi = math.atan2(cy, cx)
        cr   = safe_r3(cx, cy, cz)
        ctheta = safe_theta(cx, cy, cz)

        dx = cx - px
        dy = cy - py
        dz = cz - pz
        dphi = wrap_phi(cphi - pphi)
        dr = cr - pr
        dtheta = ctheta - ptheta

        h_dx.Fill(dx)
        h_dy.Fill(dy)
        h_dz.Fill(dz)
        h_dphi.Fill(dphi)
        h_dr.Fill(dr)
        h_dtheta.Fill(dtheta)

    n_filled += 1

# ---------------------------------
# Print summary
# ---------------------------------
print("\n========== Summary ==========")
print(f"Total events                 : {n_total}")
print(f"Events with >=1 cluster      : {n_with_cluster}")
print(f"Events with >=1 projection   : {n_with_proj}")
print(f"Events with both             : {n_with_both}")
print(f"Events with ECAL proj point  : {n_with_both - n_no_ecal_proj_point}")
print(f"Events filled in histograms  : {n_filled}")
print(f"Skipped (no ECAL proj point) : {n_no_ecal_proj_point}")

if n_filled > 0:
    print(f"Mean d3D     = {h_d3D.GetMean():.3f} mm")
    print(f"Mean dxy     = {h_dxy.GetMean():.3f} mm")
    print(f"Mean dx      = {h_dx.GetMean():.3f} mm")
    print(f"Mean dy      = {h_dy.GetMean():.3f} mm")
    print(f"Mean dz      = {h_dz.GetMean():.3f} mm")
    print(f"Mean dphi    = {h_dphi.GetMean():.6f} rad")
    print(f"Mean dr      = {h_dr.GetMean():.3f} mm")
    print(f"Mean dtheta  = {h_dtheta.GetMean():.6f} rad")

# ---------------------------------
# Draw histograms
# ---------------------------------
c1 = ROOT.TCanvas("c1", "d3D", 900, 700)
h_d3D.SetLineWidth(2)
h_d3D.Draw("HIST")
c1.SaveAs(OUTPNG1)

c2 = ROOT.TCanvas("c2", "dxy", 900, 700)
h_dxy.SetLineWidth(2)
h_dxy.Draw("HIST")
c2.SaveAs(OUTPNG2)

c3 = ROOT.TCanvas("c3", "cluster minus projection deltas", 1400, 1000)
c3.Divide(2, 2)

c3.cd(1)
h_dx.SetLineWidth(2)
h_dx.Draw("HIST")

c3.cd(2)
h_dy.SetLineWidth(2)
h_dy.Draw("HIST")

c3.cd(3)
h_dz.SetLineWidth(2)
h_dz.Draw("HIST")

c3.cd(4)
h_dphi.SetLineWidth(2)
h_dphi.Draw("HIST")

c3.SaveAs(OUTPNG3)

c4 = ROOT.TCanvas("c4", "cluster minus projection spherical deltas", 1000, 800)
c4.Divide(1, 2)

c4.cd(1)
h_dr.SetLineWidth(2)
h_dr.Draw("HIST")

c4.cd(2)
h_dtheta.SetLineWidth(2)
h_dtheta.Draw("HIST")

c4.SaveAs(OUTPNG4)

# ---------------------------------
# Save ROOT file
# ---------------------------------
out = ROOT.TFile(OUTROOT, "RECREATE")
h_d3D.Write()
h_dxy.Write()
h_dx.Write()
h_dy.Write()
h_dz.Write()
h_dphi.Write()
h_dr.Write()
h_dtheta.Write()
h_nclus.Write()
h_nproj.Write()
c1.Write()
c2.Write()
c3.Write()
c4.Write()
out.Close()

print(f"\n✅ Wrote ROOT file: {OUTROOT}")
print(f"✅ Wrote image:     {OUTPNG1}")
print(f"✅ Wrote image:     {OUTPNG2}")
print(f"✅ Wrote image:     {OUTPNG3}")
print(f"✅ Wrote image:     {OUTPNG4}")
