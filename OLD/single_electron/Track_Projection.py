import ROOT
import math

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# ---------------------------------
# Input / output
# ---------------------------------
LISTFILE = "single_electron.txt"
OUTROOT  = "track_vs_ecalclusters_first10_3Dmatch.root"
OUTPNG   = "track_vs_ecalclusters_first10_3Dmatch_xy.png"

EVENT_TREE = "events"
MAX_EVENTS = 10

# Based on what you found
ECAL_SYSTEM = 101
HCAL_SYSTEM = 111

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
# Read first 10 events
# ---------------------------------
events_data = []

global_xmin =  1e30
global_xmax = -1e30
global_ymin =  1e30
global_ymax = -1e30

for ievt, event in enumerate(chain):
    if ievt >= MAX_EVENTS:
        break

    if not hasattr(event, "CalorimeterTrackProjections"):
        continue
    if not hasattr(event, "_CalorimeterTrackProjections_points"):
        continue
    if not hasattr(event, "EcalBarrelClusters"):
        continue

    projs    = event.CalorimeterTrackProjections
    points   = event._CalorimeterTrackProjections_points
    clusters = event.EcalBarrelClusters

    if len(projs) < 1:
        continue

    proj = projs[0]
    b = int(proj.points_begin)
    e = int(proj.points_end)

    track_pts = []
    ecal_proj_pt = None

    for j in range(b, e):
        pt = points[j]
        x = float(pt.position.x)
        y = float(pt.position.y)
        z = float(pt.position.z)
        sysid = int(pt.system)

        track_pts.append((x, y, z, sysid))

        if sysid == ECAL_SYSTEM:
            ecal_proj_pt = (x, y, z)

        global_xmin = min(global_xmin, x)
        global_xmax = max(global_xmax, x)
        global_ymin = min(global_ymin, y)
        global_ymax = max(global_ymax, y)

    cluster_pts = []
    for cl in clusters:
        x = float(cl.position.x)
        y = float(cl.position.y)
        z = float(cl.position.z)
        cluster_pts.append((x, y, z))

        global_xmin = min(global_xmin, x)
        global_xmax = max(global_xmax, x)
        global_ymin = min(global_ymin, y)
        global_ymax = max(global_ymax, y)

    nearest_dist_3d = None
    nearest_dist_xy = None
    nearest_cluster = None

    if ecal_proj_pt is not None and len(cluster_pts) > 0:
        px, py, pz = ecal_proj_pt

        best_d3 = 1e30
        best_dxy = None
        best_cl = None

        for (cx, cy, cz) in cluster_pts:
            dxy = math.sqrt((cx - px)**2 + (cy - py)**2)
            d3  = math.sqrt((cx - px)**2 + (cy - py)**2 + (cz - pz)**2)

            if d3 < best_d3:
                best_d3 = d3
                best_dxy = dxy
                best_cl = (cx, cy, cz)

        nearest_dist_3d = best_d3
        nearest_dist_xy = best_dxy
        nearest_cluster = best_cl

    events_data.append({
        "event": ievt,
        "track_pts": track_pts,
        "ecal_proj_pt": ecal_proj_pt,
        "cluster_pts": cluster_pts,
        "nearest_dist_3d": nearest_dist_3d,
        "nearest_dist_xy": nearest_dist_xy,
        "nearest_cluster": nearest_cluster
    })

if len(events_data) == 0:
    raise RuntimeError("❌ No usable events found.")

margin = 20.0
global_xmin -= margin
global_xmax += margin
global_ymin -= margin
global_ymax += margin

# ---------------------------------
# Canvas
# ---------------------------------
c = ROOT.TCanvas("c", "Track projections vs ECal clusters", 2200, 900)
c.Divide(5, 2)

saved_objects = []

for ipad, ev in enumerate(events_data, start=1):
    c.cd(ipad)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTopMargin(0.10)

    ievt = ev["event"]
    track_pts = ev["track_pts"]
    cluster_pts = ev["cluster_pts"]
    nearest_cluster = ev["nearest_cluster"]
    nearest_dist_3d = ev["nearest_dist_3d"]
    nearest_dist_xy = ev["nearest_dist_xy"]

    frame = ROOT.TH2F(
        f"frame_evt{ievt}",
        f"Event {ievt};x [mm];y [mm]",
        100, global_xmin, global_xmax,
        100, global_ymin, global_ymax
    )
    frame.SetStats(0)
    frame.Draw()
    saved_objects.append(frame)

    g_ecal = ROOT.TGraph()
    g_hcal = ROOT.TGraph()
    g_line = ROOT.TGraph()
    g_clus = ROOT.TGraph()
    g_best = ROOT.TGraph()

    for i, (x, y, z, sysid) in enumerate(track_pts):
        g_line.SetPoint(i, x, y)

        if sysid == ECAL_SYSTEM:
            g_ecal.SetPoint(g_ecal.GetN(), x, y)
        elif sysid == HCAL_SYSTEM:
            g_hcal.SetPoint(g_hcal.GetN(), x, y)

    for i, (x, y, z) in enumerate(cluster_pts):
        g_clus.SetPoint(i, x, y)

    if nearest_cluster is not None:
        g_best.SetPoint(0, nearest_cluster[0], nearest_cluster[1])

    g_line.SetLineColor(ROOT.kBlue + 1)
    g_line.SetLineWidth(2)

    g_ecal.SetMarkerStyle(20)
    g_ecal.SetMarkerSize(1.2)
    g_ecal.SetMarkerColor(ROOT.kGreen + 2)

    g_hcal.SetMarkerStyle(20)
    g_hcal.SetMarkerSize(1.0)
    g_hcal.SetMarkerColor(ROOT.kRed)

    g_clus.SetMarkerStyle(24)
    g_clus.SetMarkerSize(1.0)
    g_clus.SetMarkerColor(ROOT.kBlack)

    g_best.SetMarkerStyle(29)
    g_best.SetMarkerSize(1.6)
    g_best.SetMarkerColor(ROOT.kMagenta + 1)

    g_clus.Draw("P SAME")
    g_line.Draw("L SAME")
    g_ecal.Draw("P SAME")
    g_hcal.Draw("P SAME")
    if nearest_cluster is not None:
        g_best.Draw("P SAME")

    leg = ROOT.TLegend(0.14, 0.70, 0.56, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.AddEntry(g_clus, "ECal clusters", "p")
    leg.AddEntry(g_ecal, "ECal projection point", "p")
    leg.AddEntry(g_hcal, "HCal projection points", "p")
    leg.AddEntry(g_line, "Track projection path", "l")
    if nearest_cluster is not None:
        leg.AddEntry(g_best, "Nearest cluster (3D)", "p")
    leg.Draw()

    txt = ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextSize(0.030)

    if nearest_dist_3d is not None:
        txt.DrawLatex(0.14, 0.63, f"d_{{3D}} = {nearest_dist_3d:.2f} mm")
        txt.DrawLatex(0.14, 0.58, f"d_{{xy}} = {nearest_dist_xy:.2f} mm")
        print(f"Event {ievt}: nearest cluster d3D = {nearest_dist_3d:.3f} mm, dxy = {nearest_dist_xy:.3f} mm")
    else:
        txt.DrawLatex(0.14, 0.63, "No ECAL projection point or no ECAL cluster")
        print(f"Event {ievt}: no ECAL projection point or no ECAL cluster")

    saved_objects.extend([g_ecal, g_hcal, g_line, g_clus, g_best, leg, txt])

c.SaveAs(OUTPNG)

out = ROOT.TFile(OUTROOT, "RECREATE")
c.Write()
for obj in saved_objects:
    if obj:
        obj.Write()
out.Close()

print(f"✅ Wrote image: {OUTPNG}")
print(f"✅ Wrote ROOT file: {OUTROOT}")
