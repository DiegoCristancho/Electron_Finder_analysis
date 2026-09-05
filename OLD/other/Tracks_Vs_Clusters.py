import ROOT

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
INPUTROOT = "pythia8NCDIS_26_718x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_1.0008.eicrecon.edm4eic.root"

EVENT_TREE = "events"
TARGET_EVENT = 1
OUTROOT = "one_event_track_projection_points.root"

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

if not hasattr(target, "CalorimeterTrackProjections"):
    raise RuntimeError("No CalorimeterTrackProjections branch")
if not hasattr(target, "_CalorimeterTrackProjections_track"):
    raise RuntimeError("No _CalorimeterTrackProjections_track branch")
if not hasattr(target, "_CalorimeterTrackProjections_points"):
    raise RuntimeError("No _CalorimeterTrackProjections_points branch")

projs  = target.CalorimeterTrackProjections
tracks = target._CalorimeterTrackProjections_track
points = target._CalorimeterTrackProjections_points

print(f"\n========== Event {TARGET_EVENT} ==========")
print(f"n projections = {len(projs)}")
print(f"n track refs  = {len(tracks)}")
print(f"n total points = {len(points)}\n")

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

track_data = []

for i in range(len(projs)):
    proj = projs[i]
    trk  = tracks[i]

    b = int(proj.points_begin)
    e = int(proj.points_end)

    pts = []
    print(f"Projection {i}")
    print(f"  track.index        = {int(trk.index)}")
    print(f"  track.collectionID = {int(trk.collectionID)}")
    print(f"  points_begin       = {b}")
    print(f"  points_end         = {e}")
    print(f"  n points           = {e-b}")

    for j in range(b, e):
        pt = points[j]
        x = float(pt.position.x)
        y = float(pt.position.y)
        z = float(pt.position.z)
        sysid = int(pt.system)
        surf  = int(pt.surface)

        pts.append((x, y, z, sysid, surf))

        xmin = min(xmin, x)
        xmax = max(xmax, x)
        ymin = min(ymin, y)
        ymax = max(ymax, y)
        zmin = min(zmin, z)
        zmax = max(zmax, z)

        print(
            f"    point {j}: "
            f"system={sysid}, surface={surf}, "
            f"x={x:.3f}, y={y:.3f}, z={z:.3f}"
        )

    print()
    track_data.append({
        "proj_index": i,
        "track_index": int(trk.index),
        "collectionID": int(trk.collectionID),
        "points_begin": b,
        "points_end": e,
        "pts": pts
    })

margin = 20.0
xmin -= margin
xmax += margin
ymin -= margin
ymax += margin
zmin -= margin
zmax += margin

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

leg3d = ROOT.TLegend(0.72, 0.65, 0.92, 0.90)
leg3d.SetBorderSize(0)
leg3d.SetFillStyle(0)

for i, td in enumerate(track_data):
    pts = td["pts"]
    color = colors[i % len(colors)]

    line = ROOT.TPolyLine3D(len(pts))
    mark = ROOT.TPolyMarker3D(len(pts))

    for k, (x, y, z, sysid, surf) in enumerate(pts):
        line.SetPoint(k, x, y, z)
        mark.SetPoint(k, x, y, z)

    line.SetLineColor(color)
    line.SetLineWidth(2)

    mark.SetMarkerColor(color)
    mark.SetMarkerStyle(20)
    mark.SetMarkerSize(1.0)

    line.Draw("same")
    mark.Draw("same")

    saved.extend([line, mark])
    leg3d.AddEntry(mark, f"track {td['track_index']}", "p")

leg3d.Draw()
saved.append(leg3d)

c3d.SaveAs("one_event_track_projection_points_3D.png")

# ---------------------------------
# X-Y canvas
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

legxy = ROOT.TLegend(0.72, 0.65, 0.92, 0.90)
legxy.SetBorderSize(0)
legxy.SetFillStyle(0)

for i, td in enumerate(track_data):
    pts = td["pts"]
    color = colors[i % len(colors)]

    g = ROOT.TGraph(len(pts))
    for k, (x, y, z, sysid, surf) in enumerate(pts):
        g.SetPoint(k, x, y)

    g.SetLineColor(color)
    g.SetMarkerColor(color)
    g.SetLineWidth(2)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.1)

    g.Draw("LP SAME")
    saved.append(g)
    legxy.AddEntry(g, f"track {td['track_index']}", "lp")

legxy.Draw()
saved.append(legxy)

cxy.SaveAs("one_event_track_projection_points_XY.png")

# ---------------------------------
# Save ROOT file
# ---------------------------------
out = ROOT.TFile(OUTROOT, "RECREATE")
c3d.Write()
cxy.Write()
for obj in saved:
    obj.Write()
out.Close()

print("Saved:")
print("  one_event_track_projection_points_3D.png")
print("  one_event_track_projection_points_XY.png")
print(f"  {OUTROOT}")
