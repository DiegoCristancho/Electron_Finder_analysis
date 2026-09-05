import ROOT
import math

ROOT.gROOT.SetBatch(True)

# -----------------------------
# Inputs
# -----------------------------
LISTFILE = "list_18x275_Q2_100_latest.txt"
OUTFILE  = "EcalBarrelClusters_xyz_TH3.root"

# If your clusters are stored as a separate TTree, set this to "EcalBarrelClusters".
# If they are a collection branch inside "events", set EVENT_TREE to "events" (default).
EVENT_TREE = "events"
CLUSTER_TREE_CANDIDATE = "EcalBarrelClusters"  # used only if present as a real TTree

# -----------------------------
# Build chain (events)
# -----------------------------
chain = ROOT.TChain(EVENT_TREE)
files = []
with open(LISTFILE) as f:
    for line in f:
        fn = line.strip()
        if fn:
            files.append(fn)
            chain.Add(fn)

print(f"✅ Added {chain.GetNtrees()} files, total {chain.GetEntries()} entries in '{EVENT_TREE}'")

# -----------------------------
# Helper: detect if a separate TTree exists
# (we check the first file only)
# -----------------------------
has_separate_cluster_tree = False
if len(files) > 0:
    tf = ROOT.TFile.Open(files[0])
    if tf and not tf.IsZombie():
        obj = tf.Get(CLUSTER_TREE_CANDIDATE)
        if obj and isinstance(obj, ROOT.TTree):
            has_separate_cluster_tree = True
    if tf:
        tf.Close()

print("ℹ️ Layout:",
      "separate TTree 'EcalBarrelClusters'" if has_separate_cluster_tree
      else "clusters as a collection branch inside 'events'")

# -----------------------------
# First pass: find min/max ranges (auto-range)
# -----------------------------
# Change MAX_SCAN if you want to scan fewer entries for speed.
MAX_SCAN = -1   # -1 = scan all entries
MARGIN   = 5.0  # add margin to min/max

xmin = ymin = zmin =  1e30
xmax = ymax = zmax = -1e30
n_points = 0

def update_ranges(x, y, z):
    global xmin, ymin, zmin, xmax, ymax, zmax, n_points
    if x < xmin: xmin = x
    if y < ymin: ymin = y
    if z < zmin: zmin = z
    if x > xmax: xmax = x
    if y > ymax: ymax = y
    if z > zmax: zmax = z
    n_points += 1

if has_separate_cluster_tree:
    # Build a chain for the cluster tree too (same file list)
    cchain = ROOT.TChain(CLUSTER_TREE_CANDIDATE)
    for fn in files:
        cchain.Add(fn)

    nentries = cchain.GetEntries() if MAX_SCAN < 0 else min(MAX_SCAN, cchain.GetEntries())
    for i in range(nentries):
        cchain.GetEntry(i)

        # Two common possibilities:
        # (A) branches named position.x / position.y / position.z
        # (B) a "position" object with x,y,z members
        try:
            x = float(getattr(cchain, "position.x"))
            y = float(getattr(cchain, "position.y"))
            z = float(getattr(cchain, "position.z"))
            update_ranges(x, y, z)
            continue
        except Exception:
            pass

        try:
            pos = getattr(cchain, "position")
            x = float(pos.x); y = float(pos.y); z = float(pos.z)
            update_ranges(x, y, z)
            continue
        except Exception:
            pass

else:
    # clusters are a collection inside event tree
    nentries = chain.GetEntries() if MAX_SCAN < 0 else min(MAX_SCAN, chain.GetEntries())
    for ievt, event in enumerate(chain):
        if ievt >= nentries:
            break

        # event.EcalBarrelClusters should be an iterable collection
        if not hasattr(event, "EcalBarrelClusters"):
            continue

        for cl in event.EcalBarrelClusters:
            # From your screenshot: cl.position.x, cl.position.y, cl.position.z
            x = float(cl.position.x)
            y = float(cl.position.y)
            z = float(cl.position.z)
            update_ranges(x, y, z)

if n_points == 0:
    raise RuntimeError("❌ No (x,y,z) points found. Check tree/branch names and layout.")

# Add margins
xmin -= MARGIN; ymin -= MARGIN; zmin -= MARGIN
xmax += MARGIN; ymax += MARGIN; zmax += MARGIN

print(f"✅ Found {n_points} cluster positions")
print(f"   x: [{xmin:.3f}, {xmax:.3f}]")
print(f"   y: [{ymin:.3f}, {ymax:.3f}]")
print(f"   z: [{zmin:.3f}, {zmax:.3f}]")

# -----------------------------
# Create TH3F (choose binning)
# -----------------------------
# Adjust bins to taste (more bins = more memory)
NBX, NBY, NBZ = 200, 200, 200

h3 = ROOT.TH3F(
    "h3_EcalBarrelClusters_xyz",
    "EcalBarrelClusters position; x; y; z",
    NBX, xmin, xmax,
    NBY, ymin, ymax,
    NBZ, zmin, zmax
)

# -----------------------------
# Second pass: fill histogram
# -----------------------------
if has_separate_cluster_tree:
    cchain = ROOT.TChain(CLUSTER_TREE_CANDIDATE)
    for fn in files:
        cchain.Add(fn)

    for i in range(cchain.GetEntries()):
        cchain.GetEntry(i)

        # Try same options as before
        filled = False
        try:
            x = float(getattr(cchain, "position.x"))
            y = float(getattr(cchain, "position.y"))
            z = float(getattr(cchain, "position.z"))
            h3.Fill(x, y, z)
            filled = True
        except Exception:
            pass

        if not filled:
            try:
                pos = getattr(cchain, "position")
                h3.Fill(float(pos.x), float(pos.y), float(pos.z))
            except Exception:
                pass

else:
    for event in chain:
        if not hasattr(event, "EcalBarrelClusters"):
            continue
        for cl in event.EcalBarrelClusters:
            h3.Fill(float(cl.position.x), float(cl.position.y), float(cl.position.z))

# -----------------------------
# Save output
# -----------------------------
out = ROOT.TFile(OUTFILE, "RECREATE")
h3.Write()
out.Close()

print(f"✅ Wrote TH3F to: {OUTFILE}")
print("   Histogram name: h3_EcalBarrelClusters_xyz")
