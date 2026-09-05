import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(1110)

# ---------------------------------
# Inputs
# ---------------------------------
LISTFILE = "list_18x275_Q2_100_26_02.txt"   # change if needed
EVENT_TREE = "events"

OUTROOT = "parent_status_first_status1_electron.root"
OUTPNG  = "parent_status_first_status1_electron.png"

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
# Histogram
# ---------------------------------
h_parent_status = ROOT.TH1F(
    "h_parent_status",
    "Parent status of first MC electron with PDG=11 and generatorStatus=1;Parent generator status;Events",
    60, -0.5, 59.5
)

# Optional parent PDG histogram too
h_parent_pdg = ROOT.TH1F(
    "h_parent_pdg",
    "Parent PDG of first MC electron with PDG=11 and generatorStatus=1;Parent PDG;Events",
    100, -50, 50
)

# ---------------------------------
# Counters / bookkeeping
# ---------------------------------
n_total = 0
n_found_electron = 0
n_with_parent = 0
n_no_parent = 0
n_bad_parent_index = 0

records = []

# ---------------------------------
# Event loop
# ---------------------------------
for ievt, event in enumerate(chain):
    n_total += 1

    if not hasattr(event, "MCParticles"):
        continue
    if not hasattr(event, "_MCParticles_parents"):
        continue

    mcparts = event.MCParticles
    parent_refs = event._MCParticles_parents

    electron_index = None

    # Find first electron with generatorStatus = 1 and PDG = 11
    for i, part in enumerate(mcparts):
        if int(part.generatorStatus) == 1 and int(part.PDG) == 11:
            electron_index = i
            break

    if electron_index is None:
        continue

    n_found_electron += 1
    ele = mcparts[electron_index]

    pb = int(ele.parents_begin)
    pe = int(ele.parents_end)

    if pe <= pb:
        n_no_parent += 1
        records.append({
            "event": ievt,
            "electron_index": electron_index,
            "electron_pdg": int(ele.PDG),
            "electron_status": int(ele.generatorStatus),
            "parent_index": None,
            "parent_pdg": None,
            "parent_status": None,
        })
        continue

    # Use the first parent
    pref = parent_refs[pb]
    parent_index = int(pref.index)

    if parent_index < 0 or parent_index >= len(mcparts):
        n_bad_parent_index += 1
        records.append({
            "event": ievt,
            "electron_index": electron_index,
            "electron_pdg": int(ele.PDG),
            "electron_status": int(ele.generatorStatus),
            "parent_index": parent_index,
            "parent_pdg": None,
            "parent_status": None,
        })
        continue

    parent = mcparts[parent_index]
    parent_pdg = int(parent.PDG)
    parent_status = int(parent.generatorStatus)

    h_parent_status.Fill(parent_status)
    h_parent_pdg.Fill(parent_pdg)
    n_with_parent += 1

    records.append({
        "event": ievt,
        "electron_index": electron_index,
        "electron_pdg": int(ele.PDG),
        "electron_status": int(ele.generatorStatus),
        "parent_index": parent_index,
        "parent_pdg": parent_pdg,
        "parent_status": parent_status,
    })

# ---------------------------------
# Print summary
# ---------------------------------
print("\n========== Summary ==========")
print(f"Total events                                : {n_total}")
print(f"Events with first electron (PDG=11, status=1): {n_found_electron}")
print(f"Of those, with parent                       : {n_with_parent}")
print(f"Of those, no parent                         : {n_no_parent}")
print(f"Of those, invalid parent index              : {n_bad_parent_index}")

if n_with_parent > 0:
    print(f"Mean parent status                          : {h_parent_status.GetMean():.3f}")

# Count how many are exactly status 23
bin23 = h_parent_status.FindBin(23)
n_status23 = h_parent_status.GetBinContent(bin23)
print(f"Parent status == 23 count                   : {int(n_status23)}")

if n_with_parent > 0:
    frac23 = n_status23 / n_with_parent
    print(f"Fraction with parent status == 23           : {frac23:.6f}")

# Print non-23 cases
print("\nEvents where parent status is not 23:")
found_non23 = False
for rec in records:
    if rec["parent_status"] is not None and rec["parent_status"] != 23:
        found_non23 = True
        print(
            f"  event={rec['event']}, "
            f"electron_idx={rec['electron_index']}, "
            f"parent_idx={rec['parent_index']}, "
            f"parent_pdg={rec['parent_pdg']}, "
            f"parent_status={rec['parent_status']}"
        )

if not found_non23:
    print("  None found.")

# ---------------------------------
# Draw histograms
# ---------------------------------
c1 = ROOT.TCanvas("c1", "parent status", 900, 700)
h_parent_status.SetLineWidth(2)
h_parent_status.Draw("HIST")
c1.SaveAs(OUTPNG)

c2 = ROOT.TCanvas("c2", "parent pdg", 900, 700)
h_parent_pdg.SetLineWidth(2)
h_parent_pdg.Draw("HIST")
c2.SaveAs("parent_pdg_first_status1_electron.png")

# ---------------------------------
# Save ROOT file
# ---------------------------------
out = ROOT.TFile(OUTROOT, "RECREATE")
h_parent_status.Write()
h_parent_pdg.Write()
c1.Write()
c2.Write()
out.Close()

print(f"\n✅ Wrote ROOT file: {OUTROOT}")
print(f"✅ Wrote image:     {OUTPNG}")
print(f"✅ Wrote image:     parent_pdg_first_status1_electron.png")
