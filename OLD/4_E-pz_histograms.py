import ROOT

# -----------------------------
# -----------------------------
# Load chain of files
# -----------------------------
tree = ROOT.TChain("events")
with open("list_SIDIS_18x275_q2_0to1.txt") as f:
    for line in f:
        fname = line.strip()
        if fname:
            tree.Add(fname)

if tree.GetEntries() == 0:
    raise RuntimeError("Error: No events found in the chain!")

print(f"✅ Added {tree.GetNtrees()} files with total {tree.GetEntries()} events")

outfile = ROOT.TFile("SIDIS_analysis_merged.root", "RECREATE")

Ee = 18.0                 # electron beam energy [GeV]
twoEe = 2.0 * Ee          # = 36 GeV

# PDGs to study for fractions
pdg_list = [2112, 130]    # neutron, K0L

# -----------------------------
# Reduced tree (event-wise)
# -----------------------------
pdg_vec    = ROOT.std.vector('int')()
status_vec = ROOT.std.vector('int')()
px_vec     = ROOT.std.vector('double')()
py_vec     = ROOT.std.vector('double')()
pz_vec     = ROOT.std.vector('double')()
mass_vec   = ROOT.std.vector('double')()
energy_vec = ROOT.std.vector('double')()

reduced_tree = ROOT.TTree("MCParticlesReduced", "MCParticles (vectors per event)")
reduced_tree.Branch("PDG",    pdg_vec)
reduced_tree.Branch("status", status_vec)
reduced_tree.Branch("px",     px_vec)
reduced_tree.Branch("py",     py_vec)
reduced_tree.Branch("pz",     pz_vec)
reduced_tree.Branch("mass",   mass_vec)
reduced_tree.Branch("energy", energy_vec)

# -----------------------------
# Histograms you requested
# -----------------------------
# 1) Fractions (E-Pz)/36 by PDG (status==1, MCParticles)
h_frac = {
    2112: ROOT.TH1F("E_Pz_fraction_2112", "(E-Pz)/36 per event (PDG=2112, status==1)", 100, 0.0, 1.0),
    130:  ROOT.TH1F("E_Pz_fraction_130",  "(E-Pz)/36 per event (PDG=130,  status==1)", 100, 0.0, 1.0),
}

# 2) Total (E-Pz) per event for ReconstructedParticles
h_reco_total = ROOT.TH1F("E_Pz_total_Reco", "Total (E-Pz) per event (ReconstructedParticles)", 200, 0.0, 50.0)

# 3) Total (E-Pz) per event for MCParticles with status==1
h_mcstat1_total = ROOT.TH1F("E_Pz_total_MC_status1", "Total (E-Pz) per event (MCParticles, status==1); E - p_{z} [GeV]", 200, 0.0, 50.0)

# -----------------------------
# Event loop
# -----------------------------
for event in tree:
    # ---- Reduced tree (one fill per event) ----
    pdg_vec.clear(); status_vec.clear()
    px_vec.clear();  py_vec.clear(); pz_vec.clear()
    mass_vec.clear(); energy_vec.clear()

    for part in event.MCParticles:
        pdg_vec.push_back(part.PDG)
        status_vec.push_back(part.generatorStatus)
        px_vec.push_back(part.momentum.x)
        py_vec.push_back(part.momentum.y)
        pz_vec.push_back(part.momentum.z)
        mass_vec.push_back(part.mass)
        # True energy from MC
        e = (part.momentum.x**2 + part.momentum.y**2 + part.momentum.z**2 + part.mass**2) ** 0.5
        energy_vec.push_back(e)

    reduced_tree.Fill()  # <-- exactly once per event

    # ---- Analysis histograms ----

    # A) MCParticles: status==1 totals and fractions by PDG
    sums_by_pdg = {pdg: 0.0 for pdg in pdg_list}
    mc_status1_total = 0.0

    for part in event.MCParticles:
        if part.generatorStatus != 1:
            continue
        e = (part.momentum.x**2 + part.momentum.y**2 + part.momentum.z**2 + part.mass**2) ** 0.5
        epz = e - part.momentum.z
        mc_status1_total += epz
        if part.PDG in sums_by_pdg:
            sums_by_pdg[part.PDG] += epz

    # fill total (status==1)
    h_mcstat1_total.Fill(mc_status1_total)

    # fill per-PDG fractions
    for pdg, s in sums_by_pdg.items():
        h_frac[pdg].Fill(s / twoEe)

    # B) ReconstructedParticles: total (E-Pz) per event
    reco_total = 0.0
    for part in event.ReconstructedParticles:
        epz = part.energy - part.momentum.z
        reco_total += epz
    h_reco_total.Fill(reco_total)

# -----------------------------
# Write out
# -----------------------------
outfile.cd()
reduced_tree.Write()
h_reco_total.Write()
h_mcstat1_total.Write()
for h in h_frac.values():
    h.Write()

outfile.Close()
print("✅ Wrote event-wise MCParticlesReduced and requested histograms to SIDS_analysis_merged.root")

