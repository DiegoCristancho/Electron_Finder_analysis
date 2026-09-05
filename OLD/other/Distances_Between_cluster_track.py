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


# ============================================================
# Helper functions
# ============================================================

def particle_energy(obj):
    px = float(obj.momentum.x)
    py = float(obj.momentum.y)
    pz = float(obj.momentum.z)

    m = float(obj.mass) if hasattr(obj, "mass") else 0.0

    return math.sqrt(
        px * px +
        py * py +
        pz * pz +
        m * m
    )


def particle_p(obj):
    px = float(obj.momentum.x)
    py = float(obj.momentum.y)
    pz = float(obj.momentum.z)

    return math.sqrt(
        px * px +
        py * py +
        pz * pz
    )


# ------------------------------------------------------------
# Position-based spherical coordinates
#
# IMPORTANT:
# These functions use ONLY x, y, z.
#
# They do NOT use reconstructed eta, phi, or theta.
# ------------------------------------------------------------

def position_r(x, y, z):
    """
    3D spherical radius:
        r = sqrt(x^2 + y^2 + z^2)
    """
    return math.sqrt(
        x * x +
        y * y +
        z * z
    )


def position_theta(x, y, z):
    """
    Polar angle calculated from position:
        theta = acos(z/r)
    """
    r = position_r(x, y, z)

    if r <= 0:
        return None

    cos_theta = z / r

    # numerical protection
    cos_theta = max(-1.0, min(1.0, cos_theta))

    return math.acos(cos_theta)


def position_eta(x, y, z):
    """
    Pseudorapidity calculated ONLY from the position:
        eta = -ln(tan(theta/2))
    """
    theta = position_theta(x, y, z)

    if theta is None:
        return None

    tan_half = math.tan(theta / 2.0)

    if tan_half <= 0:
        return None

    return -math.log(tan_half)


# ------------------------------------------------------------
# Parent information
# ------------------------------------------------------------

def get_first_parent_status(part, mc_particles, parent_refs):

    pb = int(part.parents_begin)
    pe = int(part.parents_end)

    if pe <= pb:
        return None

    pref = parent_refs[pb]

    parent_idx = int(pref.index)

    if parent_idx < 0 or parent_idx >= len(mc_particles):
        return None

    return int(
        mc_particles[parent_idx].generatorStatus
    )


def get_first_parent_pdg(part, mc_particles, parent_refs):

    pb = int(part.parents_begin)
    pe = int(part.parents_end)

    if pe <= pb:
        return None

    pref = parent_refs[pb]

    parent_idx = int(pref.index)

    if parent_idx < 0 or parent_idx >= len(mc_particles):
        return None

    return int(
        mc_particles[parent_idx].PDG
    )


# ------------------------------------------------------------
# Find scattered electron
# ------------------------------------------------------------

def find_first_scattered_electron_mc_index(mc_particles, parent_refs):
    """
    Scattered-electron definition:

        generatorStatus == 1
        PDG == 11
        parent_status != 2

    Take the first electron satisfying these conditions.
    """

    for i, part in enumerate(mc_particles):

        if int(part.generatorStatus) != 1:
            continue

        if int(part.PDG) != 11:
            continue

        parent_status = get_first_parent_status(
            part,
            mc_particles,
            parent_refs
        )

        if parent_status != 2:
            return i

    return None


# ------------------------------------------------------------
# MC -> reconstructed particle
# ------------------------------------------------------------

def find_reco_index_for_mc_index(
    mc_index,
    assoc_rec,
    assoc_sim
):

    if mc_index is None:
        return None

    for ia in range(len(assoc_sim)):

        sim_ref = assoc_sim[ia]

        if int(sim_ref.index) == int(mc_index):

            return int(
                assoc_rec[ia].index
            )

    return None


# ------------------------------------------------------------
# Reconstructed particle -> track
# ------------------------------------------------------------

def find_track_ref_for_reco_index(
    reco_index,
    reco_particles,
    reco_particle_tracks
):

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

    return (
        int(a.collectionID) == int(b.collectionID)
        and
        int(a.index) == int(b.index)
    )


# ============================================================
# Build chain
# ============================================================

chain = ROOT.TChain(EVENT_TREE)

if USE_LISTFILE:

    with open(LISTFILE) as f:

        for line in f:

            fn = line.strip()

            if fn:
                chain.Add(fn)

    print(
        f"✅ Added {chain.GetNtrees()} files from txt list, "
        f"total {chain.GetEntries()} events"
    )

else:

    chain.Add(INPUTROOT)

    print(
        f"✅ Added {chain.GetNtrees()} file, "
        f"total {chain.GetEntries()} events"
    )


# ============================================================
# Histograms
# ============================================================

# ------------------------------------------------------------
# d3D without E/p
# ------------------------------------------------------------

h_d3D_noEoP = ROOT.TH1F(
    "h_d3D_noEoP",
    "Nearest ECal cluster distance for scattered electron "
    "(no E/p cut);d_{3D} [mm];Events",
    100,
    0,
    200
)


# ------------------------------------------------------------
# d3D with E/p
# ------------------------------------------------------------

h_d3D_withEoP = ROOT.TH1F(
    "h_d3D_withEoP",
    f"Nearest ECal cluster distance for scattered electron "
    f"({EOP_MIN:.1f} < E/p < {EOP_MAX:.1f});"
    f"d_{{3D}} [mm];Events",
    100,
    0,
    200
)


# ============================================================
# NEW: position-based coordinate differences
# ============================================================

# ------------------------------------------------------------
# WITHOUT E/p
# ------------------------------------------------------------

h_dr_noEoP = ROOT.TH1F(
    "h_dr_noEoP",
    "Cluster - projection #Delta r (no E/p cut);"
    "r_{cluster} - r_{proj} [mm];Events",
    120,
    -150,
    150
)

h_dtheta_noEoP = ROOT.TH1F(
    "h_dtheta_noEoP",
    "Cluster - projection #Delta#theta (no E/p cut);"
    "#theta_{cluster} - #theta_{proj} [rad];Events",
    120,
    -0.2,
    0.2
)

h_deta_noEoP = ROOT.TH1F(
    "h_deta_noEoP",
    "Cluster - projection #Delta#eta (no E/p cut);"
    "#eta_{cluster} - #eta_{proj};Events",
    120,
    -0.5,
    0.5
)


# ------------------------------------------------------------
# WITH E/p
# ------------------------------------------------------------

h_dr_withEoP = ROOT.TH1F(
    "h_dr_withEoP",
    f"Cluster - projection #Delta r "
    f"({EOP_MIN:.1f} < E/p < {EOP_MAX:.1f});"
    f"r_{{cluster}} - r_{{proj}} [mm];Events",
    120,
    -150,
    150
)

h_dtheta_withEoP = ROOT.TH1F(
    "h_dtheta_withEoP",
    f"Cluster - projection #Delta#theta "
    f"({EOP_MIN:.1f} < E/p < {EOP_MAX:.1f});"
    f"#theta_{{cluster}} - #theta_{{proj}} [rad];Events",
    120,
    -0.2,
    0.2
)

h_deta_withEoP = ROOT.TH1F(
    "h_deta_withEoP",
    f"Cluster - projection #Delta#eta "
    f"({EOP_MIN:.1f} < E/p < {EOP_MAX:.1f});"
    f"#eta_{{cluster}} - #eta_{{proj}};Events",
    120,
    -0.5,
    0.5
)


# ------------------------------------------------------------
# E/p diagnostics
# ------------------------------------------------------------

h_eop_all = ROOT.TH1F(
    "h_eop_all",
    "E/p of all ECal clusters considered for scattered electron;"
    "E/p;Clusters",
    120,
    0,
    3
)

h_eop_selected = ROOT.TH1F(
    "h_eop_selected",
    f"E/p of selected clusters passing "
    f"{EOP_MIN:.1f} < E/p < {EOP_MAX:.1f};"
    f"E/p;Clusters",
    120,
    0,
    3
)


# ------------------------------------------------------------
# Fraction vs d3D cut
# ------------------------------------------------------------

NBINS_CUT = 100
DCUT_MAX = 200.0


h_total = ROOT.TH1F(
    "h_total",
    "Total scattered electrons with ECAL projection point;"
    "distance cut [mm];Events",
    NBINS_CUT,
    0,
    DCUT_MAX
)


h_pass_d3D_noEoP = ROOT.TH1F(
    "h_pass_d3D_noEoP",
    "Scattered electrons with nearest cluster within d_{3D} cut "
    "(no E/p);distance cut [mm];Events",
    NBINS_CUT,
    0,
    DCUT_MAX
)


h_pass_d3D_withEoP = ROOT.TH1F(
    "h_pass_d3D_withEoP",
    f"Scattered electrons with nearest cluster within d_{{3D}} cut "
    f"({EOP_MIN:.1f} < E/p < {EOP_MAX:.1f});"
    f"distance cut [mm];Events",
    NBINS_CUT,
    0,
    DCUT_MAX
)


# ------------------------------------------------------------
# Failed scattered-electron energies
# ------------------------------------------------------------

h_fail_recoE = ROOT.TH1F(
    "h_fail_recoE",
    f"Reco energy of scattered electrons failing "
    f"d_{{3D}} < {CHOSEN_CUT_MM:.0f} mm;"
    f"Reco electron energy [GeV];Events",
    100,
    0,
    30
)


h_fail_mcE = ROOT.TH1F(
    "h_fail_mcE",
    f"MC energy of scattered electrons failing "
    f"d_{{3D}} < {CHOSEN_CUT_MM:.0f} mm;"
    f"MC electron energy [GeV];Events",
    100,
    0,
    30
)


# ------------------------------------------------------------
# Parent-status histograms
# ------------------------------------------------------------

h_parent_status_all = ROOT.TH1F(
    "h_parent_status_all",
    "Parent generator status of all electrons with status=1 "
    "and PDG=11;"
    "Parent generator status;Electrons",
    60,
    -0.5,
    59.5
)


h_parent_status_selected = ROOT.TH1F(
    "h_parent_status_selected",
    "Parent generator status of selected scattered electrons "
    "(parent status != 2);"
    "Parent generator status;Events",
    60,
    -0.5,
    59.5
)


h_parent_pdg_all = ROOT.TH1F(
    "h_parent_pdg_all",
    "Parent PDG of all electrons with status=1 and PDG=11;"
    "Parent PDG;Electrons",
    200,
    -100,
    100
)


# ============================================================
# Counters / bookkeeping
# ============================================================

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


# ============================================================
# Event loop
# ============================================================

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

    if any(
        not hasattr(event, name)
        for name in required
    ):
        continue


    mc_particles = event.MCParticles

    parent_refs = event._MCParticles_parents

    projs = event.CalorimeterTrackProjections

    proj_tracks = event._CalorimeterTrackProjections_track

    points = event._CalorimeterTrackProjections_points

    clusters = event.EcalBarrelClusters

    reco_particles = event.ReconstructedParticles

    reco_particle_tracks = event._ReconstructedParticles_tracks

    assoc_rec = event._ReconstructedParticleAssociations_rec

    assoc_sim = event._ReconstructedParticleAssociations_sim


    # ========================================================
    # Parent-status histogram BEFORE parent_status != 2 cut
    # ========================================================

    for part in mc_particles:

        if (
            int(part.generatorStatus) == 1
            and
            int(part.PDG) == 11
        ):

            n_all_status1_electrons += 1

            parent_status = get_first_parent_status(
                part,
                mc_particles,
                parent_refs
            )

            parent_pdg = get_first_parent_pdg(
                part,
                mc_particles,
                parent_refs
            )

            if parent_status is not None:
                h_parent_status_all.Fill(
                    parent_status
                )

            if parent_pdg is not None:
                h_parent_pdg_all.Fill(
                    parent_pdg
                )


    # ========================================================
    # Step 1:
    # Find scattered electron in MC
    #
    # generatorStatus = 1
    # PDG = 11
    # parent_status != 2
    # ========================================================

    scattered_mc_index = (
        find_first_scattered_electron_mc_index(
            mc_particles,
            parent_refs
        )
    )

    if scattered_mc_index is None:
        continue


    n_events_with_scattered_electron += 1

    scattered_mc = mc_particles[
        scattered_mc_index
    ]

    mc_energy = particle_energy(
        scattered_mc
    )


    parent_status_selected = (
        get_first_parent_status(
            scattered_mc,
            mc_particles,
            parent_refs
        )
    )

    if parent_status_selected is not None:

        h_parent_status_selected.Fill(
            parent_status_selected
        )


    # ========================================================
    # Step 2:
    # MC scattered electron -> ReconstructedParticle
    # ========================================================

    reco_index = find_reco_index_for_mc_index(
        scattered_mc_index,
        assoc_rec,
        assoc_sim
    )


    if (
        reco_index is None
        or
        reco_index < 0
        or
        reco_index >= len(reco_particles)
    ):

        failed_records.append({
            "event": ievt,
            "reason": "no_reco_match",
            "nearest_d3_noEoP": None,
            "nearest_d3_withEoP": None,
            "reco_energy": None,
            "mc_energy": mc_energy,
            "parent_status": parent_status_selected,
        })

        h_fail_mcE.Fill(
            mc_energy
        )

        continue


    n_scattered_with_reco += 1


    reco_obj = reco_particles[
        reco_index
    ]

    reco_energy = particle_energy(
        reco_obj
    )

    reco_p = particle_p(
        reco_obj
    )


    # ========================================================
    # Step 3:
    # ReconstructedParticle -> reconstructed track
    # ========================================================

    target_track_ref = (
        find_track_ref_for_reco_index(
            reco_index,
            reco_particles,
            reco_particle_tracks
        )
    )


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

        h_fail_recoE.Fill(
            reco_energy
        )

        h_fail_mcE.Fill(
            mc_energy
        )

        continue


    n_scattered_with_track += 1


    # ========================================================
    # Step 4:
    # Find corresponding calorimeter projection
    # ========================================================

    matched_proj_index = None


    for i in range(len(projs)):

        if same_track_ref(
            proj_tracks[i],
            target_track_ref
        ):

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

        h_fail_recoE.Fill(
            reco_energy
        )

        h_fail_mcE.Fill(
            mc_energy
        )

        continue


    n_scattered_with_projection += 1


    # ========================================================
    # Step 5:
    # Find ECAL projection point (system == 101)
    # ========================================================

    proj = projs[
        matched_proj_index
    ]


    b = int(
        proj.points_begin
    )

    e = int(
        proj.points_end
    )


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

        h_fail_recoE.Fill(
            reco_energy
        )

        h_fail_mcE.Fill(
            mc_energy
        )


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


    # ========================================================
    # Calculate projection r, theta, eta FROM POSITION ONLY
    # ========================================================

    proj_r = position_r(
        px,
        py,
        pz
    )

    proj_theta = position_theta(
        px,
        py,
        pz
    )

    proj_eta = position_eta(
        px,
        py,
        pz
    )


    # ========================================================
    # Step 6a:
    # Nearest ECAL cluster WITHOUT E/p cut
    # ========================================================

    best_d3_noEoP = None

    best_dxy_noEoP = None

    best_cluster_noEoP = None


    if len(clusters) > 0:

        best_d3_tmp = 1e30

        best_dxy_tmp = 1e30

        found_noEoP = False


        for cl in clusters:

            cx = float(
                cl.position.x
            )

            cy = float(
                cl.position.y
            )

            cz = float(
                cl.position.z
            )


            dxy = math.sqrt(
                (cx - px) ** 2 +
                (cy - py) ** 2
            )


            d3 = math.sqrt(
                (cx - px) ** 2 +
                (cy - py) ** 2 +
                (cz - pz) ** 2
            )


            if d3 < best_d3_tmp:

                best_d3_tmp = d3

                best_dxy_tmp = dxy

                best_cluster_noEoP = (
                    cx,
                    cy,
                    cz
                )

                found_noEoP = True


        if found_noEoP:

            best_d3_noEoP = (
                best_d3_tmp
            )

            best_dxy_noEoP = (
                best_dxy_tmp
            )


            n_scattered_with_cluster_noEoP += 1


            h_d3D_noEoP.Fill(
                best_d3_noEoP
            )


            # ------------------------------------------------
            # Position-based cluster coordinates
            # ------------------------------------------------

            cx_best, cy_best, cz_best = (
                best_cluster_noEoP
            )


            cluster_r = position_r(
                cx_best,
                cy_best,
                cz_best
            )


            cluster_theta = position_theta(
                cx_best,
                cy_best,
                cz_best
            )


            cluster_eta = position_eta(
                cx_best,
                cy_best,
                cz_best
            )


            # ------------------------------------------------
            # Fill NEW histograms
            # ------------------------------------------------

            h_dr_noEoP.Fill(
                cluster_r - proj_r
            )


            if (
                proj_theta is not None
                and
                cluster_theta is not None
            ):

                h_dtheta_noEoP.Fill(
                    cluster_theta -
                    proj_theta
                )


            if (
                proj_eta is not None
                and
                cluster_eta is not None
            ):

                h_deta_noEoP.Fill(
                    cluster_eta -
                    proj_eta
                )


    # ========================================================
    # Step 6b:
    # Nearest ECAL cluster WITH E/p cut
    # ========================================================

    best_d3_withEoP = None

    best_dxy_withEoP = None

    best_eop = None

    best_cluster_withEoP = None


    if (
        len(clusters) > 0
        and
        reco_p > 0
    ):

        best_d3_tmp = 1e30

        best_dxy_tmp = 1e30

        best_eop_tmp = None

        found_withEoP = False


        for cl in clusters:

            cx = float(
                cl.position.x
            )

            cy = float(
                cl.position.y
            )

            cz = float(
                cl.position.z
            )

            cE = float(
                cl.energy
            )


            dxy = math.sqrt(
                (cx - px) ** 2 +
                (cy - py) ** 2
            )


            d3 = math.sqrt(
                (cx - px) ** 2 +
                (cy - py) ** 2 +
                (cz - pz) ** 2
            )


            eop = cE / reco_p


            h_eop_all.Fill(
                eop
            )


            # ------------------------------------------------
            # E/p cut
            # ------------------------------------------------

            if not (
                EOP_MIN < eop < EOP_MAX
            ):
                continue


            if d3 < best_d3_tmp:

                best_d3_tmp = d3

                best_dxy_tmp = dxy

                best_eop_tmp = eop

                best_cluster_withEoP = (
                    cx,
                    cy,
                    cz
                )

                found_withEoP = True


        if found_withEoP:

            best_d3_withEoP = (
                best_d3_tmp
            )

            best_dxy_withEoP = (
                best_dxy_tmp
            )

            best_eop = (
                best_eop_tmp
            )


            n_scattered_with_cluster_withEoP += 1


            h_d3D_withEoP.Fill(
                best_d3_withEoP
            )


            h_eop_selected.Fill(
                best_eop
            )


            # ------------------------------------------------
            # Position-based cluster coordinates
            # ------------------------------------------------

            cx_best, cy_best, cz_best = (
                best_cluster_withEoP
            )


            cluster_r = position_r(
                cx_best,
                cy_best,
                cz_best
            )


            cluster_theta = position_theta(
                cx_best,
                cy_best,
                cz_best
            )


            cluster_eta = position_eta(
                cx_best,
                cy_best,
                cz_best
            )


            # ------------------------------------------------
            # Fill NEW histograms
            # ------------------------------------------------

            h_dr_withEoP.Fill(
                cluster_r -
                proj_r
            )


            if (
                proj_theta is not None
                and
                cluster_theta is not None
            ):

                h_dtheta_withEoP.Fill(
                    cluster_theta -
                    proj_theta
                )


            if (
                proj_eta is not None
                and
                cluster_eta is not None
            ):

                h_deta_withEoP.Fill(
                    cluster_eta -
                    proj_eta
                )


    # ========================================================
    # Record event results
    # ========================================================

    scattered_records.append({

        "event":
            ievt,

        "nearest_d3_noEoP":
            best_d3_noEoP,

        "nearest_d3_withEoP":
            best_d3_withEoP,

        "nearest_dxy_noEoP":
            best_dxy_noEoP,

        "nearest_dxy_withEoP":
            best_dxy_withEoP,

        "has_cluster_noEoP":
            best_d3_noEoP is not None,

        "has_cluster_withEoP":
            best_d3_withEoP is not None,

        "has_ecal_point":
            True,

        "reco_energy":
            reco_energy,

        "mc_energy":
            mc_energy,

        "parent_status":
            parent_status_selected,

        "best_eop":
            best_eop,
    })


    # ========================================================
    # Fraction vs distance cut
    # ========================================================

    for ibin in range(
        1,
        NBINS_CUT + 1
    ):

        dcut = h_total.GetBinCenter(
            ibin
        )


        h_total.Fill(
            dcut
        )


        if (
            best_d3_noEoP is not None
            and
            best_d3_noEoP < dcut
        ):

            h_pass_d3D_noEoP.Fill(
                dcut
            )


        if (
            best_d3_withEoP is not None
            and
            best_d3_withEoP < dcut
        ):

            h_pass_d3D_withEoP.Fill(
                dcut
            )


    # ========================================================
    # Failure at chosen cut WITH E/p
    # ========================================================

    passed_cut_withEoP = (
        best_d3_withEoP is not None
        and
        best_d3_withEoP < CHOSEN_CUT_MM
    )


    if not passed_cut_withEoP:

        if best_d3_withEoP is None:

            reason = "no_cluster_after_EoP"

        else:

            reason = "cluster_too_far_after_EoP"


        failed_records.append({

            "event":
                ievt,

            "reason":
                reason,

            "nearest_d3_noEoP":
                best_d3_noEoP,

            "nearest_d3_withEoP":
                best_d3_withEoP,

            "reco_energy":
                reco_energy,

            "mc_energy":
                mc_energy,

            "parent_status":
                parent_status_selected,

            "best_eop":
                best_eop,
        })


        h_fail_recoE.Fill(
            reco_energy
        )

        h_fail_mcE.Fill(
            mc_energy
        )


# ============================================================
# Build fraction histograms
# ============================================================

h_frac_d3D_noEoP = (
    h_pass_d3D_noEoP.Clone(
        "h_frac_d3D_noEoP"
    )
)

h_frac_d3D_noEoP.SetTitle(
    "Fraction of scattered electrons with nearby reconstructed "
    "ECal cluster (no E/p cut);"
    "d_{3D} cut [mm];Fraction"
)

h_frac_d3D_noEoP.Divide(
    h_total
)


h_frac_d3D_withEoP = (
    h_pass_d3D_withEoP.Clone(
        "h_frac_d3D_withEoP"
    )
)

h_frac_d3D_withEoP.SetTitle(
    f"Fraction of scattered electrons with nearby reconstructed "
    f"ECal cluster ({EOP_MIN:.1f} < E/p < {EOP_MAX:.1f});"
    f"d_{{3D}} cut [mm];Fraction"
)

h_frac_d3D_withEoP.Divide(
    h_total
)


# ============================================================
# Plateau information
# ============================================================

plateau_bin_noEoP = (
    h_frac_d3D_noEoP.GetMaximumBin()
)

plateau_fraction_noEoP = (
    h_frac_d3D_noEoP.GetBinContent(
        plateau_bin_noEoP
    )
)

plateau_cut_noEoP = (
    h_frac_d3D_noEoP.GetBinCenter(
        plateau_bin_noEoP
    )
)


first_plateau_bin_noEoP = None


for ibin in range(
    1,
    h_frac_d3D_noEoP.GetNbinsX() + 1
):

    if abs(
        h_frac_d3D_noEoP.GetBinContent(ibin)
        -
        plateau_fraction_noEoP
    ) < 1e-12:

        first_plateau_bin_noEoP = ibin

        break


first_plateau_cut_noEoP = (
    h_frac_d3D_noEoP.GetBinCenter(
        first_plateau_bin_noEoP
    )
    if first_plateau_bin_noEoP
    else None
)


# ------------------------------------------------------------

plateau_bin_withEoP = (
    h_frac_d3D_withEoP.GetMaximumBin()
)

plateau_fraction_withEoP = (
    h_frac_d3D_withEoP.GetBinContent(
        plateau_bin_withEoP
    )
)

plateau_cut_withEoP = (
    h_frac_d3D_withEoP.GetBinCenter(
        plateau_bin_withEoP
    )
)


first_plateau_bin_withEoP = None


for ibin in range(
    1,
    h_frac_d3D_withEoP.GetNbinsX() + 1
):

    if abs(
        h_frac_d3D_withEoP.GetBinContent(ibin)
        -
        plateau_fraction_withEoP
    ) < 1e-12:

        first_plateau_bin_withEoP = ibin

        break


first_plateau_cut_withEoP = (
    h_frac_d3D_withEoP.GetBinCenter(
        first_plateau_bin_withEoP
    )
    if first_plateau_bin_withEoP
    else None
)


# ============================================================
# Chosen d3D cut
# ============================================================

n_pass_cut_noEoP = 0
n_fail_cut_noEoP = 0

n_pass_cut_withEoP = 0
n_fail_cut_withEoP = 0


for rec in scattered_records:

    if not rec["has_ecal_point"]:

        n_fail_cut_noEoP += 1

        n_fail_cut_withEoP += 1

        continue


    d3_noEoP = (
        rec["nearest_d3_noEoP"]
    )

    d3_withEoP = (
        rec["nearest_d3_withEoP"]
    )


    if (
        d3_noEoP is not None
        and
        d3_noEoP < CHOSEN_CUT_MM
    ):

        n_pass_cut_noEoP += 1

    else:

        n_fail_cut_noEoP += 1


    if (
        d3_withEoP is not None
        and
        d3_withEoP < CHOSEN_CUT_MM
    ):

        n_pass_cut_withEoP += 1

    else:

        n_fail_cut_withEoP += 1


fraction_at_cut_noEoP = (
    n_pass_cut_noEoP /
    n_scattered_with_ecal_point

    if n_scattered_with_ecal_point > 0

    else 0.0
)


fraction_at_cut_withEoP = (
    n_pass_cut_withEoP /
    n_scattered_with_ecal_point

    if n_scattered_with_ecal_point > 0

    else 0.0
)


# ============================================================
# Print summary
# ============================================================

print(
    "\n========== Summary =========="
)

print(
    f"Total events                                 : "
    f"{n_total_events}"
)

print(
    f"All electrons with status=1 and PDG=11       : "
    f"{n_all_status1_electrons}"
)

print(
    f"Events with selected scattered electron      : "
    f"{n_events_with_scattered_electron}"
)

print(
    f"Scattered electrons with reco match          : "
    f"{n_scattered_with_reco}"
)

print(
    f"Scattered electrons with reco track          : "
    f"{n_scattered_with_track}"
)

print(
    f"Scattered electrons with projection          : "
    f"{n_scattered_with_projection}"
)

print(
    f"Scattered electrons with ECAL projection     : "
    f"{n_scattered_with_ecal_point}"
)

print(
    f"Scattered electrons with cluster (no E/p)    : "
    f"{n_scattered_with_cluster_noEoP}"
)

print(
    f"Scattered electrons with cluster (with E/p)  : "
    f"{n_scattered_with_cluster_withEoP}"
)


if h_d3D_noEoP.GetEntries() > 0:

    print(
        f"Mean nearest d3D (no E/p)                    : "
        f"{h_d3D_noEoP.GetMean():.3f} mm"
    )


if h_d3D_withEoP.GetEntries() > 0:

    print(
        f"Mean nearest d3D (with E/p)                  : "
        f"{h_d3D_withEoP.GetMean():.3f} mm"
    )


# ------------------------------------------------------------
# NEW coordinate-difference summary
# ------------------------------------------------------------

print(
    "\n========== Position-based differences =========="
)

if h_dr_noEoP.GetEntries() > 0:

    print(
        f"Mean delta r     (no E/p)                    : "
        f"{h_dr_noEoP.GetMean():.3f} mm"
    )

    print(
        f"Mean delta theta (no E/p)                    : "
        f"{h_dtheta_noEoP.GetMean():.6f} rad"
    )

    print(
        f"Mean delta eta   (no E/p)                    : "
        f"{h_deta_noEoP.GetMean():.6f}"
    )


if h_dr_withEoP.GetEntries() > 0:

    print(
        f"Mean delta r     (with E/p)                  : "
        f"{h_dr_withEoP.GetMean():.3f} mm"
    )

    print(
        f"Mean delta theta (with E/p)                  : "
        f"{h_dtheta_withEoP.GetMean():.6f} rad"
    )

    print(
        f"Mean delta eta   (with E/p)                  : "
        f"{h_deta_withEoP.GetMean():.6f}"
    )


# ============================================================
# Parent-status summary
# ============================================================

print(
    "\n========== Parent-status summary =========="
)


status_counts = {}


for ibin in range(
    1,
    h_parent_status_all.GetNbinsX() + 1
):

    count = int(
        h_parent_status_all.GetBinContent(
            ibin
        )
    )

    if count > 0:

        status_val = int(
            h_parent_status_all.GetBinCenter(
                ibin
            )
        )

        status_counts[
            status_val
        ] = count


for status_val in sorted(
    status_counts
):

    print(
        f"Parent generator status "
        f"{status_val:2d} : "
        f"{status_counts[status_val]}"
    )


bin2_all = (
    h_parent_status_all.FindBin(2)
)


n_parent2_before = int(
    h_parent_status_all.GetBinContent(
        bin2_all
    )
)


print(
    f"Electrons removed by parent_status != 2 cut : "
    f"{n_parent2_before}"
)


# ============================================================
# Plateau summary
# ============================================================

print(
    "\n========== Plateau info =========="
)


print(
    f"No E/p cut  - plateau fraction             : "
    f"{plateau_fraction_noEoP:.6f}"
)

print(
    f"No E/p cut  - first plateau cut            : "
    f"{first_plateau_cut_noEoP:.3f} mm"
)

print(
    f"With E/p    - plateau fraction             : "
    f"{plateau_fraction_withEoP:.6f}"
)

print(
    f"With E/p    - first plateau cut            : "
    f"{first_plateau_cut_withEoP:.3f} mm"
)


# ============================================================
# Chosen cut summary
# ============================================================

print(
    f"\n========== At chosen cut = "
    f"{CHOSEN_CUT_MM:.1f} mm =========="
)


print(
    f"No E/p cut  - passing scattered electrons  : "
    f"{n_pass_cut_noEoP}"
)

print(
    f"No E/p cut  - failing scattered electrons  : "
    f"{n_fail_cut_noEoP}"
)

print(
    f"No E/p cut  - fraction at chosen cut       : "
    f"{fraction_at_cut_noEoP:.6f}"
)


print(
    f"With E/p    - passing scattered electrons  : "
    f"{n_pass_cut_withEoP}"
)

print(
    f"With E/p    - failing scattered electrons  : "
    f"{n_fail_cut_withEoP}"
)

print(
    f"With E/p    - fraction at chosen cut       : "
    f"{fraction_at_cut_withEoP:.6f}"
)


# ============================================================
# Canvases
#
# ROOT only.
# No PNG files are created.
# ============================================================

# ------------------------------------------------------------
# d3D without E/p
# ------------------------------------------------------------

c1 = ROOT.TCanvas(
    "c1",
    "d3D without E/p",
    900,
    700
)

h_d3D_noEoP.SetLineWidth(2)

h_d3D_noEoP.Draw(
    "HIST"
)


# ------------------------------------------------------------
# d3D with E/p
# ------------------------------------------------------------

c2 = ROOT.TCanvas(
    "c2",
    "d3D with E/p",
    900,
    700
)

h_d3D_withEoP.SetLineWidth(2)

h_d3D_withEoP.SetLineColor(
    ROOT.kBlue + 1
)

h_d3D_withEoP.Draw(
    "HIST"
)


# ------------------------------------------------------------
# Fraction vs d3D cut
# ------------------------------------------------------------

c3 = ROOT.TCanvas(
    "c3",
    "fraction vs distance cut",
    900,
    700
)


h_frac_d3D_noEoP.SetLineWidth(
    2
)

h_frac_d3D_noEoP.SetLineColor(
    ROOT.kBlack
)

h_frac_d3D_noEoP.SetMinimum(
    0.0
)

h_frac_d3D_noEoP.SetMaximum(
    1.05
)

h_frac_d3D_noEoP.Draw(
    "HIST"
)


h_frac_d3D_withEoP.SetLineWidth(
    2
)

h_frac_d3D_withEoP.SetLineColor(
    ROOT.kRed
)

h_frac_d3D_withEoP.Draw(
    "HIST SAME"
)


leg = ROOT.TLegend(
    0.52,
    0.73,
    0.88,
    0.88
)

leg.SetBorderSize(0)

leg.SetFillStyle(0)

leg.AddEntry(
    h_frac_d3D_noEoP,
    "No E/p cut",
    "l"
)

leg.AddEntry(
    h_frac_d3D_withEoP,
    f"{EOP_MIN:.1f} < E/p < {EOP_MAX:.1f}",
    "l"
)

leg.Draw()


# ------------------------------------------------------------
# Failed electron energies
# ------------------------------------------------------------

c4 = ROOT.TCanvas(
    "c4",
    "failing scattered electron energies",
    1200,
    500
)

c4.Divide(
    2,
    1
)

c4.cd(1)

h_fail_recoE.SetLineWidth(
    2
)

h_fail_recoE.Draw(
    "HIST"
)


c4.cd(2)

h_fail_mcE.SetLineWidth(
    2
)

h_fail_mcE.Draw(
    "HIST"
)


# ------------------------------------------------------------
# Parent status
# ------------------------------------------------------------

c5 = ROOT.TCanvas(
    "c5",
    "parent status histograms",
    1200,
    500
)

c5.Divide(
    2,
    1
)


c5.cd(1)

h_parent_status_all.SetLineWidth(
    2
)

h_parent_status_all.Draw(
    "HIST"
)


c5.cd(2)

h_parent_status_selected.SetLineWidth(
    2
)

h_parent_status_selected.SetLineColor(
    ROOT.kBlue + 1
)

h_parent_status_selected.Draw(
    "HIST"
)


# ------------------------------------------------------------
# Parent PDG
# ------------------------------------------------------------

c6 = ROOT.TCanvas(
    "c6",
    "parent PDG all electrons",
    900,
    700
)

h_parent_pdg_all.SetLineWidth(
    2
)

h_parent_pdg_all.Draw(
    "HIST"
)


# ------------------------------------------------------------
# E/p
# ------------------------------------------------------------

c7 = ROOT.TCanvas(
    "c7",
    "E over p",
    1200,
    500
)

c7.Divide(
    2,
    1
)


c7.cd(1)

h_eop_all.SetLineWidth(
    2
)

h_eop_all.Draw(
    "HIST"
)


c7.cd(2)

h_eop_selected.SetLineWidth(
    2
)

h_eop_selected.SetLineColor(
    ROOT.kBlue + 1
)

h_eop_selected.Draw(
    "HIST"
)


# ============================================================
# NEW CANVAS:
# r, theta, eta from positions
# ============================================================

c8 = ROOT.TCanvas(
    "c8",
    "position based cluster projection differences",
    1500,
    900
)

c8.Divide(
    3,
    2
)


# ------------------------------------------------------------
# Top row: NO E/p
# ------------------------------------------------------------

c8.cd(1)

h_dr_noEoP.SetLineWidth(
    2
)

h_dr_noEoP.Draw(
    "HIST"
)


c8.cd(2)

h_dtheta_noEoP.SetLineWidth(
    2
)

h_dtheta_noEoP.Draw(
    "HIST"
)


c8.cd(3)

h_deta_noEoP.SetLineWidth(
    2
)

h_deta_noEoP.Draw(
    "HIST"
)


# ------------------------------------------------------------
# Bottom row: WITH E/p
# ------------------------------------------------------------

c8.cd(4)

h_dr_withEoP.SetLineWidth(
    2
)

h_dr_withEoP.SetLineColor(
    ROOT.kBlue + 1
)

h_dr_withEoP.Draw(
    "HIST"
)


c8.cd(5)

h_dtheta_withEoP.SetLineWidth(
    2
)

h_dtheta_withEoP.SetLineColor(
    ROOT.kBlue + 1
)

h_dtheta_withEoP.Draw(
    "HIST"
)


c8.cd(6)

h_deta_withEoP.SetLineWidth(
    2
)

h_deta_withEoP.SetLineColor(
    ROOT.kBlue + 1
)

h_deta_withEoP.Draw(
    "HIST"
)


# ============================================================
# Save ROOT file
# ============================================================

out = ROOT.TFile(
    OUTROOT,
    "RECREATE"
)


# d3D
h_d3D_noEoP.Write()

h_d3D_withEoP.Write()


# E/p
h_eop_all.Write()

h_eop_selected.Write()


# fractions
h_total.Write()

h_pass_d3D_noEoP.Write()

h_pass_d3D_withEoP.Write()

h_frac_d3D_noEoP.Write()

h_frac_d3D_withEoP.Write()


# energies
h_fail_recoE.Write()

h_fail_mcE.Write()


# parent information
h_parent_status_all.Write()

h_parent_status_selected.Write()

h_parent_pdg_all.Write()


# ------------------------------------------------------------
# NEW position-based histograms
# ------------------------------------------------------------

h_dr_noEoP.Write()

h_dtheta_noEoP.Write()

h_deta_noEoP.Write()


h_dr_withEoP.Write()

h_dtheta_withEoP.Write()

h_deta_withEoP.Write()


# canvases
c1.Write()

c2.Write()

c3.Write()

c4.Write()

c5.Write()

c6.Write()

c7.Write()

c8.Write()


out.Close()


print(
    f"\n✅ Wrote ROOT file: "
    f"{OUTROOT}"
)
