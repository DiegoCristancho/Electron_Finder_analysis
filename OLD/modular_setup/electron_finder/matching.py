"""Associations and nearest-cluster matching, independent of ROOT I/O."""

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class ClusterMatch:
    index: int
    position: tuple
    distance_mm: float
    transverse_distance_mm: float
    energy: float
    eop: object


def find_projection_index(projections, track_refs, target_track):
    """Use the first projection for the same (collectionID, index) track."""
    for index in range(len(projections)):
        if same_track_ref(track_refs[index], target_track):
            return index
    return None


def find_ecal_point(projection, points, system):
    """Return (point index, XYZ) for the first point in the requested system."""
    for index in range(int(projection.points_begin), int(projection.points_end)):
        point = points[index]
        if int(point.system) == system:
            return index, (float(point.position.x), float(point.position.y), float(point.position.z))
    return None, None


def match_clusters(clusters, position, momentum, eop_min, eop_max):
    """Find nearest in 3D both before and after the strict E/p window.

    Filtering by E/p precedes choosing the nearest qualifying cluster. Equal
    distances retain the first cluster, as in the original analysis. The third
    return value contains E/p for every cluster, for the inclusive histogram.
    """
    best = best_eop = None
    best_distance = best_eop_distance = 1e30
    eop_values = []
    px, py, pz = position
    for index, cluster in enumerate(clusters):
        cx, cy, cz = float(cluster.position.x), float(cluster.position.y), float(cluster.position.z)
        dxy = math.sqrt((cx - px) ** 2 + (cy - py) ** 2)
        d3 = math.sqrt((cx - px) ** 2 + (cy - py) ** 2 + (cz - pz) ** 2)
        energy = float(cluster.energy)
        eop = energy / momentum if momentum > 0 else None
        candidate = ClusterMatch(index, (cx, cy, cz), d3, dxy, energy, eop)
        if d3 < best_distance:
            best, best_distance = candidate, d3
        if eop is not None:
            eop_values.append(eop)
            if eop_min < eop < eop_max and d3 < best_eop_distance:
                best_eop, best_eop_distance = candidate, d3
    return best, best_eop, eop_values


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
    return int(a.collectionID) == int(b.collectionID) and int(a.index) == int(b.index)


def find_reco_particle_for_track(track_ref, reco_particles, reco_particle_tracks):
    """
    Given one track reference from _CalorimeterTrackProjections_track,
    find which ReconstructedParticle owns that track.
    """
    trk_coll = int(track_ref.collectionID)
    trk_idx = int(track_ref.index)
    for (ireco, reco) in enumerate(reco_particles):
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
