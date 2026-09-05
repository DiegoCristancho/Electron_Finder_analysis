"""Analysis, physics calculations, selections, cuts, and printed reports.

Edit AnalysisConfig below for the default cuts and collections.
Edit ElectronAnalysis.process_event for the full-file analysis sequence.
Edit find_first_scattered_electron_mc_index or match_clusters for selection rules.
Edit inspect_event for the single-event data to inspect.
Histogram definitions, figures and ROOT file writing are in histograms.py.
"""

from collections import Counter
from dataclasses import asdict, dataclass, fields
import json
import math
from pathlib import Path

from analysis_io import EventWriter, event_source, iter_events
from histograms import book_histograms, save_analysis_results



# ========================================================================
# SETTINGS: edit the defaults here; CLI options can override them
# ========================================================================

@dataclass(frozen=True)
class AnalysisConfig:
    event_tree: str = "events"
    cluster_collection: str = "EcalBarrelClusters"
    ecal_system: int = 101
    distance_cut_mm: float = 87.0
    eop_min: float = 0.7
    eop_max: float = 1.3
    cut_scan_bins: int = 100
    cut_scan_max_mm: float = 200.0

    def __post_init__(self):
        for name in ("event_tree", "cluster_collection"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{name} must be a nonempty string")
        for name in ("ecal_system", "cut_scan_bins"):
            if type(getattr(self, name)) is not int:
                raise ValueError(f"{name} must be an integer")
        if self.ecal_system < 0 or self.cut_scan_bins < 1:
            raise ValueError("ecal_system must be >= 0 and cut_scan_bins must be >= 1")
        for name in ("distance_cut_mm", "eop_min", "eop_max", "cut_scan_max_mm"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f"{name} must be a finite number")
        if self.distance_cut_mm <= 0 or self.cut_scan_max_mm <= 0:
            raise ValueError("Distance cut and scan maximum must be positive")
        if not 0 <= self.eop_min < self.eop_max:
            raise ValueError("Require 0 <= eop_min < eop_max")

    @classmethod
    def load(cls, filename=None, **overrides):
        values = {}
        if filename is not None:
            values = json.loads(Path(filename).read_text())
            if not isinstance(values, dict):
                raise ValueError("Configuration must be a JSON object")
        unknown = set(values) - {field.name for field in fields(cls)}
        if unknown:
            raise ValueError(f"Unknown configuration settings: {', '.join(sorted(unknown))}")
        values.update({key: value for key, value in overrides.items() if value is not None})
        return cls(**values)

    def to_dict(self):
        return asdict(self)


# ========================================================================
# FULL-FILE ANALYSIS: the per-event sequence, counts and fractions
# ========================================================================

def required_branches(config):
    return (
        "MCParticles", "_MCParticles_parents", "CalorimeterTrackProjections",
        "_CalorimeterTrackProjections_track", "_CalorimeterTrackProjections_points",
        config.cluster_collection, "ReconstructedParticles", "_ReconstructedParticles_tracks",
        "_ReconstructedParticleAssociations_rec", "_ReconstructedParticleAssociations_sim",
    )


class ElectronAnalysis:
    """Create one instance per run. process_event returns only copied Python values."""

    def __init__(self, config=None):
        self.config = config or AnalysisConfig()
        self.histograms = book_histograms(self.config)
        self._finalized = False
        self.counts = Counter({name: 0 for name in (
            "events_processed", "status1_electrons", "selected_electrons", "with_reco",
            "with_track", "with_projection", "with_ecal_point", "with_cluster_no_eop",
            "with_cluster_with_eop", "pass_distance_no_eop", "pass_distance_with_eop",
        )})
        self.failures = Counter()

    def _fail(self, record, reason):
        record["reason"] = reason
        self.failures[reason] += 1
        if record["mc_energy"] is not None:
            self.histograms["h_fail_mcE"].Fill(record["mc_energy"])
        if record["reco_energy"] is not None:
            self.histograms["h_fail_recoE"].Fill(record["reco_energy"])
        return record

    def _fill_match(self, match, projection, suffix):
        hist = self.histograms
        hist[f"h_d3D_{suffix}"].Fill(match.distance_mm)
        hist[f"h_dr_{suffix}"].Fill(position_r(*match.position) - position_r(*projection))
        for name, function in (("dtheta", position_theta), ("deta", position_eta)):
            cluster_value, projection_value = function(*match.position), function(*projection)
            if cluster_value is not None and projection_value is not None:
                hist[f"h_{name}_{suffix}"].Fill(cluster_value - projection_value)

    def process_event(self, event, event_index=0):
        """Apply the original selections, fill histograms, and return an event record."""
        if self._finalized:
            raise RuntimeError("This run is finalized. Create a new ElectronAnalysis for another run.")
        missing = [name for name in required_branches(self.config) if not hasattr(event, name)]
        if missing:
            raise RuntimeError(f"Event {event_index}: missing branches: {', '.join(missing)}")
        self.counts["events_processed"] += 1
        record = {
            "event": event_index, "reason": None, "selected": False,
            "mc_index": None, "reco_index": None, "track_index": None,
            "track_collection_id": None, "projection_index": None, "ecal_point_index": None,
            "projection_xyz_mm": None, "has_ecal_point": False,
            "mc_energy": None, "reco_energy": None, "reco_momentum": None,
            "parent_status": None, "parent_pdg": None,
            "cluster_index_noEoP": None, "cluster_index_withEoP": None,
            "nearest_d3_noEoP": None, "nearest_d3_withEoP": None,
            "nearest_dxy_noEoP": None, "nearest_dxy_withEoP": None,
            "has_cluster_noEoP": False, "has_cluster_withEoP": False,
            "best_eop": None, "passed_noEoP": False, "passed_withEoP": False,
        }
        mc = event.MCParticles
        parents = event._MCParticles_parents
        hist = self.histograms
        for particle in mc:
            if int(particle.generatorStatus) == 1 and int(particle.PDG) == 11:
                self.counts["status1_electrons"] += 1
                status = get_first_parent_status(particle, mc, parents)
                pdg = get_first_parent_pdg(particle, mc, parents)
                if status is not None:
                    hist["h_parent_status_all"].Fill(status)
                if pdg is not None:
                    hist["h_parent_pdg_all"].Fill(pdg)

        mc_index = find_first_scattered_electron_mc_index(mc, parents)
        if mc_index is None:
            return self._fail(record, "no_scattered_electron")
        particle = mc[mc_index]
        self.counts["selected_electrons"] += 1
        record.update(selected=True, mc_index=mc_index, mc_energy=particle_energy(particle),
                      parent_status=get_first_parent_status(particle, mc, parents),
                      parent_pdg=get_first_parent_pdg(particle, mc, parents))
        if record["parent_status"] is not None:
            hist["h_parent_status_selected"].Fill(record["parent_status"])

        reco_index = find_reco_index_for_mc_index(
            mc_index, event._ReconstructedParticleAssociations_rec,
            event._ReconstructedParticleAssociations_sim,
        )
        record["reco_index"] = reco_index
        reco_particles = event.ReconstructedParticles
        if reco_index is None or not 0 <= reco_index < len(reco_particles):
            return self._fail(record, "no_reco_match")
        self.counts["with_reco"] += 1
        reco = reco_particles[reco_index]
        momentum = particle_p(reco)
        record.update(reco_energy=particle_energy(reco), reco_momentum=momentum)

        track = find_track_ref_for_reco_index(reco_index, reco_particles, event._ReconstructedParticles_tracks)
        if track is None:
            return self._fail(record, "no_track_on_reco_particle")
        self.counts["with_track"] += 1
        record.update(track_index=int(track.index), track_collection_id=int(track.collectionID))

        projections = event.CalorimeterTrackProjections
        projection_index = find_projection_index(projections, event._CalorimeterTrackProjections_track, track)
        if projection_index is None:
            return self._fail(record, "no_projection_for_scattered_electron_track")
        self.counts["with_projection"] += 1
        record["projection_index"] = projection_index
        point_index, position = find_ecal_point(
            projections[projection_index], event._CalorimeterTrackProjections_points, self.config.ecal_system,
        )
        if position is None:
            return self._fail(record, "no_ecal_projection_point")
        self.counts["with_ecal_point"] += 1
        record.update(has_ecal_point=True, ecal_point_index=point_index, projection_xyz_mm=position)

        nearest, nearest_eop, eop_values = match_clusters(
            getattr(event, self.config.cluster_collection), position, momentum,
            self.config.eop_min, self.config.eop_max,
        )
        for eop in eop_values:
            hist["h_eop_all"].Fill(eop)
        for suffix, count_suffix, match in (("noEoP", "no_eop", nearest), ("withEoP", "with_eop", nearest_eop)):
            if match is None:
                continue
            self.counts[f"with_cluster_{count_suffix}"] += 1
            self._fill_match(match, position, suffix)
            record[f"cluster_index_{suffix}"] = match.index
            record[f"nearest_d3_{suffix}"] = match.distance_mm
            record[f"nearest_dxy_{suffix}"] = match.transverse_distance_mm
            record[f"has_cluster_{suffix}"] = True
            record[f"passed_{suffix}"] = match.distance_mm < self.config.distance_cut_mm
            if record[f"passed_{suffix}"]:
                self.counts[f"pass_distance_{count_suffix}"] += 1
        if nearest_eop is not None:
            record["best_eop"] = nearest_eop.eop
            hist["h_eop_selected"].Fill(nearest_eop.eop)

        # Only electrons that reached an ECAL projection enter this denominator.
        for bin_index in range(1, self.config.cut_scan_bins + 1):
            cut = hist["h_total"].GetBinCenter(bin_index)
            hist["h_total"].Fill(cut)
            for suffix, match in (("noEoP", nearest), ("withEoP", nearest_eop)):
                if match is not None and match.distance_mm < cut:
                    hist[f"h_pass_d3D_{suffix}"].Fill(cut)
        if not record["passed_withEoP"]:
            reason = "no_cluster_after_EoP" if nearest_eop is None else "cluster_too_far_after_EoP"
            return self._fail(record, reason)
        return record

    def run(self, chain, output_dir, max_events=None, progress_every=1000, *, event_index=None):
        """Stream all records to disk; return a compact summary of this run."""
        with EventWriter(output_dir) as writer:
            events = iter_events(chain, max_events, event_index=event_index)
            for processed, (index, event) in enumerate(events, start=1):
                record = self.process_event(event, index)
                record.update(event_source(chain, index))
                writer.write(record)
                if progress_every and processed % progress_every == 0:
                    print(f"Processed {processed} events", flush=True)
        self.finalize()
        return self.summary()

    def finalize(self):
        if self._finalized:
            return
        for suffix in ("noEoP", "withEoP"):
            name = f"h_frac_d3D_{suffix}"
            if name in self.histograms:
                continue
            histogram = self.histograms[f"h_pass_d3D_{suffix}"].Clone(name)
            histogram.SetDirectory(0)
            label = "no E/p cut" if suffix == "noEoP" else f"{self.config.eop_min:g} < E/p < {self.config.eop_max:g}"
            histogram.SetTitle(
                f"Fraction of scattered electrons with nearby reconstructed ECal cluster ({label});"
                "d_{3D} cut [mm];Fraction"
            )
            histogram.Divide(self.histograms["h_total"])
            self.histograms[name] = histogram
        self._finalized = True

    def summary(self):
        self.finalize()
        mean_names = [f"h_{quantity}_{suffix}"
                      for quantity in ("d3D", "dr", "dtheta", "deta")
                      for suffix in ("noEoP", "withEoP")]
        histogram_means = {
            name: float(self.histograms[name].GetMean()) if self.histograms[name].GetEntries() else None
            for name in mean_names
        }
        parent_histogram = self.histograms["h_parent_status_all"]
        parent_status_counts = {}
        for index in range(1, parent_histogram.GetNbinsX() + 1):
            count = int(parent_histogram.GetBinContent(index))
            if count:
                parent_status_counts[str(int(parent_histogram.GetBinCenter(index)))] = count
        projected = self.counts["with_ecal_point"]
        at_cut = {}
        for suffix in ("no_eop", "with_eop"):
            passed = self.counts[f"pass_distance_{suffix}"]
            at_cut[suffix] = {
                "passed": passed, "failed_among_ecal_projected": projected - passed,
                "denominator_ecal_projected": projected,
                "fraction": passed / projected if projected else 0.0,
                # The original printed fail count also included missing ECAL points.
                "legacy_failed_count": projected - passed + self.failures["no_ecal_projection_point"],
            }
        plateau = {}
        for suffix in ("noEoP", "withEoP"):
            histogram = self.histograms[f"h_frac_d3D_{suffix}"]
            maximum = histogram.GetBinContent(histogram.GetMaximumBin())
            first = next((i for i in range(1, histogram.GetNbinsX() + 1)
                          if abs(histogram.GetBinContent(i) - maximum) < 1e-12), None)
            plateau[suffix] = {"fraction": maximum, "first_cut_mm": histogram.GetBinCenter(first) if first else None}
        return {
            "cutflow": dict(self.counts), "failure_reasons": dict(self.failures),
            "distance_cut_mm": self.config.distance_cut_mm,
            "at_distance_cut": at_cut, "plateau": plateau,
            "histogram_means": histogram_means,
            "parent_status_counts": parent_status_counts,
            "parent_status_2_electrons": int(parent_histogram.GetBinContent(parent_histogram.FindBin(2))),
            "selection": "First MC electron with PDG=11, generatorStatus=1, first-parent status != 2; missing parent passes",
            "fraction_denominator": "Selected electrons with a reconstructed ECAL projection point",
        }

    def print_summary(self):
        """Print the detailed Trial.py-style summary for this completed run."""
        print(format_analysis_summary(self.summary()), end="", flush=True)

    def save_results(self, directory, make_plots=True):
        """Send the finished histograms and report to histograms.py for output."""
        summary = self.summary()
        save_analysis_results(self.histograms, self.config, summary,
                              format_analysis_summary(summary), directory, make_plots)


# ========================================================================
# PHYSICS: momentum/energy in GeV; positions in mm
# ========================================================================

def particle_energy(obj):
    px = float(obj.momentum.x)
    py = float(obj.momentum.y)
    pz = float(obj.momentum.z)
    m = float(obj.mass) if hasattr(obj, 'mass') else 0.0
    return math.sqrt(px * px + py * py + pz * pz + m * m)


def particle_p(obj):
    px = float(obj.momentum.x)
    py = float(obj.momentum.y)
    pz = float(obj.momentum.z)
    return math.sqrt(px * px + py * py + pz * pz)


def position_r(x, y, z):
    """
    3D spherical radius:
        r = sqrt(x^2 + y^2 + z^2)
    """
    return math.sqrt(x * x + y * y + z * z)


def position_theta(x, y, z):
    """
    Polar angle calculated from position:
        theta = acos(z/r)
    """
    r = position_r(x, y, z)
    if r <= 0:
        return None
    cos_theta = z / r
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


# ========================================================================
# ELECTRON SELECTION: first status-1 electron, first-parent status != 2
# ========================================================================

# A missing/invalid parent index passes the original veto.

def get_first_parent_status(part, mc_particles, parent_refs):
    pb = int(part.parents_begin)
    pe = int(part.parents_end)
    if pe <= pb:
        return None
    pref = parent_refs[pb]
    parent_idx = int(pref.index)
    if parent_idx < 0 or parent_idx >= len(mc_particles):
        return None
    return int(mc_particles[parent_idx].generatorStatus)


def get_first_parent_pdg(part, mc_particles, parent_refs):
    pb = int(part.parents_begin)
    pe = int(part.parents_end)
    if pe <= pb:
        return None
    pref = parent_refs[pb]
    parent_idx = int(pref.index)
    if parent_idx < 0 or parent_idx >= len(mc_particles):
        return None
    return int(mc_particles[parent_idx].PDG)


def find_first_scattered_electron_mc_index(mc_particles, parent_refs):
    """
    Scattered-electron definition:

        generatorStatus == 1
        PDG == 11
        parent_status != 2

    Take the first electron satisfying these conditions.
    """
    for (i, part) in enumerate(mc_particles):
        if int(part.generatorStatus) != 1:
            continue
        if int(part.PDG) != 11:
            continue
        parent_status = get_first_parent_status(part, mc_particles, parent_refs)
        if parent_status != 2:
            return i
    return None


# ========================================================================
# TRACK / CLUSTER MATCHING: associations, projections and E/p cuts
# ========================================================================

@dataclass(frozen=True)
class ClusterMatch:
    index: int
    position: tuple
    distance_mm: float
    transverse_distance_mm: float
    energy: float
    eop: object


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


# ========================================================================
# SINGLE-EVENT INSPECTION: all projected tracks and ECal clusters
# ========================================================================

PROJECTION_BRANCHES = (
    "CalorimeterTrackProjections", "_CalorimeterTrackProjections_track",
    "_CalorimeterTrackProjections_points",
)

def display_branches(config, tracks_only=False):
    if tracks_only:
        return PROJECTION_BRANCHES
    return PROJECTION_BRANCHES + (
        config.cluster_collection, "ReconstructedParticles", "_ReconstructedParticles_tracks",
        "_ReconstructedParticleAssociations_rec", "_ReconstructedParticleAssociations_sim", "MCParticles",
    )


def inspect_event(event, index, config, tracks_only=False):
    """Copy positions and associations out of ROOT's current event view."""
    tracks = []
    points = event._CalorimeterTrackProjections_points
    for projection_index, projection in enumerate(event.CalorimeterTrackProjections):
        track = event._CalorimeterTrackProjections_track[projection_index]
        reco_index = truth_pdg = None
        if not tracks_only:
            reco_index = find_reco_particle_for_track(
                track, event.ReconstructedParticles, event._ReconstructedParticles_tracks,
            )
            truth_pdg = find_mc_pdg_for_reco_index(
                reco_index, event._ReconstructedParticleAssociations_rec,
                event._ReconstructedParticleAssociations_sim, event.MCParticles,
            )
        positions = []
        for point_index in range(int(projection.points_begin), int(projection.points_end)):
            point = points[point_index]
            x, y, z = float(point.position.x), float(point.position.y), float(point.position.z)
            positions.append({
                "point_index": point_index, "x": x, "y": y, "z": z,
                "rho": math.sqrt(x*x + y*y), "system": int(point.system), "surface": int(point.surface),
            })
        tracks.append({
            "projection_index": projection_index, "track_index": int(track.index),
            "collection_id": int(track.collectionID), "reco_index": reco_index,
            "truth_pdg": truth_pdg, "points_begin": int(projection.points_begin),
            "points_end": int(projection.points_end), "points": positions,
        })
    clusters = []
    if not tracks_only:
        for cluster_index, cluster in enumerate(getattr(event, config.cluster_collection)):
            x, y, z = float(cluster.position.x), float(cluster.position.y), float(cluster.position.z)
            clusters.append({"index": cluster_index, "x": x, "y": y, "z": z,
                             "rho": math.sqrt(x*x + y*y), "energy": float(cluster.energy)})
    counts = {
        "projections": len(event.CalorimeterTrackProjections),
        "track_refs": len(event._CalorimeterTrackProjections_track),
        "projection_points": len(points),
    }
    if not tracks_only:
        counts.update(
            clusters=len(clusters), reconstructed_particles=len(event.ReconstructedParticles),
            reconstructed_particle_tracks=len(event._ReconstructedParticles_tracks),
            reco_mc_associations=len(event._ReconstructedParticleAssociations_rec),
            mc_particles=len(event.MCParticles),
        )
    return {"event": index, "tracks_only": tracks_only, "tracks": tracks,
            "cluster_collection": config.cluster_collection, "clusters": clusters, "counts": counts}


# ========================================================================
# PRINTED REPORTS: Trial.py summary and Analysis.py event details
# ========================================================================

def format_analysis_summary(summary):
    """Render a completed summary without ROOT or access to the input files."""
    lines = []

    def heading(title):
        lines.extend(["", f"========== {title} =========="])

    def row(label, value):
        lines.append(f"{label:<46}: {value}")

    def number(value, digits=3, unit=""):
        return "N/A (no entries)" if value is None else f"{value:.{digits}f}{unit}"

    counts = summary["cutflow"]
    means = summary["histogram_means"]
    heading("Summary")
    for label, key in (
        ("Total events", "events_processed"),
        ("All electrons with status=1 and PDG=11", "status1_electrons"),
        ("Events with selected scattered electron", "selected_electrons"),
        ("Scattered electrons with reco match", "with_reco"),
        ("Scattered electrons with reco track", "with_track"),
        ("Scattered electrons with projection", "with_projection"),
        ("Scattered electrons with ECAL projection", "with_ecal_point"),
        ("Scattered electrons with cluster (no E/p)", "with_cluster_no_eop"),
        ("Scattered electrons with cluster (with E/p)", "with_cluster_with_eop"),
    ):
        row(label, counts[key])
    row("Mean nearest d3D (no E/p)", number(means["h_d3D_noEoP"], unit=" mm"))
    row("Mean nearest d3D (with E/p)", number(means["h_d3D_withEoP"], unit=" mm"))

    heading("Position-based differences")
    for label, suffix in (("no E/p", "noEoP"), ("with E/p", "withEoP")):
        row(f"Mean delta r ({label})", number(means[f"h_dr_{suffix}"], unit=" mm"))
        row(f"Mean delta theta ({label})", number(means[f"h_dtheta_{suffix}"], 6, " rad"))
        row(f"Mean delta eta ({label})", number(means[f"h_deta_{suffix}"], 6))

    heading("Parent-status summary")
    statuses = summary["parent_status_counts"]
    for status in sorted(statuses, key=int):
        row(f"Parent generator status {int(status):2d}", statuses[status])
    if not statuses:
        lines.append("No parent statuses recorded.")
    row("Electrons removed by parent_status != 2 cut", summary["parent_status_2_electrons"])

    heading("Plateau info")
    for label, suffix in (("No E/p cut", "noEoP"), ("With E/p", "withEoP")):
        plateau = summary["plateau"][suffix]
        row(f"{label:<11}- plateau fraction", number(plateau["fraction"], 6))
        row(f"{label:<11}- first plateau cut", number(plateau["first_cut_mm"], unit=" mm"))

    heading(f"At chosen cut = {summary['distance_cut_mm']:g} mm")
    for label, suffix in (("No E/p cut", "no_eop"), ("With E/p", "with_eop")):
        result = summary["at_distance_cut"][suffix]
        row(f"{label:<11}- passing scattered electrons", result["passed"])
        row(f"{label:<11}- failing scattered electrons", result["legacy_failed_count"])
        row(f"{label:<11}- fraction at chosen cut", number(result["fraction"], 6))
    lines.append(f"Fraction denominator: {counts['with_ecal_point']} scattered electrons with an ECAL projection point.")
    lines.append("Failing counts follow Trial.py: they also include missing ECAL points.")

    heading("Failure reasons")
    for reason, count in sorted(summary["failure_reasons"].items()):
        row(reason, count)
    if not summary["failure_reasons"]:
        lines.append("None.")
    return "\n".join(lines) + "\n"


def format_event_details(record):
    """The per-projection/point/cluster inspection originally printed by Analysis.py."""
    lines = [f"\n========== Event {record['event']} =========="]
    if record.get("source_file") is not None:
        lines.append(f"Source file: {record['source_file']}")
        lines.append(f"Entry in source file: {record['source_entry']}")
    for label, key in (
        ("n projections", "projections"), ("n track refs", "track_refs"),
        ("n total projection points", "projection_points"), ("n ECal clusters", "clusters"),
        ("n ReconstructedParticles", "reconstructed_particles"),
        ("n _ReconstructedParticles_tracks", "reconstructed_particle_tracks"),
        ("n Reco->MC associations", "reco_mc_associations"), ("n MCParticles", "mc_particles"),
    ):
        if key in record["counts"]:
            lines.append(f"{label:<32} = {record['counts'][key]}")
    for track in record["tracks"]:
        lines.extend([
            "", f"Projection {track['projection_index']}",
            f"  track.index        = {track['track_index']}",
            f"  track.collectionID = {track['collection_id']}",
            f"  points_begin       = {track['points_begin']}",
            f"  points_end         = {track['points_end']}",
            f"  n points           = {len(track['points'])}",
        ])
        if not record["tracks_only"]:
            lines.extend([f"  matched reco index = {track['reco_index']}",
                          f"  truth PDG          = {track['truth_pdg']}"])
        for point in track["points"]:
            lines.append(
                f"    point {point['point_index']}: system={point['system']}, surface={point['surface']}, "
                f"x={point['x']:.3f}, y={point['y']:.3f}, z={point['z']:.3f}, r={point['rho']:.3f}"
            )
    if not record["tracks_only"]:
        lines.extend(["", record["cluster_collection"]])
        for cluster in record["clusters"]:
            lines.append(
                f"  cluster {cluster['index']}: x={cluster['x']:.3f}, y={cluster['y']:.3f}, "
                f"z={cluster['z']:.3f}, r={cluster['rho']:.3f}, E={cluster['energy']:.3f} GeV"
            )
    return "\n".join(lines) + "\n"
