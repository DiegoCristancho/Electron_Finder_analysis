"""State for one matching study; physics and matching helpers remain reusable."""

from collections import Counter
from pathlib import Path

from .config import AnalysisConfig
from .histograms import book_histograms
from .io import EventWriter, event_source, iter_events, root_output, write_json
from .matching import (
    find_ecal_point, find_projection_index, find_reco_index_for_mc_index,
    find_track_ref_for_reco_index, match_clusters,
)
from .physics import particle_energy, particle_p, position_eta, position_r, position_theta
from .reporting import format_analysis_summary
from .truth import (
    find_first_scattered_electron_mc_index, get_first_parent_pdg, get_first_parent_status,
)


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
        self.finalize()
        directory = Path(directory)
        canvases = {}
        if make_plots:
            from .plotting import plot_analysis
            canvases = plot_analysis(self.histograms, self.config, directory / "plots")
        with root_output(directory / "analysis.root"):
            for histogram in self.histograms.values():
                histogram.Write()
            for canvas in canvases.values():
                canvas.Write()
        summary = self.summary()
        write_json(directory / "summary.json", summary)
        (directory / "summary.txt").write_text(format_analysis_summary(summary))
