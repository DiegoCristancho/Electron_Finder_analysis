"""Readable terminal reports, using the same sections as the original Trial.py."""


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
