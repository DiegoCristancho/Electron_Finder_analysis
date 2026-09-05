"""MC selection conventions preserved from the original matching study.

Only the first parent is examined. A missing/invalid parent returns None, which
also passes the original parent_status != 2 selection. This is a truth-seeded
matching study; it is not yet a reconstruction-only electron identification.
"""


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
