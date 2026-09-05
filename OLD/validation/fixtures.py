"""Small event records covering each branch of the original selection."""

from types import SimpleNamespace as NS


class EventChain:
    """A readable chain of fixture events that records which entries were loaded."""

    def __init__(self, events):
        self.events = events
        self.read_entries = []
        self.index = None

    def GetEntries(self):
        return len(self.events)

    def GetEntry(self, index):
        self.read_entries.append(index)
        self.index = index
        return 1

    def GetCurrentFile(self):
        return NS(GetName=lambda: "fixture.root")

    def GetTree(self):
        return NS(GetReadEntry=lambda: self.index)

    def __getattr__(self, name):
        return getattr(self.events[self.index], name)


def vector(x=0.0, y=0.0, z=0.0):
    return NS(x=x, y=y, z=z)


def cluster(distance, eop=1.0):
    return NS(position=vector(1000 + distance, 0, 100), energy=10.0 * eop)


def base_event():
    electron = NS(PDG=11, generatorStatus=1, parents_begin=0, parents_end=1,
                  momentum=vector(10.0), mass=0.000511)
    parent = NS(PDG=11, generatorStatus=23, parents_begin=0, parents_end=0,
                momentum=vector(10.0), mass=0.000511)
    reco = NS(momentum=vector(10.0), mass=0.000511, tracks_begin=0, tracks_end=1)
    return NS(
        MCParticles=[electron, parent], _MCParticles_parents=[NS(index=1)],
        ReconstructedParticles=[reco], _ReconstructedParticles_tracks=[NS(collectionID=7, index=3)],
        _ReconstructedParticleAssociations_rec=[NS(index=0)],
        _ReconstructedParticleAssociations_sim=[NS(index=0)],
        CalorimeterTrackProjections=[NS(points_begin=0, points_end=1)],
        _CalorimeterTrackProjections_track=[NS(collectionID=7, index=3)],
        _CalorimeterTrackProjections_points=[NS(position=vector(1000, 0, 100), system=101, surface=1)],
        EcalBarrelClusters=[cluster(10), cluster(30)],
    )


def selection_cases():
    cases = [base_event() for _ in range(16)]
    cases[1].MCParticles[0].PDG = 211  # no scattered electron
    cases[2].MCParticles[1].generatorStatus = 2  # excluded first parent
    cases[3]._ReconstructedParticleAssociations_sim = []  # no reco
    cases[4].ReconstructedParticles[0].tracks_end = 0  # no track
    cases[5]._CalorimeterTrackProjections_track[0].collectionID = 8  # same index, different collection
    cases[6]._CalorimeterTrackProjections_points[0].system = 102  # no ECAL point
    cases[7].EcalBarrelClusters = []
    cases[8].EcalBarrelClusters = [cluster(1, 0.7), cluster(2, 1.3)]  # strict E/p endpoints
    cases[9].EcalBarrelClusters = [cluster(87)]  # strict distance endpoint
    cases[10].EcalBarrelClusters = [cluster(1, 0.5), cluster(30, 1.0)]  # E/p first, then nearest
    cases[11].ReconstructedParticles[0].momentum = vector()  # no division by zero
    cases[12].MCParticles[0].parents_end = 0  # missing parent passes the legacy veto
    cases[13]._MCParticles_parents[0].index = 99  # invalid parent index also passes
    cases[14].EcalBarrelClusters = [cluster(250, 1.0)]  # histogram overflow
    cases[15].EcalBarrelClusters = [cluster(15, 0.9), cluster(15, 1.1)]  # tie retains first
    return cases
