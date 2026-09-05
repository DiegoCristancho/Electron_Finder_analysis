"""A minimal histogram recorder for selection parity tests without PyROOT.

This verifies histogram fills and fractions, not ROOT serialization or drawing.
"""

from copy import deepcopy
from types import SimpleNamespace


class Histogram:
    def __init__(self, name, title, bins, low, high):
        self.name, self.title = name, title
        self.bins, self.low, self.high = bins, low, high
        self.contents = [0.0] * (bins + 2)
        self.entries = 0
        self.fills = []

    def SetDirectory(self, *_):
        pass

    def FindBin(self, value):
        if value < self.low:
            return 0
        if value >= self.high:
            return self.bins + 1
        return int((value - self.low) / (self.high - self.low) * self.bins) + 1

    def Fill(self, value):
        self.contents[self.FindBin(value)] += 1
        self.entries += 1
        self.fills.append(value)

    def GetBinCenter(self, index):
        return self.low + (index - 0.5) * (self.high - self.low) / self.bins

    def GetBinContent(self, index):
        return self.contents[index]

    def GetNbinsX(self):
        return self.bins

    def GetEntries(self):
        return self.entries

    def GetMean(self):
        # Match ROOT's default exclusion of underflow/overflow from statistics.
        values = [value for value in self.fills if self.low <= value < self.high]
        return sum(values) / len(values) if values else 0.0

    def GetMaximumBin(self):
        return max(range(1, self.bins + 1), key=lambda i: self.contents[i])

    def Clone(self, name):
        clone = deepcopy(self)
        clone.name = name
        return clone

    def SetTitle(self, title):
        self.title = title

    def Divide(self, denominator):
        self.contents = [num / den if den else 0.0 for num, den in zip(self.contents, denominator.contents)]


ROOT = SimpleNamespace(
    TH1F=Histogram,
    PyConfig=SimpleNamespace(IgnoreCommandLineOptions=True),
    gROOT=SimpleNamespace(SetBatch=lambda _: None),
    gStyle=SimpleNamespace(SetOptStat=lambda _: None),
)
