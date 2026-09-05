"""Settings shared by the analysis and diagnostics (no ROOT dependency)."""

from dataclasses import asdict, dataclass, fields
import json
import math
from pathlib import Path


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
