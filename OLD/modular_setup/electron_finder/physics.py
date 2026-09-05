"""Kinematics in GeV and position geometry in mm (independent of ROOT)."""

import math


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
