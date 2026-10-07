"""Diffusion-theory reflectance (Farrell 1992; Kienle & Patterson 1997) for checking the Monte Carlo."""
from __future__ import annotations

import numpy as np


def internal_reflection_A(n: float) -> float:
    """Boundary parameter A = (1 + r_d) / (1 - r_d) for a tissue/air interface (Groenhuis empirical r_d)."""
    rd = -1.440 / n ** 2 + 0.710 / n + 0.668 + 0.0636 * n
    return (1 + rd) / (1 - rd)


def _parts(rho, mua: float, musp: float, n: float):
    rho = np.asarray(rho, float)
    mut = mua + musp
    z0 = 1.0 / mut
    d = 1.0 / (3 * mut)
    zb = 2 * internal_reflection_A(n) * d
    mueff = np.sqrt(3 * mua * mut)
    r1 = np.sqrt(z0 ** 2 + rho ** 2)
    r2 = np.sqrt((z0 + 2 * zb) ** 2 + rho ** 2)
    phi = (np.exp(-mueff * r1) / r1 - np.exp(-mueff * r2) / r2) / (4 * np.pi * d)
    flux = (z0 * (mueff + 1 / r1) * np.exp(-mueff * r1) / r1 ** 2
            + (z0 + 2 * zb) * (mueff + 1 / r2) * np.exp(-mueff * r2) / r2 ** 2) / (4 * np.pi)
    return phi, flux


def reflectance_flux(rho, mua: float, musp: float, n: float = 1.4):
    """Farrell's diffuse flux J(rho) (mm^-2) at distance rho (mm) from a pencil beam."""
    return _parts(rho, mua, musp, n)[1]


def reflectance(rho, mua: float, musp: float, n: float = 1.4):
    """Escaping reflectance (mm^-2), Kienle & Patterson partial-current form (coefficients for n = 1.4)."""
    if abs(n - 1.4) > 1e-9:
        raise ValueError("the C1/C2 coefficients here are tabulated for n = 1.4")
    phi, flux = _parts(rho, mua, musp, n)
    return 0.118 * phi + 0.306 * flux
