"""Optical properties of skin and tattoo ink (mm, 1/mm, nm). Scattering and melanin after Jacques 2013."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

GREEN, RED, IR = 530.0, 660.0, 940.0
WAVELENGTHS = (GREEN, RED, IR)

# molar extinction (cm^-1 / M): oxy- and deoxy-haemoglobin
_EPS = {GREEN: (39956.0, 39036.0), RED: (319.6, 3226.56), IR: (1214.0, 693.44)}
_HB_MOL = 150.0 / 64500.0                      # mol/L
_WATER = {GREEN: 0.00004, RED: 0.00036, IR: 0.027}   # mm^-1 (approximate)


def reduced_scattering(wl: float) -> float:
    """Generic skin mu_s' in mm^-1."""
    return 4.6 * (wl / 500.0) ** -1.421


def scattering(wl: float, g: float = 0.9) -> float:
    """Scattering coefficient mu_s in mm^-1 for anisotropy g."""
    return reduced_scattering(wl) / (1.0 - g)


def melanin(wl: float) -> float:
    """mu_a of melanosome interior, mm^-1 (multiply by the epidermal melanosome volume fraction)."""
    return 51.9 * (wl / 500.0) ** -3.5


def blood(wl: float, spo2: float) -> float:
    """mu_a of whole blood, mm^-1, at oxygen saturation spo2 (0-1)."""
    if wl not in _EPS:
        raise ValueError(f"haemoglobin data only for {sorted(_EPS)} nm")
    e_oxy, e_deoxy = _EPS[wl]
    eps = spo2 * e_oxy + (1.0 - spo2) * e_deoxy
    return 2.303 * _HB_MOL * eps / 10.0        # cm^-1 -> mm^-1


def water(wl: float) -> float:
    return _WATER[wl]


@dataclass(frozen=True)
class Ink:
    """A parametric tattoo ink."""

    name: str
    peak: float | None          # nm; None = broadband carbon black
    width: float = 45.0         # nm (Gaussian sigma)
    baseline: float = 0.05      # fraction of the peak absorbed everywhere

    def mu_a(self, wl: float, density: float) -> float:
        if self.peak is None:
            return density * (wl / 500.0) ** -1.0
        return density * (self.baseline + np.exp(-0.5 * ((wl - self.peak) / self.width) ** 2))


INKS = {
    "black": Ink("black", None),
    "red": Ink("red", 530.0),       # absorbs green light (looks red)
    "green": Ink("green", 650.0),   # absorbs red light (looks green)
    "blue": Ink("blue", 610.0),     # absorbs orange/red light (looks blue)
}
