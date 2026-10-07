"""Wrist PPG with and without tattoo ink, from a white Monte Carlo run.

Skin model (depths in mm, all configurable through `Skin`):
  epidermis   0 - 0.1     melanin (melanosome volume fraction f_mel)
  dermis      0.1 - 2.0   blood volume fraction f_blood, water
  ink sheet   ink_top - ink_top + ink_thickness, inside the dermis (pigment adds to the dermal absorption)
  subcutis    2.0 - 8.0   lower blood fraction
The heartbeat changes the dermal blood volume by a fraction `pulse`; to first order the PPG modulation is

    AC/DC = pulse * mu_a,blood * <L_blood>,   <L_blood> = absorption-weighted mean path in the perfused layers

and the shot-noise-limited SNR of the pulse is AC/DC * sqrt(N_detected).
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from . import optics
from .mc import MCResult


@dataclass(frozen=True)
class Skin:
    f_mel: float = 0.03            # epidermal melanosome volume fraction (light skin ~0.01-0.06, dark ~0.2-0.4)
    f_blood: float = 0.02          # dermal blood volume fraction
    f_blood_sub: float = 0.005     # subcutis blood volume fraction
    spo2: float = 0.98
    pulse: float = 0.02            # fractional blood-volume change per heartbeat
    epidermis: float = 0.1
    dermis_bottom: float = 2.0
    ink: str | None = None
    ink_density: float = 0.0       # mm^-1 at the ink's absorption peak
    ink_top: float = 0.4
    ink_thickness: float = 0.2
    background: float = 0.002      # bloodless tissue absorption, mm^-1 (assumption)

    def with_ink(self, ink: str | None, density: float = 0.0) -> "Skin":
        return replace(self, ink=ink, ink_density=density)


def absorption_profile(skin: Skin, wl: float, edges: np.ndarray, blood_scale: float = 1.0):
    """mu_a per depth bin (mm^-1) and the blood-only part (for the pulsatile term)."""
    mid = 0.5 * (edges[:-1] + edges[1:])
    mb = optics.blood(wl, skin.spo2) * blood_scale
    epi = mid < skin.epidermis
    derm = (mid >= skin.epidermis) & (mid < skin.dermis_bottom)
    sub = mid >= skin.dermis_bottom
    blood = np.where(derm, skin.f_blood * mb, 0.0) + np.where(sub, skin.f_blood_sub * mb, 0.0)
    mu = np.full(mid.shape, skin.background) + blood + (derm | sub) * 0.7 * optics.water(wl)
    mu = mu + epi * skin.f_mel * optics.melanin(wl)
    if skin.ink is not None and skin.ink_density > 0:
        lo, hi = skin.ink_top, skin.ink_top + skin.ink_thickness
        frac = np.clip((np.minimum(edges[1:], hi) - np.maximum(edges[:-1], lo)) / np.diff(edges), 0, 1)
        mu = mu + frac * optics.INKS[skin.ink].mu_a(wl, skin.ink_density)
    return mu, blood


@dataclass
class PPG:
    reflectance: float      # detected fraction per mm^2 of detector (DC)
    ac_dc: float            # pulsatile modulation (perfusion index / 100)
    path_blood: float       # mean path through blood-weighted tissue, mm
    n_detected: float       # photons per sample for the given photon budget
    snr: float              # shot-noise-limited SNR of the pulse

    @property
    def snr_db(self) -> float:
        return 20 * np.log10(max(self.snr, 1e-12))


def measure(mcr: MCResult, skin: Skin, wl: float, rho: float, ring: float = 0.5,
            photons: float = 1e11) -> PPG:
    """PPG at source-detector separation rho (mm), detector ring width `ring` (mm).

    photons: emitted photons per sample. The default is an assumption chosen so that light skin without ink
    gives a pulse SNR of the order of 30-40 dB at 530 nm and 2 mm; compare configurations relatively.
    """
    sel = np.abs(mcr.rho - rho) <= ring / 2
    L = mcr.path[sel].astype(np.float64)
    mu, blood = absorption_profile(skin, wl, mcr.edges)
    w = np.exp(-L @ mu)
    area = np.pi * ((rho + ring / 2) ** 2 - max(rho - ring / 2, 0) ** 2)
    refl = w.sum() / mcr.n_launched / area
    # first-order PPG: d ln R = - sum_k d mu_k <L_k>,  d mu = pulse * blood part
    lb = (w[:, None] * L).sum(0) / max(w.sum(), 1e-300)          # mean path per bin
    ac_dc = skin.pulse * float(blood @ lb)
    path_blood = float(blood @ lb / max(blood.max(), 1e-12))
    n_det = photons * refl
    return PPG(refl, ac_dc, path_blood, n_det, ac_dc * np.sqrt(n_det))


def ratio_of_ratios(mc_red: MCResult, mc_ir: MCResult, skin: Skin, rho: float) -> float:
    """Pulse-oximetry ratio R = (AC/DC)_red / (AC/DC)_IR."""
    return measure(mc_red, skin, optics.RED, rho).ac_dc / measure(mc_ir, skin, optics.IR, rho).ac_dc


def spo2_calibration(mc_red: MCResult, mc_ir: MCResult, skin: Skin, rho: float,
                     levels=np.linspace(0.70, 1.00, 31)):
    """Ratio R for each true SpO2 on ink-free skin: the oximeter's calibration curve."""
    base = skin.with_ink(None)
    r = np.array([ratio_of_ratios(mc_red, mc_ir, replace(base, spo2=float(s)), rho) for s in levels])
    return np.asarray(levels), r


def apparent_spo2(r: float, calib) -> float:
    """Invert the calibration curve (R decreases with SpO2)."""
    levels, rr = calib
    order = np.argsort(rr)
    return float(np.interp(r, rr[order], levels[order]))
