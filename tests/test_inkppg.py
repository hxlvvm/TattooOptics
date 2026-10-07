"""Tests of the transport code and the ink physics."""
import numpy as np
import pytest

from inkppg import Skin, measure, mc, optics, ratio_of_ratios
from inkppg.diffusion import reflectance


@pytest.fixture(scope="module")
def runs():
    return {wl: mc.run(optics.scattering(wl), n_photons=30_000, seed=7) for wl in optics.WAVELENGTHS}


def test_fresnel_normal_incidence():
    assert mc._fresnel(np.array([1.0]), 1.4, 1.0)[0] == pytest.approx((0.4 / 2.4) ** 2, rel=1e-6)
    assert mc._fresnel(np.array([0.5]), 1.4, 1.0)[0] == 1.0          # beyond the critical angle


def test_henyey_greenstein_mean_cosine():
    xi = np.random.default_rng(0).random(200_000)
    assert mc._hg_cos(0.9, xi).mean() == pytest.approx(0.9, abs=0.005)


def test_photon_bookkeeping(runs):
    r = runs[optics.RED]
    assert 0.5 * r.n_launched < len(r.rho) <= r.n_launched
    assert np.all(r.path >= 0) and np.all(r.path.sum(1) > 0)


def test_monte_carlo_matches_diffusion_theory(runs):
    """White-MC reweighting of a homogeneous absorber vs diffusion theory (Kienle & Patterson) at 1-4 mm."""
    wl, mua = optics.IR, 0.01
    r = runs[wl]
    musp = optics.reduced_scattering(wl)
    w = np.exp(-mua * r.path.sum(1))
    for rho in (1.0, 2.0, 3.0, 4.0):
        sel = np.abs(r.rho - rho) <= 0.25
        area = np.pi * ((rho + 0.25) ** 2 - (rho - 0.25) ** 2)
        mc_r = w[sel].sum() / r.n_launched / area
        assert mc_r == pytest.approx(reflectance(rho, mua, musp), rel=0.15)


def test_black_ink_darkens_every_wavelength(runs):
    for wl in optics.WAVELENGTHS:
        clean = measure(runs[wl], Skin(), wl, 2.0)
        inked = measure(runs[wl], Skin().with_ink("black", 2.0), wl, 2.0)
        assert inked.reflectance < 0.7 * clean.reflectance and inked.snr < clean.snr


def test_red_ink_hits_green_light_much_more_than_infrared(runs):
    loss = {wl: measure(runs[wl], Skin().with_ink("red", 2.0), wl, 2.0).reflectance
            / measure(runs[wl], Skin(), wl, 2.0).reflectance for wl in (optics.GREEN, optics.IR)}
    assert loss[optics.GREEN] < 0.5 * loss[optics.IR]


def test_coloured_ink_shifts_the_oximetry_ratio(runs):
    clean = ratio_of_ratios(runs[optics.RED], runs[optics.IR], Skin(), 3.0)
    inked = ratio_of_ratios(runs[optics.RED], runs[optics.IR], Skin().with_ink("green", 2.0), 3.0)
    assert abs(inked / clean - 1) > 0.01
