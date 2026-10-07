"""White Monte Carlo photon transport in a scattering slab; absorption is applied afterwards by reweighting."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

DEPTH_EDGES = np.concatenate([np.arange(0.0, 3.0 + 1e-9, 0.05), [8.0]])   # mm; 60 fine bins + deep bin


@dataclass
class MCResult:
    """Exiting photons of one transport run."""

    rho: np.ndarray            # (M,) exit radius, mm
    path: np.ndarray           # (M, n_bins) path length per depth bin, mm
    n_launched: int
    mus: float
    g: float
    edges: np.ndarray = field(default_factory=lambda: DEPTH_EDGES.copy())

    def save(self, path: str) -> None:
        np.savez_compressed(path, rho=self.rho, path=self.path, n_launched=self.n_launched, mus=self.mus,
                            g=self.g, edges=self.edges)

    @classmethod
    def load(cls, path: str) -> "MCResult":
        z = np.load(path)
        return cls(z["rho"], z["path"], int(z["n_launched"]), float(z["mus"]), float(z["g"]), z["edges"])


def _fresnel(cos_i: np.ndarray, n1: float, n2: float) -> np.ndarray:
    """Unpolarised Fresnel reflectance from medium n1 into n2 at incidence cosine cos_i."""
    sin_t = n1 / n2 * np.sqrt(np.clip(1.0 - cos_i ** 2, 0.0, 1.0))
    tir = sin_t >= 1.0
    cos_t = np.sqrt(np.clip(1.0 - sin_t ** 2, 0.0, 1.0))
    rs = ((n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)) ** 2
    rp = ((n1 * cos_t - n2 * cos_i) / (n1 * cos_t + n2 * cos_i)) ** 2
    return np.where(tir, 1.0, 0.5 * (rs + rp))


def _hg_cos(g: float, xi: np.ndarray) -> np.ndarray:
    if g == 0:
        return 2 * xi - 1
    return (1 + g * g - ((1 - g * g) / (1 - g + 2 * g * xi)) ** 2) / (2 * g)


def _scatter(d: np.ndarray, cos_t: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Rotate unit directions d (N, 3) by polar cosine cos_t and azimuth phi (MCML formulas)."""
    sin_t = np.sqrt(np.clip(1 - cos_t ** 2, 0, 1))
    cp, sp = np.cos(phi), np.sin(phi)
    ux, uy, uz = d[:, 0], d[:, 1], d[:, 2]
    out = np.empty_like(d)
    vert = np.abs(uz) > 0.99999
    nv = ~vert
    den = np.sqrt(1 - uz[nv] ** 2)
    out[nv, 0] = sin_t[nv] * (ux[nv] * uz[nv] * cp[nv] - uy[nv] * sp[nv]) / den + ux[nv] * cos_t[nv]
    out[nv, 1] = sin_t[nv] * (uy[nv] * uz[nv] * cp[nv] + ux[nv] * sp[nv]) / den + uy[nv] * cos_t[nv]
    out[nv, 2] = -sin_t[nv] * cp[nv] * den + uz[nv] * cos_t[nv]
    out[vert, 0] = sin_t[vert] * cp[vert]
    out[vert, 1] = sin_t[vert] * sp[vert]
    out[vert, 2] = np.sign(uz[vert]) * cos_t[vert]
    return out / np.linalg.norm(out, axis=1, keepdims=True)


def _add_path(acc: np.ndarray, rows: np.ndarray, z0: np.ndarray, z1: np.ndarray, s: np.ndarray,
              edges: np.ndarray) -> None:
    """Add segment lengths s (from depth z0 to z1) to the depth bins of acc[rows]."""
    nb = len(edges) - 1
    b0 = np.clip(np.searchsorted(edges, z0, side="right") - 1, 0, nb - 1)
    b1 = np.clip(np.searchsorted(edges, z1, side="right") - 1, 0, nb - 1)
    same = b0 == b1
    np.add.at(acc, (rows[same], b0[same]), s[same])          # most steps stay inside one 50 um bin
    cross = np.flatnonzero(~same)
    if cross.size:
        lo = np.minimum(z0[cross], z1[cross])[:, None]
        hi = np.maximum(z0[cross], z1[cross])[:, None]
        overlap = np.clip(np.minimum(hi, edges[None, 1:]) - np.maximum(lo, edges[None, :-1]), 0, None)
        acc[rows[cross]] += (s[cross, None] * overlap / (hi - lo)).astype(acc.dtype)


def run(mus: float, g: float = 0.9, n_photons: int = 100_000, n_tissue: float = 1.4, zmax: float = 8.0,
        rho_max: float = 8.0, max_steps: int = 6000, batch: int = 50_000, seed: int = 0,
        edges: np.ndarray = DEPTH_EDGES) -> MCResult:
    """Transport n_photons through a non-absorbing slab of thickness zmax (mm) with scattering mus (mm^-1)."""
    rng = np.random.default_rng(seed)
    rhos, paths = [], []
    nb = len(edges) - 1
    for start in range(0, n_photons, batch):
        n = min(batch, n_photons - start)
        pos = np.zeros((n, 3))
        d = np.tile([0.0, 0.0, 1.0], (n, 1))           # z points into the tissue
        acc = np.zeros((n, nb), np.float32)
        active = np.arange(n)
        for _ in range(max_steps):
            if active.size == 0:
                break
            p, u = pos[active], d[active]
            s = -np.log(rng.random(active.size)) / mus
            new = p + s[:, None] * u
            out_top = new[:, 2] < 0
            out_bot = new[:, 2] > zmax
            # truncate steps that hit a boundary at the boundary
            s_eff = s.copy()
            s_eff[out_top] = -p[out_top, 2] / u[out_top, 2]
            s_eff[out_bot] = (zmax - p[out_bot, 2]) / u[out_bot, 2]
            new = p + s_eff[:, None] * u
            _add_path(acc, active, p[:, 2], new[:, 2], s_eff, edges)
            pos[active] = new
            # top surface: Fresnel reflection or escape
            top = np.flatnonzero(out_top)
            escaped = np.zeros(active.size, bool)
            if top.size:
                refl = rng.random(top.size) < _fresnel(np.abs(u[top, 2]), n_tissue, 1.0)
                esc = top[~refl]
                escaped[esc] = True
                ref = top[refl]
                d[active[ref], 2] *= -1
                pos[active[ref], 2] = 0.0
            done = escaped | out_bot
            # scatter the photons that are still inside and did not just reflect
            inside = ~(done | out_top)
            idx = active[inside]
            if idx.size:
                cos_t = _hg_cos(g, rng.random(idx.size))
                d[idx] = _scatter(d[idx], cos_t, 2 * np.pi * rng.random(idx.size))
            ex = active[escaped]
            r = np.hypot(pos[ex, 0], pos[ex, 1])
            keep = r <= rho_max
            rhos.append(r[keep])
            paths.append(acc[ex][keep])
            active = active[~done]
    return MCResult(np.concatenate(rhos), np.concatenate(paths), n_photons, mus, g, edges)
