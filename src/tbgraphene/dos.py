"""Density of states: numerical (k-mesh histogram) and analytic NN reference.

DOS is per spin per unit cell (2 atoms), normalised so that integral over E = 2 (two bands).
"""
from __future__ import annotations

import numpy as np
from scipy.special import ellipk

from .graphene import Graphene


def dos_numeric(gr: Graphene, nk: int, edges: np.ndarray, chunk: int = 400) -> np.ndarray:
    """Histogram DOS on an nk x nk mesh of the rhombic BZ, evaluated chunk-wise."""
    counts = np.zeros(len(edges) - 1)
    s = (np.arange(nk) + 0.5) / nk
    for i0 in range(0, nk, chunk):
        u = s[i0:i0 + chunk]
        U, V = np.meshgrid(u, s, indexing="ij")
        k = U[..., None] * gr.b1 + V[..., None] * gr.b2
        e = gr.bands(k).reshape(-1)
        counts += np.histogram(e, bins=edges)[0]
    dE = np.diff(edges)
    return counts / (nk * nk) / dE


def dos_analytic_nn(E: np.ndarray, t: float) -> np.ndarray:
    """Hobson & Nierenberg (1953) / Castro Neto et al. RMP 2009 Eq. (14), NN only, t'=0.

    Castro Neto's expression counts spin (and gives 4 states/cell over the band);
    here it is divided by 2 to give per-spin DOS (integral 2).
    """
    x = np.abs(np.asarray(E, dtype=float)) / t
    rho = np.zeros_like(x)
    inner = (x > 0) & (x < 1)
    outer = (x > 1) & (x < 3)
    F = lambda y: (1 + y) ** 2 - (y ** 2 - 1) ** 2 / 4.0
    z0 = np.where(inner, F(x), np.where(outer, 4 * x, 1.0))
    z1 = np.where(inner, 4 * x, np.where(outer, F(x), 0.0))
    m = np.clip(z1 / z0, 0, 1 - 1e-16)          # scipy ellipk takes parameter m = k^2
    val = (4.0 / np.pi ** 2) * x / t / np.sqrt(z0) * ellipk(m)
    rho = np.where(inner | outer, val, 0.0)
    return rho / 2.0
