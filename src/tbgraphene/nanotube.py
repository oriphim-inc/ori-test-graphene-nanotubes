"""Zone-folded single-wall carbon nanotubes (Saito, Dresselhaus & Dresselhaus 1998, ch. 3-4).

Two independent routes to the tube bands:
  1. cutting lines in the graphene BZ: k = mu K1 + kappa K2/|K2|;
  2. an explicit real-space translational unit cell (2N atoms on the unrolled sheet, wrapped
     periodically along C_h), whose 2N x 2N Bloch Hamiltonian is diagonalised directly.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import gcd

import numpy as np
from scipy.optimize import minimize_scalar

from .graphene import Graphene


@dataclass
class Tube:
    n: int
    m: int
    gr: Graphene

    def __post_init__(self):
        n, m, gr = self.n, self.m, self.gr
        self.Ch = n * gr.a1 + m * gr.a2
        self.L = np.linalg.norm(self.Ch)                 # circumference, angstrom
        self.d = self.L / np.pi                          # diameter, angstrom
        # chiral angle w.r.t. a1 (zigzag = 0, armchair = 30 deg)
        self.theta = np.degrees(np.arccos(self.Ch @ gr.a1 / (self.L * np.linalg.norm(gr.a1))))
        self.dR = gcd(2 * m + n, 2 * n + m)
        self.t1 = (2 * m + n) // self.dR
        self.t2 = -(2 * n + m) // self.dR
        self.T = self.t1 * gr.a1 + self.t2 * gr.a2
        self.Tlen = np.linalg.norm(self.T)
        self.N = 2 * (n * n + n * m + m * m) // self.dR  # hexagons per cell
        self.K1 = (-self.t2 * gr.b1 + self.t1 * gr.b2) / self.N
        self.K2 = (m * gr.b1 - n * gr.b2) / self.N
        self.K2hat = self.K2 / np.linalg.norm(self.K2)

    # ---- geometry checks -----------------------------------------------
    def k_on_cutting_line(self) -> bool:
        """True if graphene's K point lies on a cutting line (k.C_h = 2 pi * integer)."""
        K = self.gr.high_symmetry()["K"]
        x = K @ self.Ch / (2 * np.pi)
        return abs(x - round(x)) < 1e-9

    # ---- route 1: cutting lines ----------------------------------------
    def kvec(self, mu, kappa):
        mu = np.asarray(mu, dtype=float)
        kappa = np.asarray(kappa, dtype=float)
        return mu[..., None] * self.K1 + kappa[..., None] * self.K2hat

    def bands_cutting(self, kappa: np.ndarray) -> np.ndarray:
        """All 2N subbands at 1D wavevectors kappa, shape (len(kappa), 2N) sorted."""
        mu = np.arange(self.N)
        MU, KA = np.meshgrid(mu, kappa, indexing="ij")
        e = self.gr.bands(self.kvec(MU, KA))            # (N, nk, 2)
        e = np.transpose(e, (1, 0, 2)).reshape(len(kappa), -1)
        return np.sort(e, axis=1)

    def gap_cutting(self, nk: int = 2001):
        """Global gap (min E_c - max E_v) over all cutting lines, with bounded refinement."""
        kmax = np.pi / self.Tlen
        kap = np.linspace(-kmax, kmax, nk)
        mu = np.arange(self.N)
        MU, KA = np.meshgrid(mu, kap, indexing="ij")
        e = self.gr.bands(self.kvec(MU, KA))            # (N, nk, 2)
        ev, ec = e[..., 0], e[..., 1]
        dk = kap[1] - kap[0]

        def refine(arr, sign, band):
            i_mu, i_k = np.unravel_index(np.argmax(sign * arr), arr.shape)
            k0 = kap[i_k]
            res = minimize_scalar(
                lambda q: -sign * self.gr.bands(self.kvec(np.array(i_mu), np.array(q)))[band],
                bounds=(max(-kmax, k0 - dk), min(kmax, k0 + dk)), method="bounded",
                options={"xatol": 1e-12})
            val = -sign * res.fun
            best = sign * max(sign * val, sign * arr[i_mu, i_k])
            return best, int(i_mu), float(res.x)

        vbm, mu_v, k_v = refine(ev, +1, 0)
        cbm, mu_c, k_c = refine(ec, -1, 1)
        # raw value: negative would mean band overlap (semimetal), not clipped
        return cbm - vbm, dict(vbm=vbm, cbm=cbm, mu_v=mu_v, k_v=k_v, mu_c=mu_c, k_c=k_c)

    # ---- route 2: explicit real-space unit cell ------------------------
    def build_cell(self):
        gr = self.gr
        B = np.column_stack([self.Ch, self.T])           # columns: C_h, T
        Binv = np.linalg.inv(B)
        corners = np.array([[0, 0], [self.n, self.m], [self.t1, self.t2],
                            [self.n + self.t1, self.m + self.t2]])
        lo, hi = corners.min(0) - 2, corners.max(0) + 2
        I, J = np.meshgrid(np.arange(lo[0], hi[0] + 1), np.arange(lo[1], hi[1] + 1), indexing="ij")
        R = I.reshape(-1, 1) * gr.a1 + J.reshape(-1, 1) * gr.a2
        offB = gr.nn_vectors[0]
        eps = 1e-9
        sites = []
        for sub, off in ((0, np.zeros(2)), (1, offB)):
            fr = (R + off) @ Binv.T
            keep = np.all((fr >= -eps) & (fr < 1 - eps), axis=1)
            for p in fr[keep]:
                sites.append((sub, p))
        self.sub = np.array([s for s, _ in sites])
        self.frac = np.array([p for _, p in sites])
        self.pos = self.frac @ B.T
        return self.sub, self.frac

    def _index(self, frac_wrapped):
        d = self.frac - frac_wrapped
        d -= np.round(d)
        i = np.argmin(np.sum(d ** 2, axis=1))
        if np.sum(d[i] ** 2) > 1e-14:
            raise RuntimeError("neighbour not found in unit cell")
        return i

    def _hoppings(self):
        """List of (i, j, shift_along_T, amplitude)."""
        if not hasattr(self, "frac"):
            self.build_cell()
        gr = self.gr
        Binv = np.linalg.inv(np.column_stack([self.Ch, self.T]))
        hops = []
        for i, (s, r) in enumerate(zip(self.sub, self.pos)):
            nn = gr.nn_vectors if s == 0 else -gr.nn_vectors
            for v in nn:
                fr = (r + v) @ Binv.T
                shift = int(np.floor(fr[1] + 1e-9))
                fw = fr - np.floor(fr + 1e-9)
                hops.append((i, self._index(fw), shift, -gr.t))
            if gr.tp != 0.0:
                for v in gr.nnn_vectors:
                    fr = (r + v) @ Binv.T
                    shift = int(np.floor(fr[1] + 1e-9))
                    fw = fr - np.floor(fr + 1e-9)
                    hops.append((i, self._index(fw), shift, -gr.tp))
        return hops

    def _blocks(self):
        """H(k) = sum_s H_s exp(i k s |T|): dict shift -> real matrix."""
        if not hasattr(self, "_H"):
            hops = self._hoppings()
            nat = len(self.sub)
            Hs = {}
            for i, j, s, amp in hops:
                Hs.setdefault(s, np.zeros((nat, nat)))[i, j] += amp
            Hs.setdefault(0, np.zeros((nat, nat)))
            Hs[0] = Hs[0] + np.diag(np.where(self.sub == 0, self.gr.delta, -self.gr.delta))
            self._H = Hs
        return self._H

    def bands_realspace(self, kappa: np.ndarray) -> np.ndarray:
        Hs = self._blocks()
        kappa = np.atleast_1d(np.asarray(kappa, dtype=float))
        H = sum(M[None] * np.exp(1j * kappa * s * self.Tlen)[:, None, None] for s, M in Hs.items())
        return np.linalg.eigvalsh(H)

    def gap_realspace(self, nk: int = 401):
        """Gap from the explicit-cell Hamiltonian: grid search + bounded refinement."""
        kmax = np.pi / self.Tlen
        kap = np.linspace(-kmax, kmax, nk)
        e = self.bands_realspace(kap)
        h = e.shape[1] // 2
        dk = kap[1] - kap[0]

        def ext(band, sign):
            i = np.argmax(sign * e[:, band])
            res = minimize_scalar(lambda q: -sign * self.bands_realspace(q)[0, band],
                                  bounds=(max(-kmax, kap[i] - dk), min(kmax, kap[i] + dk)),
                                  method="bounded", options={"xatol": 1e-12})
            return sign * max(sign * (-sign * res.fun), sign * e[i, band])

        return ext(h, -1) - ext(h - 1, +1)


def tubes_in_range(dmin_nm: float, dmax_nm: float, gr: Graphene):
    out = []
    for n in range(1, 60):
        for m in range(0, n + 1):
            d = gr.a * np.sqrt(n * n + n * m + m * m) / np.pi / 10.0
            if dmin_nm <= d <= dmax_nm:
                out.append((n, m))
    return out
