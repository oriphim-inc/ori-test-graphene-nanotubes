"""Pi-band tight-binding model of graphene on the honeycomb lattice.

Conventions (Castro Neto et al., Rev. Mod. Phys. 81, 109 (2009), with a sublattice mass added):
    H(k) = [[+Delta, -t f(k)], [-t f(k)^*, -Delta]] - t' g(k) * 1
    f(k) = sum_{3 NN delta} exp(i k.delta),   g(k) = sum_{6 NNN} exp(i k.d) = |f|^2 - 3
    E_pm(k) = -t' g(k) +- sqrt(Delta^2 + t^2 |f(k)|^2)
Lengths in angstrom, energies in eV, wavevectors in 1/angstrom.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

HBAR_EV_S = 6.582119569e-16  # CODATA 2018, eV s
A_CC = 1.42                  # C-C distance, angstrom (Castro Neto 2009; Saito et al. 1998)


@dataclass(frozen=True)
class Graphene:
    t: float = 2.7       # NN hopping, eV (user-specified)
    tp: float = 0.0      # NNN hopping t', eV
    delta: float = 0.0   # sublattice mass, eV (+delta on A, -delta on B)
    acc: float = A_CC

    # ---- lattice -------------------------------------------------------
    @property
    def a(self) -> float:
        return np.sqrt(3.0) * self.acc

    @property
    def a1(self) -> np.ndarray:
        return self.acc / 2.0 * np.array([3.0, np.sqrt(3.0)])

    @property
    def a2(self) -> np.ndarray:
        return self.acc / 2.0 * np.array([3.0, -np.sqrt(3.0)])

    @property
    def nn_vectors(self) -> np.ndarray:
        """Vectors from an A atom to its three B neighbours."""
        c = self.acc
        return np.array([[c / 2, c * np.sqrt(3) / 2], [c / 2, -c * np.sqrt(3) / 2], [-c, 0.0]])

    @property
    def nnn_vectors(self) -> np.ndarray:
        a1, a2 = self.a1, self.a2
        v = np.array([a1, a2, a1 - a2])
        return np.vstack([v, -v])

    @property
    def b1(self) -> np.ndarray:
        return 2 * np.pi / (3 * self.acc) * np.array([1.0, np.sqrt(3.0)])

    @property
    def b2(self) -> np.ndarray:
        return 2 * np.pi / (3 * self.acc) * np.array([1.0, -np.sqrt(3.0)])

    @property
    def cell_area(self) -> float:
        a1, a2 = self.a1, self.a2
        return abs(a1[0] * a2[1] - a1[1] * a2[0])

    def high_symmetry(self) -> dict[str, np.ndarray]:
        c = self.acc
        return {
            "G": np.zeros(2),
            "K": np.array([2 * np.pi / (3 * c), 2 * np.pi / (3 * np.sqrt(3) * c)]),
            "Kp": np.array([2 * np.pi / (3 * c), -2 * np.pi / (3 * np.sqrt(3) * c)]),
            "M": np.array([2 * np.pi / (3 * c), 0.0]),
        }

    # ---- Hamiltonian ---------------------------------------------------
    def f(self, k: np.ndarray) -> np.ndarray:
        k = np.asarray(k, dtype=float)
        ph = k @ self.nn_vectors.T          # (..., 3)
        return np.exp(1j * ph).sum(axis=-1)

    def g(self, k: np.ndarray) -> np.ndarray:
        k = np.asarray(k, dtype=float)
        ph = k @ self.nnn_vectors.T          # (..., 6)
        return np.cos(ph).sum(axis=-1)       # real: +d and -d pairs

    def hamiltonian(self, k: np.ndarray) -> np.ndarray:
        k = np.asarray(k, dtype=float)
        fk = self.f(k)
        gk = self.g(k)
        H = np.zeros(k.shape[:-1] + (2, 2), dtype=complex)
        H[..., 0, 0] = self.delta - self.tp * gk
        H[..., 1, 1] = -self.delta - self.tp * gk
        H[..., 0, 1] = -self.t * fk
        H[..., 1, 0] = -self.t * np.conj(fk)
        return H

    def bands(self, k: np.ndarray) -> np.ndarray:
        """Closed-form eigenvalues, shape (..., 2) sorted (valence, conduction)."""
        fk = self.f(k)
        r = np.sqrt(self.delta ** 2 + self.t ** 2 * np.abs(fk) ** 2)
        e0 = -self.tp * self.g(k)
        return np.stack([e0 - r, e0 + r], axis=-1)

    def bands_numeric(self, k: np.ndarray) -> np.ndarray:
        """Eigenvalues by diagonalising H(k) (independent of the closed form)."""
        return np.linalg.eigvalsh(self.hamiltonian(k))

    # ---- k-space utilities ---------------------------------------------
    def kmesh(self, n: int, shift: bool = True) -> np.ndarray:
        """Uniform n x n Monkhorst-Pack-like mesh over the rhombic BZ, shape (n*n, 2)."""
        s = (np.arange(n) + (0.5 if shift else 0.0)) / n
        u, v = np.meshgrid(s, s, indexing="ij")
        return (u[..., None] * self.b1 + v[..., None] * self.b2).reshape(-1, 2)

    def path(self, labels=("G", "K", "M", "G"), n_per_seg: int = 300):
        hs = self.high_symmetry()
        pts = [hs[l] for l in labels]
        ks, dist, ticks = [], [], [0.0]
        s0 = 0.0
        for i in range(len(pts) - 1):
            p, q = pts[i], pts[i + 1]
            x = np.linspace(0, 1, n_per_seg, endpoint=(i == len(pts) - 2))
            seg = p + x[:, None] * (q - p)
            L = np.linalg.norm(q - p)
            ks.append(seg)
            dist.append(s0 + x * L)
            s0 += L
            ticks.append(s0)
        return np.vstack(ks), np.concatenate(dist), np.array(ticks)


# ---- Dirac points -----------------------------------------------------------
def fold_to_bz(gr: Graphene, k: np.ndarray) -> np.ndarray:
    """Map k to the reciprocal-lattice image closest to Gamma (first BZ, Wigner-Seitz)."""
    best, bestn = k, np.linalg.norm(k)
    for i in range(-2, 3):
        for j in range(-2, 3):
            q = k + i * gr.b1 + j * gr.b2
            nq = np.linalg.norm(q)
            if nq < bestn - 1e-12:
                best, bestn = q, nq
    return best


def find_dirac_points(gr: Graphene, n_seed: int = 24, tol_e: float = 1e-6):
    """Numerically locate zeros of |f(k)| in the BZ (independent of the analytic K)."""
    def obj(k):
        return np.abs(gr.f(k)) ** 2

    seeds = gr.kmesh(n_seed, shift=True)
    vals = obj(seeds)
    order = np.argsort(vals)[: max(12, n_seed)]
    found: list[np.ndarray] = []
    for s in seeds[order]:
        res = minimize(obj, s, method="Nelder-Mead",
                       options={"xatol": 1e-13, "fatol": 1e-28, "maxiter": 4000})
        if gr.t * np.sqrt(res.fun) > tol_e:
            continue
        k = fold_to_bz(gr, res.x)
        # identify equivalence modulo reciprocal lattice
        dup = False
        for q in found:
            d = k - q
            # coefficients in b1,b2 basis
            B = np.column_stack([gr.b1, gr.b2])
            c = np.linalg.solve(B, d)
            if np.allclose(c, np.round(c), atol=1e-6):
                dup = True
                break
        if not dup:
            found.append(k)
    return found


def min_direct_gap(gr: Graphene, n_seed: int = 24):
    """Numerically minimise E_c(k) - E_v(k) over the BZ; returns (gap, k_at_min list)."""
    def obj(k):
        e = gr.bands(k)
        return e[..., 1] - e[..., 0]

    seeds = gr.kmesh(n_seed)
    vals = obj(seeds)
    best = []
    for s in seeds[np.argsort(vals)[:12]]:
        res = minimize(lambda k: obj(k) ** 2, s, method="Nelder-Mead",
                       options={"xatol": 1e-13, "fatol": 1e-30, "maxiter": 4000})
        best.append((float(obj(res.x)), fold_to_bz(gr, res.x)))
    g = min(b[0] for b in best)
    ks = [b[1] for b in best if b[0] < g + 1e-8]
    return g, ks


def fermi_velocity(gr: Graphene, k0: np.ndarray, qmax: float = 2e-3, nq: int = 20,
                   ndir: int = 36):
    """Fit E_c(k0+q) - E_D = hbar v |q| for small q along ndir directions.

    Returns v_F (m/s) mean, per-direction array, and the fitted Dirac energy.
    Uses the conduction band minus Dirac energy; for Delta != 0 this is not meaningful.
    """
    ED = gr.bands(k0).mean()
    qs = np.linspace(qmax / nq, qmax, nq)
    vs = []
    for th in np.linspace(0, 2 * np.pi, ndir, endpoint=False):
        u = np.array([np.cos(th), np.sin(th)])
        kk = k0 + qs[:, None] * u
        e = gr.bands(kk)
        half = (e[:, 1] - e[:, 0]) / 2.0          # symmetric part: hbar v q + O(q^2)
        slope = np.polyfit(qs, half, 2)[1]        # linear coefficient (eV angstrom)
        vs.append(slope * 1e-10 / HBAR_EV_S)
    vs = np.array(vs)
    return vs.mean(), vs, ED
