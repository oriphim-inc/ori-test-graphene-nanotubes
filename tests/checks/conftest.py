"""Sceptic's acceptance checks: shared fixtures and an INDEPENDENT graphene implementation.

Target directory: results/ (authoritative). Override for sanity tests with
    CHECKS_RESULTS_DIR=results_pilot .venv/bin/python -m pytest tests/checks -q
Everything skips if the target directory has no summary.json.

The helpers below deliberately do not import src/: lattice vectors, f(k) (in a different
gauge, f = 1 + e^{-ik.a1} + e^{-ik.a2}) and constants are re-derived here from PLAN.md, so
agreement with the run tests the code, not the code against itself.
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / os.environ.get("CHECKS_RESULTS_DIR", "results")
PILOT = ROOT / "results_pilot"
AUTHORITATIVE = RES.resolve() == (ROOT / "results").resolve()

# Constants from PLAN.md (user / Castro Neto 2009 / CODATA 2018), not from src/.
T = 2.7
ACC = 1.42
A = np.sqrt(3) * ACC
HBAR = 6.582119569e-16      # eV s
E_CHARGE = 1.602176634e-19  # C
M_E = 9.1093837015e-31      # kg
HV = 1.5 * T * ACC          # hbar v_F in eV*angstrom (analytic NN result)
VF = HV * 1e-10 / HBAR      # m/s

A1 = ACC / 2 * np.array([3.0, np.sqrt(3)])
A2 = ACC / 2 * np.array([3.0, -np.sqrt(3)])
B1 = 2 * np.pi / (3 * ACC) * np.array([1.0, np.sqrt(3)])
B2 = 2 * np.pi / (3 * ACC) * np.array([1.0, -np.sqrt(3)])
K_PT = np.array([2 * np.pi / (3 * ACC), 2 * np.pi / (3 * np.sqrt(3) * ACC)])
KP_PT = np.array([2 * np.pi / (3 * ACC), -2 * np.pi / (3 * np.sqrt(3) * ACC)])
GAMMA = np.zeros(2)
M_PT = np.array([2 * np.pi / (3 * ACC), 0.0])
CELL_AREA = abs(A1[0] * A2[1] - A1[1] * A2[0])


def f_abs(k):
    k = np.asarray(k, float)
    return np.abs(1 + np.exp(-1j * (k @ A1)) + np.exp(-1j * (k @ A2)))


def bands(k, tp=0.0, delta=0.0):
    fa = f_abs(k)
    e0 = -tp * (fa ** 2 - 3)
    r = np.sqrt(delta ** 2 + (T * fa) ** 2)
    return np.stack([e0 - r, e0 + r], -1)


def equiv_mod_G(k, k0, tol=1e-6):
    c = np.linalg.solve(np.column_stack([B1, B2]), np.asarray(k) - k0)
    return bool(np.all(np.abs(c - np.round(c)) < tol))


def tube_geometry(n, m):
    from math import gcd
    L = A * np.sqrt(n * n + n * m + m * m)
    dR = gcd(2 * m + n, 2 * n + m)
    N = 2 * (n * n + n * m + m * m) // dR
    return dict(L=L, d_nm=L / np.pi / 10, theta=np.degrees(np.arctan2(np.sqrt(3) * m, 2 * n + m)),
                N=N, Tlen=np.sqrt(3) * L / dR)


def tube_min_f(n, m):
    """Independent zone folding: min |f| over all k with k.C_h = 2 pi mu (mu integer).

    Lines are parametrised as k = (2 pi mu / L) C^ + s T^ over a region larger than the BZ,
    so no cutting-line bookkeeping (K1, K2, N) from src/ is reused. Dense scan + bounded
    refinement around the best few samples.
    """
    from scipy.optimize import minimize_scalar
    Ch = n * A1 + m * A2
    L = np.linalg.norm(Ch)
    c = Ch / L
    tv = np.array([-c[1], c[0]])
    kmax = 2.0 * np.linalg.norm(K_PT)
    mus = np.arange(-int(kmax * L / (2 * np.pi)) - 1, int(kmax * L / (2 * np.pi)) + 2)
    s = np.linspace(-kmax, kmax, 40001)
    ds = s[1] - s[0]
    best = np.inf
    for mu in mus:
        k = (2 * np.pi * mu / L) * c + s[:, None] * tv
        fa = f_abs(k)
        for i in np.argsort(fa)[:4]:
            res = minimize_scalar(lambda x: f_abs((2 * np.pi * mu / L) * c + x * tv),
                                  bounds=(s[i] - ds, s[i] + ds), method="bounded",
                                  options={"xatol": 1e-13})
            best = min(best, res.fun, fa[i])
    return best


def zigzag_gap(n):
    """Exact zone-folding gap of (n,0): 2t min_q |1 + 2 cos(q pi / n)| (direct, at k=0)."""
    q = np.arange(2 * n)
    return 2 * T * np.min(np.abs(1 + 2 * np.cos(q * np.pi / n)))


def tubes_in_window(dmin=0.5, dmax=2.5):
    out = []
    for n in range(1, 80):
        for m in range(0, n + 1):
            d = A * np.sqrt(n * n + n * m + m * m) / np.pi / 10
            if dmin <= d <= dmax:
                out.append((n, m))
    return out


def _need(path):
    if not path.exists():
        pytest.skip(f"{path.relative_to(ROOT)} not present yet")
    return path


@pytest.fixture(scope="session")
def S():
    return json.loads(_need(RES / "summary.json").read_text())


@pytest.fixture(scope="session")
def S_pilot():
    if AUTHORITATIVE is False:
        pytest.skip("convergence vs pilot only checked for results/")
    return json.loads(_need(PILOT / "summary.json").read_text())


@pytest.fixture(scope="session")
def dos_csv():
    return np.loadtxt(_need(RES / "dos.csv"), delimiter=",", skiprows=1)


@pytest.fixture(scope="session")
def tubes_csv():
    with open(_need(RES / "tubes.csv")) as f:
        return list(csv.DictReader(f))


@pytest.fixture(scope="session")
def sweep_csv():
    with open(_need(RES / "gap_sweep.csv")) as f:
        return list(csv.DictReader(f))
