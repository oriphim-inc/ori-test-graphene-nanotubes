"""Acceptance checks for the main run (results/). Written before the main run; tolerances fixed here.

Tolerance philosophy:
  * Quantities that are exact in the model (Dirac energy, 2*Delta gap, 3t' shift, closed form vs
    diagonalisation, cutting-line vs real-space eigenvalues) get near-machine tolerances
    (1e-6 eV or tighter): any larger miss is a bug, not a numerical limitation.
  * Quantities limited by a k-mesh (DOS histogram, mesh band edges) get tolerances derived from
    the mesh/bin size, plus a convergence requirement against the 1/8-size pilot.
  * Physical approximations (E_g d ~ 2 a_cc t) get tolerances from the known size of the
    trigonal-warping correction (order a_cc/d), stated at each check.
"""
from __future__ import annotations

import numpy as np
import pytest

from conftest import (A, ACC, AUTHORITATIVE, CELL_AREA, E_CHARGE, GAMMA, HV, K_PT, KP_PT, M_E,
                      M_PT, T, VF, bands, equiv_mod_G, f_abs, tube_geometry, tube_min_f,
                      tubes_in_window, zigzag_gap)

TUBES = [(5, 5), (9, 0), (10, 0), (8, 4), (13, 0)]
REF_EGD = 2 * ACC * T / 10  # eV nm, 2 a_cc t (Saito 1998; White & Mintmire 1998)


# ---------------------------------------------------------------- run metadata
def test_parameters_as_agreed(S):
    assert S["t_eV"] == 2.7 and S["acc_angstrom"] == 1.42
    if AUTHORITATIVE:
        assert S["size"] == 1.0, "authoritative results/ must be the full-size run"


# ---------------------------------------------------------------- graphene bands
def test_closed_form_matches_diagonalisation(S):
    # Same 2x2 problem two ways: anything above round-off (1e-10 eV) is a bug.
    assert S["path_closed_form_vs_diag_max_diff_eV"] < 1e-10


def test_high_symmetry_energies(S):
    # NN model: Gamma = -+3t, K = 0, M = -+t, exactly (PLAN.md; Castro Neto 2009).
    pe = S["path_energies_eV"]
    np.testing.assert_allclose(pe["G"], [-3 * T, 3 * T], atol=1e-9)
    np.testing.assert_allclose(pe["K"], [0, 0], atol=1e-9)
    np.testing.assert_allclose(pe["M"], [-T, T], atol=1e-9)


def test_trigonal_warping_numbers_independent(S):
    # Recomputed with the independent f(k) (other gauge): must agree to round-off, and
    # E toward M < E toward Gamma at finite q (known warping direction, saddle at M).
    for key, e in S["trigonal_warping_Ec_eV"].items():
        q = float(key.split("=")[1])
        for name, tgt in (("toward_Gamma", GAMMA), ("toward_M", M_PT)):
            u = (tgt - K_PT) / np.linalg.norm(tgt - K_PT)
            assert e[name] == pytest.approx(bands(K_PT + q * u)[1], abs=1e-9)
        assert e["toward_M"] < e["toward_Gamma"]


# ---------------------------------------------------------------- Dirac points
def test_two_inequivalent_dirac_points_at_K_Kprime(S):
    dps = S["dirac_points"]
    assert len(dps) == 2, "exactly two inequivalent Dirac points per BZ"
    hits = []
    for d in dps:
        k = np.array(d["k_inv_angstrom"])
        # PLAN criterion |E| < 1e-6 eV; recompute with independent bands, not only the reported value
        assert d["E_c_minus_E_v_eV"] < 1e-6
        assert np.ptp(bands(k)) < 1e-6
        # location modulo reciprocal lattice (a WS-folded corner may be any of 3 images).
        # 1e-6 eV gap <-> |q| ~ 1e-7 1/A; allow 1e-6 in reduced coordinates.
        hits.append("K" if equiv_mod_G(k, K_PT) else "Kp" if equiv_mod_G(k, KP_PT) else None)
    assert sorted(hits) == ["K", "Kp"], f"Dirac points not at K and K' (got {hits})"


def test_fermi_velocity(S):
    # v_F = 3 t a_cc / (2 hbar) = 8.737e5 m/s. The fit uses |q| <= 2e-3 1/A where curvature
    # corrections are O((q a_cc)^2) ~ 1e-5, so require 1e-4 relative agreement and an
    # isotropy spread below 1e-4 (cone isotropic to that order).
    fv = S["fermi_velocity"]
    assert fv["analytic_3t_acc_over_2hbar_m_s"] == pytest.approx(VF, rel=1e-9)
    assert 8.6e5 < VF < 8.8e5  # units sanity (m/s, not cm/s or A/s)
    for key, v in fv["per_dirac_point"].items():
        assert v["v_mean_m_s"] == pytest.approx(VF, rel=1e-4), key
        assert v["rel_anisotropy"] < 1e-4, key
        assert abs(v["E_D_eV"]) < 1e-6, key


# ---------------------------------------------------------------- DOS
def test_dos_normalisation(S, dos_csv):
    # Two bands, per spin per cell: integral = 2 (histogram over a window containing all states).
    E, dn, da = dos_csv.T
    dE = E[1] - E[0]
    assert S["dos"]["integral"] == pytest.approx(2.0, abs=1e-9)
    assert dn.sum() * dE == pytest.approx(2.0, abs=1e-6)
    # analytic column independently normalised (log singularities: bin-centre sampling, 1%)
    assert da.sum() * dE == pytest.approx(2.0, rel=1e-2)
    # band edges at +-3t: no states outside
    assert np.all(dn[np.abs(E) > 3 * T + dE] == 0)
    # electron-hole symmetric for NN (allow a few states hopping across float-rounded bin edges)
    np.testing.assert_allclose(dn, dn[::-1], atol=1e-3 * dn.max())


def test_dos_matches_analytic(S, dos_csv):
    # Regular-mesh histogram error scales ~1/nk. The 1/8-size pilot (nk=375) had pointwise 42%,
    # 11-bin-smoothed 3.1%, L1 0.16; at nk=3000 expect ~8x smaller. Tolerances with ~2-3x margin:
    # pointwise <= 15%, smoothed (0.11 eV window) <= 2%, L1 <= 0.05 states/cell.
    # Region: 0.3 eV < |E| < 2.9t, more than 0.1 eV from the van Hove points.
    E, dn, da = dos_csv.T
    away = (np.abs(np.abs(E) - T) > 0.1) & (np.abs(E) < 2.9 * T) & (np.abs(E) > 0.3)
    rel = np.abs(dn - da)[away] / da[away]
    assert rel.max() <= 0.15
    assert S["dos"]["max_rel_dev_from_analytic_away_from_vHs"] <= 0.15
    k = np.ones(11) / 11
    sm = np.abs(np.convolve(dn, k, "same") - np.convolve(da, k, "same")) / np.convolve(da, k, "same")
    assert sm[away].max() <= 0.02
    assert np.abs(dn - da).sum() * (E[1] - E[0]) <= 0.05


def test_dos_van_hove_and_linear_onset(S, dos_csv):
    E, dn, _ = dos_csv.T
    dE = E[1] - E[0]
    # van Hove peaks at +-t to within 2 bins
    neg, pos = S["dos"]["vhs_peak_energies_eV"]
    assert neg == pytest.approx(-T, abs=2 * dE) and pos == pytest.approx(T, abs=2 * dE)
    assert E[E > 0][np.argmax(dn[E > 0])] == pytest.approx(T, abs=2 * dE)
    # linear onset: slope A_c / (pi (hbar v_F)^2) per spin per cell, with ANALYTIC v_F.
    # Fit window up to 0.5 eV: warping correction to the slope is O((E/t)^2) ~ 3%; allow 4%.
    slope_ref = CELL_AREA / (np.pi * HV ** 2)
    assert S["dos"]["low_E_slope_numeric_per_eV2"] == pytest.approx(slope_ref, rel=0.04)
    assert S["dos"]["low_E_slope_from_vF_Ac_over_pi_hbar2_vF2"] == pytest.approx(slope_ref, rel=1e-3)
    # DOS vanishes at E=0 (bin straddling 0 must be tiny: < slope * 2 bins)
    assert dn[np.argmin(np.abs(E))] < slope_ref * 2 * dE + 1e-12


# ---------------------------------------------------------------- nanotubes
def _row(rows, n, m):
    r = [x for x in rows if int(x["n"]) == n and int(x["m"]) == m]
    assert len(r) == 1, (n, m)
    return r[0]


@pytest.mark.parametrize("nm", TUBES)
def test_tube_geometry(S, nm):
    r = _row(S["tubes"], *nm)
    g = tube_geometry(*nm)
    assert r["d_nm"] == pytest.approx(g["d_nm"], rel=1e-9)
    assert r["theta_deg"] == pytest.approx(g["theta"], abs=1e-6)
    assert r["N_hex"] == g["N"] and r["atoms"] == 2 * g["N"]
    assert r["T_angstrom"] == pytest.approx(g["Tlen"], rel=1e-9)


@pytest.mark.parametrize("nm", TUBES)
def test_tube_metallicity_and_gap_independent(S, nm):
    # Metallic <=> (n-m) mod 3 == 0 <=> K on a cutting line. Gap (t'=0, Delta=0) = 2t min|f| over
    # allowed k, recomputed by an independent line scan; agreement to 1e-5 eV (refinement-limited).
    n, m = nm
    r = _row(S["tubes"], n, m)
    metal = (n - m) % 3 == 0
    assert r["metallic"] is metal and r["K_on_cutting_line"] is metal
    gap_ind = 2 * T * tube_min_f(n, m)
    if metal:
        assert gap_ind < 1e-5
        assert abs(r["gap_raw_eV"]) < 1e-5 and abs(r["gap_realspace_eV"]) < 1e-5
        assert r["gap_eV"] == 0.0
    else:
        assert r["gap_eV"] == pytest.approx(gap_ind, abs=1e-5)
        assert r["gap_eV"] > 0.3  # sub-nm semiconducting tubes: gap ~ 0.7-1 eV, not meV
    if m == 0:
        expect = 0.0 if metal else zigzag_gap(n)
        assert r["gap_raw_eV"] == pytest.approx(expect, abs=1e-5)


@pytest.mark.parametrize("nm", TUBES)
def test_tube_two_routes_agree(S, nm):
    # Cutting lines vs explicit 2N-atom cell: same eigenvalues (exact equivalence) -> 1e-8 eV;
    # gaps from separate refinements -> 1e-5 eV.
    r = _row(S["tubes"], *nm)
    assert r["realspace_vs_cutting_max_diff_eV"] < 1e-8
    assert r["gap_realspace_eV"] == pytest.approx(r["gap_raw_eV"], abs=1e-5)


@pytest.mark.parametrize("nm", TUBES)
def test_tube_gap_times_d(S, nm):
    # E_g d ~ 2 a_cc t = 0.767 eV nm; zone-folding warping correction O(a_cc/d) <~ 5% for d > 0.5
    # nm (e.g. exact (7,0): -4.7%). Allow 10%.
    n, m = nm
    r = _row(S["tubes"], n, m)
    if (n - m) % 3 == 0:
        pytest.skip("metallic")
    assert r["gap_times_d_eV_nm"] == pytest.approx(REF_EGD, rel=0.10)
    assert r["ref_2acc_t_over_d_eV"] * r["d_nm"] == pytest.approx(REF_EGD, rel=1e-9)


def test_tubes_csv_matches_summary(S, tubes_csv):
    for nm in TUBES:
        c, s = _row(tubes_csv, *nm), _row(S["tubes"], *nm)
        assert float(c["gap_eV"]) == pytest.approx(s["gap_eV"], abs=1e-12)


# ---------------------------------------------------------------- sweep
def test_sweep_population(S, sweep_csv):
    expect = set(tubes_in_window(0.5, 2.5))
    got = {(int(r["n"]), int(r["m"])) for r in sweep_csv}
    assert got == expect, "sweep must contain every (n,m), 0<=m<=n, 0.5<=d<=2.5 nm"
    assert S["gap_sweep"]["n_tubes"] == len(expect)


def test_sweep_metallicity_rule(S, sweep_csv):
    for r in sweep_csv:
        n, m, g = int(r["n"]), int(r["m"]), float(r["gap_eV"])
        metal = (n - m) % 3 == 0
        assert (abs(g) < 1e-5) is metal, (n, m, g)
        assert (r["K_on_cutting_line"] == "True") is metal, (n, m)
    gs = S["gap_sweep"]
    assert gs["metallic_iff_n_minus_m_mod3_eq_0"] and gs["metallic_iff_K_on_line"]


def test_sweep_zigzag_exact(sweep_csv):
    # every zigzag tube in the sweep against the exact zone-folding formula (1e-5 eV)
    for r in sweep_csv:
        n, m = int(r["n"]), int(r["m"])
        if m == 0:
            exp = 0.0 if n % 3 == 0 else zigzag_gap(n)
            assert float(r["gap_eV"]) == pytest.approx(exp, abs=1e-5), (n, m)


def test_sweep_chiral_sample_independent(sweep_csv):
    # 12 chiral semiconducting tubes spread over the window vs independent line scan (1e-5 eV)
    rows = [r for r in sweep_csv if int(r["m"]) not in (0, int(r["n"])) and (int(r["n"]) - int(r["m"])) % 3]
    rows.sort(key=lambda r: float(r["d_nm"]))
    for r in rows[:: max(1, len(rows) // 12)]:
        n, m = int(r["n"]), int(r["m"])
        assert float(r["gap_eV"]) == pytest.approx(2 * T * tube_min_f(n, m), abs=1e-5), (n, m)


def test_sweep_gap_vs_diameter(S, sweep_csv):
    # E_g d vs 2 a_cc t: all semiconducting tubes within 10%; d > 1.5 nm within 5% (O(a_cc/d));
    # mean within 2%.
    gd, d = [], []
    for r in sweep_csv:
        if (int(r["n"]) - int(r["m"])) % 3:
            gd.append(float(r["gap_eV"]) * float(r["d_nm"])); d.append(float(r["d_nm"]))
    gd, d = np.array(gd), np.array(d)
    assert np.all(np.abs(gd / REF_EGD - 1) < 0.10)
    assert np.all(np.abs(gd[d > 1.5] / REF_EGD - 1) < 0.05)
    assert gd.mean() == pytest.approx(REF_EGD, rel=0.02)
    assert S["gap_sweep"]["gap_times_d_eV_nm_mean"] == pytest.approx(gd.mean(), rel=1e-9)
    assert S["gap_sweep"]["n_semiconducting"] == len(gd)


def test_sweep_consistent_with_requested_tubes(S, sweep_csv):
    for n, m in TUBES:
        r = _row(sweep_csv, n, m)
        assert float(r["gap_eV"]) == pytest.approx(_row(S["tubes"], n, m)["gap_raw_eV"], abs=1e-5)


# ---------------------------------------------------------------- t' and Delta (graphene)
def _variants(S):
    return list(S["tp_delta"].values())


def test_tp_delta_variants_present(S):
    v = _variants(S)
    combos = {(r["tp"], r["delta"]) for r in v}
    for need in [(0.0, 0.0), (0.27, 0.0), (0.54, 0.0), (0.0, 0.05), (0.0, 0.1), (0.0, 0.25), (0.27, 0.1)]:
        assert any(abs(a - need[0]) < 1e-9 and abs(b - need[1]) < 1e-9 for a, b in combos), need


def test_tp_delta_dirac_point_energy_and_gap(S):
    # Exact: min direct gap = 2 Delta at K/K'; Dirac/mid-gap energy = 3 t' (tolerance 1e-6 eV).
    for r in _variants(S):
        tp, de = r["tp"], r["delta"]
        assert r["min_direct_gap_eV"] == pytest.approx(2 * de, abs=1e-6)
        assert np.mean(r["E_at_min"]) == pytest.approx(3 * tp, abs=1e-6)
        k = np.array(r["k_at_min"])
        assert equiv_mod_G(k, K_PT, 1e-6) or equiv_mod_G(k, KP_PT, 1e-6)
        np.testing.assert_allclose(r["E_Gamma"], bands(GAMMA, tp, de), atol=1e-9)
        np.testing.assert_allclose(r["E_M"], bands(M_PT, tp, de), atol=1e-9)


def test_tp_delta_band_edges(S):
    # Exact band extrema as functions of |f| in [0,3]; mesh (shifted MP, nk>=400) misses Gamma by
    # ~|b|/800 -> error ~1e-3 eV; allow 1e-2 eV and require reported <= exact (mesh is a subset).
    x = np.linspace(0, 3, 300001)
    for r in _variants(S):
        tp, de = r["tp"], r["delta"]
        e0 = -tp * (x ** 2 - 3)
        root = np.sqrt(de ** 2 + (T * x) ** 2)
        lo, hi = (e0 - root).min(), (e0 + root).max()
        assert r["band_bottom_eV"] == pytest.approx(lo, abs=1e-2) and r["band_bottom_eV"] >= lo - 1e-9
        assert r["band_top_eV"] == pytest.approx(hi, abs=1e-2) and r["band_top_eV"] <= hi + 1e-9


def test_tp_keeps_vF_and_breaks_eh_symmetry(S):
    # (E_c - E_v)/2 = t|f| does not depend on t' -> v_F identical to NN (1e-4 rel).
    # e-h asymmetry at q = 0.2 1/A toward Gamma = (E_c+E_v) - 2E_D = -2 t' |f|^2 exactly.
    u = (GAMMA - K_PT) / np.linalg.norm(GAMMA - K_PT)
    f2 = f_abs(K_PT + 0.2 * u) ** 2
    for r in _variants(S):
        if r["delta"] != 0:
            continue
        assert r["v_F_m_s"] == pytest.approx(VF, rel=1e-4)
        assert r["E_D_eV"] == pytest.approx(3 * r["tp"], abs=1e-6)
        assert r["eh_asym_at_q0.2_toward_G_eV"] == pytest.approx(-2 * r["tp"] * f2, abs=1e-9)


def test_delta_effective_mass(S):
    # Massive Dirac: E_c = sqrt(Delta^2 + (hbar v q)^2) - t'(1.5 a_cc q)^2 + ... at K, so the band-edge
    # curvature gives hbar^2/m* = (hbar v)^2/Delta - 2 t' (1.5 a_cc)^2. A band-edge mass is a q->0
    # property; require 1% agreement (a fit window reaching hbar v q ~ Delta biases it).
    for r in _variants(S):
        de, tp = r["delta"], r["tp"]
        if de == 0:
            continue
        c2 = HV ** 2 / (2 * de) - tp * (1.5 * ACC) ** 2          # eV A^2
        mstar = (6.582119569e-16 ** 2 * E_CHARGE / (2 * c2 * 1e-20)) / M_E
        assert r["m_star_over_me"] == pytest.approx(mstar, rel=0.01), (de, tp)


def test_dos_gap_with_delta(S):
    # DOS inside the mass gap is exactly zero; at a gapless Dirac point it is ~0 (one 15 meV bin:
    # < slope * 1 bin = 0.05*0.015 ~ 1e-3, allow 2e-3).
    for r in _variants(S):
        if r["delta"] > 0:
            assert r["dos_at_E_D"] == 0.0
        else:
            assert r["dos_at_E_D"] < 2e-3


# ---------------------------------------------------------------- t' and Delta (tubes)
def test_tubes_with_tp_delta(S):
    # E_c, E_v are monotone in |f| (|f| small at the edges, t' < t/6), so the band edges sit at the
    # min-|f| point of the closest line: gap = sqrt(E_g0^2 + 4 Delta^2), independent of t'.
    # Metallic tubes: t' keeps them metallic; Delta opens exactly 2 Delta. Tolerance 1e-5 eV.
    g0 = {(r["n"], r["m"]): (0.0 if r["metallic"] else r["gap_raw_eV"]) for r in S["tubes"]}
    rows = S["tubes_tp_delta"]
    assert len(rows) == 15
    for r in rows:
        model = r["model"]
        de = 0.10 if "Delta" in model else 0.0
        expect = np.sqrt(g0[(r["n"], r["m"])] ** 2 + 4 * de ** 2)
        assert r["gap_eV"] == pytest.approx(expect, abs=1e-5), r
        assert r["gap_realspace_eV"] == pytest.approx(r["gap_eV"], abs=1e-5), r


# ---------------------------------------------------------------- convergence vs pilot
def test_converged_vs_pilot(S, S_pilot):
    # Tube gaps are refined (grid-independent): pilot and full must agree to 1e-6 eV.
    for a, b in zip(S["tubes"], S_pilot["tubes"]):
        assert a["gap_raw_eV"] == pytest.approx(b["gap_raw_eV"], abs=1e-6) or (a["metallic"] and b["metallic"])
    # DOS must improve with the mesh.
    assert S["dos"]["max_rel_dev_from_analytic_away_from_vHs"] < S_pilot["dos"]["max_rel_dev_from_analytic_away_from_vHs"]
    assert S["fermi_velocity"]["per_dirac_point"]["dirac_0"]["v_mean_m_s"] == pytest.approx(
        S_pilot["fermi_velocity"]["per_dirac_point"]["dirac_0"]["v_mean_m_s"], rel=1e-6)
