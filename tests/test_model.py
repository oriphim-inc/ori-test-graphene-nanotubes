"""Unit tests (Ori's own). The Sceptic's acceptance checks live in tests/checks/."""
import numpy as np
import pytest

from tbgraphene import (Graphene, Tube, dos_analytic_nn, dos_numeric, fermi_velocity,
                        find_dirac_points, min_direct_gap, HBAR_EV_S)


def test_reciprocal_lattice():
    gr = Graphene()
    A = np.array([gr.a1, gr.a2])
    B = np.array([gr.b1, gr.b2])
    assert np.allclose(A @ B.T, 2 * np.pi * np.eye(2))


def test_closed_form_equals_diagonalisation():
    rng = np.random.default_rng(0)
    for tp, de in [(0, 0), (0.27, 0), (0, 0.1), (0.54, 0.25)]:
        gr = Graphene(tp=tp, delta=de)
        k = rng.uniform(-3, 3, size=(500, 2))
        assert np.allclose(gr.bands(k), gr.bands_numeric(k), atol=1e-12)


def test_high_symmetry_energies():
    gr = Graphene(tp=0.27)
    hs = gr.high_symmetry()
    assert np.allclose(gr.bands(hs["G"]), [-3 * 2.7 - 6 * 0.27, 3 * 2.7 - 6 * 0.27])
    assert np.allclose(gr.bands(hs["M"]), [2 * 0.27 - 2.7, 2 * 0.27 + 2.7])
    assert np.allclose(gr.bands(hs["K"]), [3 * 0.27, 3 * 0.27])


def test_dirac_points_found():
    gr = Graphene()
    ks = find_dirac_points(gr)
    assert len(ks) == 2
    hs = gr.high_symmetry()
    for k in ks:
        assert min(np.linalg.norm(k - hs["K"]), np.linalg.norm(k - hs["Kp"])) < 1e-6


def test_fermi_velocity():
    gr = Graphene()
    v, vs, _ = fermi_velocity(gr, gr.high_symmetry()["K"])
    assert v == pytest.approx(3 * 2.7 * 1.42e-10 / (2 * HBAR_EV_S), rel=1e-4)
    assert vs.std() / v < 1e-4


def test_mass_gap():
    gr = Graphene(delta=0.1)
    g, _ = min_direct_gap(gr)
    assert g == pytest.approx(0.2, abs=1e-8)


def test_dos_normalisation_and_shape():
    gr = Graphene()
    edges = np.linspace(-8.2, 8.2, 165)
    d = dos_numeric(gr, 300, edges)
    assert (d * np.diff(edges)).sum() == pytest.approx(2.0)
    E = np.linspace(-8.2, 8.2, 40001)
    assert np.trapezoid(dos_analytic_nn(E, 2.7), E) == pytest.approx(2.0, abs=2e-3)


@pytest.mark.parametrize("nm,N", [((5, 5), 10), ((9, 0), 18), ((10, 0), 20), ((8, 4), 56), ((13, 0), 26)])
def test_tube_cell_counts(nm, N):
    tb = Tube(*nm, Graphene())
    assert tb.N == N
    tb.build_cell()
    assert len(tb.sub) == 2 * N and (tb.sub == 0).sum() == N
    assert abs(tb.Ch @ tb.T) < 1e-9


@pytest.mark.parametrize("nm", [(5, 5), (8, 4), (10, 0)])
@pytest.mark.parametrize("tp,de", [(0, 0), (0.27, 0.1)])
def test_realspace_equals_cutting_lines(nm, tp, de):
    tb = Tube(*nm, Graphene(tp=tp, delta=de))
    kap = np.linspace(-np.pi / tb.Tlen, np.pi / tb.Tlen, 9)
    assert np.allclose(tb.bands_realspace(kap), tb.bands_cutting(kap), atol=1e-10)
