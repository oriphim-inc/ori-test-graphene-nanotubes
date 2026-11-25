> **Graphene and carbon nanotube bands.** A project Ori built and ran on a problem from Oriphim's held-out test set, which Ori was never tuned on. Short paper: [paper/paper.pdf](paper/paper.pdf). The Sceptic's review is in [notes/review.md](notes/review.md); limits are in [LIMITS.md](LIMITS.md). [oriphim.com](https://oriphim.com)

# Tight-binding graphene and zone-folded carbon nanotubes

By **Oriphim**, with Ori. Asked: build a nearest-neighbour (NN) tight-binding model of graphene
(t = 2.7 eV) and compute its bands, Dirac points, Fermi velocity and DOS. Roll it into the tubes
(5,5), (9,0), (10,0), (8,4) and (13,0), find which are metallic, and plot gap vs diameter. Then add
NNN hopping t′ and a sublattice mass Δ. The agreed plan is in [PLAN.md](PLAN.md).

## Answer

- **Dirac points**: found numerically (not assumed) at exactly the two inequivalent BZ corners,
  K = (2π/3a_cc, ±2π/(3√3 a_cc)) = (1.4749, ±0.8515) Å⁻¹, where E_c − E_v < 1e-14 eV.
- **Fermi velocity**: ħv_F = 5.751 eV·Å, so **v_F = 8.737 × 10⁵ m/s**. It is isotropic to 2e-7 and
  matches 3ta_cc/2ħ. This is the model value for t = 2.7 eV. Measured graphene values are
  ≈ 1 × 10⁶ m/s (see LIMITS.md).
- **DOS**: vanishes linearly at E = 0 (slope 0.051 states/eV²/cell/spin) and has logarithmic
  van Hove peaks at ±t = ±2.70 eV (the M point). Bands end at ±3t = ±8.1 eV. The numerical DOS
  on a 3000×3000 k-mesh matches the analytic Hobson–Nierenberg DOS to within 2% per bin, away from
  the singularities.
- **Nanotubes** (zone folding, NN, t = 2.7 eV):

  | (n,m) | d (nm) | θ (°) | type | gap (eV) | E_g·d (eV·nm) |
  |---|---|---|---|---|---|
  | (5,5)  | 0.678 | 30.0 | metallic (armchair) | 0 | – |
  | (9,0)  | 0.705 | 0    | metallic (zone folding) | 0 | – |
  | (10,0) | 0.783 | 0    | semiconducting | 0.948 | 0.742 |
  | (8,4)  | 0.829 | 19.1 | semiconducting | 0.908 | 0.752 |
  | (13,0) | 1.018 | 0    | semiconducting | 0.735 | 0.748 |

  Metallic ⇔ a cutting line passes through K ⇔ (n−m) mod 3 = 0. This holds for all 314 tubes
  with 0.5–2.5 nm diameter. The semiconducting gaps follow **E_g ≈ 2a_cc t / d = 0.767 eV·nm / d**:
  the mean E_g·d is 0.766 eV·nm and the range is 0.731–0.795. The spread comes from trigonal
  warping and splits the two (n−m) mod 3 families at small d.
- **t′ (NNN)**: the cones survive. The Dirac point moves rigidly to E = 3t′ (0.81 eV for
  t′ = 0.1t = 0.27 eV), and v_F is unchanged to first order. Electron–hole symmetry breaks
  (−0.11 eV asymmetry at q = 0.2 Å⁻¹ for t′ = 0.27 eV), and the band runs from −3t−6t′ up to 3t−6t′.
  Tube gaps are unchanged by t′.
- **Δ (sublattice mass)**: opens a direct gap of exactly **2Δ** at K and K′. The cones become
  massive Dirac hyperbolas with m* = Δ/v_F² (0.023 mₑ for Δ = 0.1 eV), and the DOS is zero inside
  the gap. In the tubes the metallic (5,5) and (9,0) open a 2Δ gap. Semiconducting gaps become
  √(E_g0² + 4Δ²), e.g. (10,0): 0.948 → 0.969 eV for Δ = 0.1 eV.

Note that the t′ and Δ effects on the tubes follow exactly from the model, because both depend
only on |f(k)| at fixed k. They are consistency results, not discoveries.

## Evidence

- `results/fig_band_surface.png`: both bands over the extended zone, plus a contour of E_c
  with the BZ hexagon and the numerically found Dirac points.
- `results/fig_bands_GKMG.png`: Γ–K–M–Γ (Γ: ±8.1 eV, K: 0, M: ±2.7 eV).
- `results/fig_dirac_cone.png`: the cone at K, with E_c along K→Γ and K→M against ħv_F q
  (trigonal warping shows above ~0.1 Å⁻¹).
- `results/fig_dos.png`, `results/dos.csv`: numerical vs analytic DOS, and a low-energy zoom.
- `results/fig_tube_bands.png`: 1D bands of the five tubes. Black lines are cutting lines and red
  dots are an independent explicit 2N-atom real-space unit cell; the two agree to 1e-14 eV.
- `results/fig_gap_vs_diameter.png`, `results/gap_sweep.csv`, `results/tubes.csv`: gap vs
  diameter for every tube from 0.5 to 2.5 nm, with the requested tubes starred and the 2a_cc t/d curve.
- `results/fig_tp_delta.png`: bands, near-K zoom and DOS with t′ and Δ. Display caveat: the DOS
  window is ±10.5 eV, so the bottom of the t′ = 0.54 eV band (−11.34 eV) is cut off. Some curves
  also coincide and hide each other (NN-only lies under Δ = 0.10 away from K; t′ = 0.27 lies under
  t′ = 0.27 with Δ = 0.10). The numbers in summary.json are not affected.
- `results/summary.json`: every number above, as written by the code.
- Checks: 18 unit tests (`tests/test_model.py`) and 44 independent acceptance checks written by
  the Sceptic before the main run (`tests/checks/`). Its review is in `notes/review.md`. All pass
  (60 passed, 2 skipped by design).

## What this does not show

See **[LIMITS.md](LIMITS.md)**. The main points:
- Zone folding of flat graphene ignores curvature and σ–π hybridisation. Real (9,0) tubes have a
  small gap (~80 meV), and the gaps of sub-nm tubes are reduced in reality.
- The model is orthogonal NN tight binding with no overlap integral.
- The gaps are single-particle values, with no excitonic or many-body corrections.
- The t′ and Δ values are illustrative assumptions.

## Re-run

```bash
uv venv .venv && uv pip install --python .venv/bin/python -e '.[test]'   # or: numpy scipy matplotlib pytest
.venv/bin/python run.py                 # full run, ~20 s; writes results/
.venv/bin/python run.py --size 0.125 --out results_pilot   # quick pilot
.venv/bin/python -m pytest tests -q     # unit tests + Sceptic acceptance checks
```
