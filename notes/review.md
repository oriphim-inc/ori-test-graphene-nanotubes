
## Sceptic checks (2026-11-14T01:03, after 2 runs)

44 acceptance checks in tests/checks/ (conftest.py with an independent graphene/zone-folding implementation, test_acceptance.py). They skip without results/; against results_pilot (CHECKS_RESULTS_DIR=results_pilot) 39 pass, 3 skip, 2 fail: DOS pointwise deviation (mesh-limited, expected to pass at full size) and the Delta effective mass (a real fit-window bias that full size will not fix).

- check: Run metadata: t=2.7 eV, a_cc=1.42 A; results/ must be size=1.0
- check: Closed form vs 2x2 diagonalisation on G-K-M-G path: max diff < 1e-10 eV
- check: High-symmetry energies: Gamma=-+3t, K=0, M=-+t to 1e-9 eV
- check: Trigonal-warping energies recomputed with an independent f(k) (other gauge) to 1e-9 eV; E toward M < toward Gamma
- check: Exactly 2 Dirac points, |E_c-E_v|<1e-6 eV (recomputed independently), one equivalent to K and one to K' modulo reciprocal lattice (1e-6 reduced coords)
- check: v_F within 1e-4 rel of 3 t a_cc/(2 hbar)=8.737e5 m/s (units sanity 8.6e5..8.8e5), anisotropy < 1e-4, E_D=0 to 1e-6
- check: DOS: integral 2 (1e-6), analytic column integral 2 (1%), zero outside +-3t, e-h symmetric (1e-3 of max)
- check: DOS vs analytic Hobson-Nierenberg in 0.3 eV<|E|<2.9t, >0.1 eV from vHs: pointwise <=15%, 11-bin smoothed <=2%, L1 <=0.05 states/cell (pilot nk=375 gave 42%/3.1%/0.16; ~1/nk scaling)
- check: van Hove peaks at +-t within 2 bins; low-E slope within 4% of A_c/(pi (hbar v_F)^2) using analytic v_F; DOS(0) bin < slope*2*dE
- check: Tube geometry (d, theta, N, 2N atoms, |T|) vs independent Saito formulas to 1e-9
- check: Five tubes: metallic iff (n-m) mod 3 = 0 iff K on cutting line; gaps vs independent line scan of 2t min|f| over k.C_h=2 pi mu, 1e-5 eV; zigzag vs exact 2t min_q|1+2cos(q pi/n)|; metallic gaps <1e-5 eV both routes
- check: Cutting-line vs real-space-cell eigenvalues < 1e-8 eV; gaps agree to 1e-5 eV
- check: E_g*d within 10% of 2 a_cc t = 0.7668 eV nm for semiconducting requested tubes
- check: Sweep: exactly the set of (n,m), 0<=m<=n, 0.5<=d<=2.5 nm; metallicity rule for every tube; every zigzag exact (1e-5 eV); ~12 chiral semiconducting tubes vs independent line scan (1e-5 eV); E_g d within 10% all, 5% for d>1.5 nm, mean within 2%; summary consistent with CSV
- check: Graphene t'/Delta: min direct gap = 2 Delta and mid-gap energy = 3t' to 1e-6 eV, at K/K'; E(Gamma), E(M) exact; band bottom/top vs exact extrema over |f| in [0,3] within 1e-2 eV (and on the correct side)
- check: t': v_F unchanged (1e-4 rel), E_D=3t', e-h asymmetry = -2 t'|f(q)|^2 exactly (1e-9)
- check: Delta: band-edge effective mass m* within 1% of the q->0 curvature hbar^2/m* = (hbar v)^2/Delta - 2t'(1.5 a_cc)^2
- check: DOS at E_D: exactly 0 in a mass gap; < 2e-3 at gapless Dirac points
- check: Tubes with t'/Delta: gap = sqrt(E_g0^2 + 4 Delta^2), independent of t' (exact, 1e-5 eV); real-space route agrees to 1e-5
- check: Convergence vs pilot: tube gaps equal to 1e-6 eV, DOS deviation decreases, v_F equal to 1e-6 rel
- suspect: Effective mass under Delta (run.py part_tp_delta) is fitted over q up to 4e-3 1/A, where hbar v q ~ 0.023 eV is not small vs Delta; pilot m*/me = 0.01260 vs q->0 value 0.01152 for Delta=0.05 eV (+9%), +2.4% for 0.10 eV. Size-independent, so it will fail the 1% check in the main run unless the fit window is scaled to q << Delta/(hbar v) (or curvature taken analytically/by finite difference at K).
- note: m_star_check_Delta_over_vF2_over_me for the t'+Delta variant ignores the t' curvature term (-2t'(1.5a_cc)^2, ~0.7%); the check uses the full expression.
- note: Pilot DOS pointwise deviation from analytic is 42% at nk=375 (low-E bins with few states); acceptable only if the full mesh brings it under 15% (expected ~5%).
- note: With t' the zone-folded tube gaps are exactly unchanged and with Delta they are exactly sqrt(Eg0^2+4Delta^2), because both bands are monotone in |f|; this is expected, not a finding to present as new.

## Sceptic review (after full run, size 1.0)

Tests: 60 passed, 2 skipped (metallic tubes have no E_g*d), 0 failed. The earlier suspect is fixed (m* fit window q<=0.02 Delta/hbar v): m* agrees with the analytic q->0 curvature to 2e-4 including the t' term. No hard-coded results found; the metallic gap is set to 0 only via METAL_TOL=1e-5 eV, and gap_raw_eV is reported alongside. Units were checked: v_F 8.737e5 m/s, m* in m_e, DOS per spin per cell. Figures checked: BZ hexagon and Dirac points correct, (10,0) flat bands at +-t, gap-vs-d is consistent with 2 a_cc t/d.
No suspect findings. Notes: fig_tp_delta DOS truncates the t'=0.54 band bottom (-11.34 eV, range +-10.5 eV), and coincident curves hide each other; README.md is still a stub. Limits: see LIMITS.md.

## Sceptic review (2026-11-14T01:08, after 4 runs)

Full run passes all acceptance checks (60 passed, 2 skipped by design for metallic tubes). The earlier m* suspect is fixed: m* now agrees with the analytic q->0 curvature, including the t' term, to 2e-4. Code, units and figures were checked: v_F = 8.737e5 m/s matches 3ta_cc/2ħ; DOS integrates to 2 per spin per cell with vHs at ±t; the BZ hexagon has the Dirac points at its corners; tube gaps agree between cutting lines and the real-space cell to 1e-14 eV; E_g·d lies 0.731–0.795 eV nm against 0.7668; t'/Δ effects are exact (E_D = 3t', gap 2Δ, tube gap sqrt(Eg0²+4Δ²)). No hard-coded answers. No suspect findings. LIMITS.md written.

- limit: Zone folding of flat graphene ignores curvature and σ–π hybridisation. (9,0) is reported exactly metallic, while experiment shows a curvature gap of about 80 meV (Ouyang 2001). Gaps of the small-diameter tubes (10,0) 0.948 eV and (8,4) 0.908 eV are zone-folded TB values and would be reduced by hybridisation. These are single-particle gaps: no excitons and no GW.
- limit: Orthogonal NN π TB without overlap s, so the bands are e-h symmetric by construction. v_F = 8.74e5 m/s follows directly from the input t = 2.7 eV and sits below the measured ~1e6 m/s.
- limit: t' (0.27 and 0.54 eV) and Δ (0.05–0.25 eV) are assumed illustrative values. All Δ values exceed the literature hBN masses of ~15–30 meV. The sign of t' is convention-dependent: H contains -t' g, as in Castro Neto, so E_D = +3t', and only the e-h asymmetry is physical. At t' = 0.2t the conduction-band maximum moves off Γ, which is a model artefact. Δ on nanotubes is a formal exercise. The t'-invariance of tube gaps is exact within zone folding only.
- limit: The DOS is per spin per cell; the analytic Castro Neto Eq. 14 curve was halved to match. Van Hove peak heights depend on the bin width; only their positions are meaningful. The low-E slope is 1.2% above A_c/(πħ²v_F²) because of trigonal warping in the 0.05–0.5 eV fit window, not because v_F is wrong.
- limit: Metallic classification uses METAL_TOL = 1e-5 eV (raw (5,5) gap 2.7e-8 eV is set to 0; gap_raw_eV is reported). The rule metallic iff (n-m) mod 3 = 0 is verified for all 314 tubes with 0.5 ≤ d ≤ 2.5 nm only.
- note: fig_tp_delta.png DOS panel: the histogram range is ±10.5 eV, but the t' = 0.54 eV band bottom is at -11.34 eV, so the green curve is truncated. Coincident curves (NN-only under Δ = 0.10; t' = 0.27 under t' = 0.27+Δ) hide each other. Cosmetic: widen edges to ±12 eV and use dashed or offset styles. A rerun is optional.
- note: README.md is a 3-line stub. It needs the question, the key numbers (v_F, the tube table, E_g·d, the t'/Δ effects) and a pointer to LIMITS.md.

## Sceptic draft (2026-11-24T03:19, after 4 runs)

The draft is mostly sound. Every result number in paper.md matches results/summary.json through placeholders, the credit line is exact, the PDF is 3 pages, and the Figure 1 caption is accurate. One suspect: the test counts are misreported (44 checks + 18 unit tests = 62, not 60). Following the paper's own reproduce steps also gives 59 passed and 3 skipped, because results_pilot/ is gitignored, absent, and not made by `python run.py`, so the pilot-convergence claim cannot be re-run as written. The rest are limits from LIMITS.md that the paper leaves out, plus wording notes (m* formula and its regime, 'rigid' t′ shift, the reason for t′-invariance, Figure 2 caption, false precision).

- suspect: Sec. 4 test counts are wrong: '60 tests pass (the 44 checks plus 18 unit tests), 0 fail, and 2 are skipped'. 44+18 = 62. The 2 skips are 2 of the 44 checks (test_acceptance.py:187, metallic tubes). Replace with: 'On the full run, 62 tests were collected: 60 pass (all 18 unit tests and 42 of the 44 checks), 0 fail, and 2 checks are skipped by design (gap × diameter is undefined for metals).'
- suspect: The reproduce steps do not reproduce the stated test result or the convergence claim. results_pilot/ is gitignored, is absent from the tree, and is not written by `python run.py`, so test_converged_vs_pilot skips. A fresh run of the Sec. 6 commands gives '59 passed, 3 skipped' (Sceptic re-ran it in a temporary venv). Fix: add `.venv/bin/python run.py --size 0.125 --out results_pilot   # pilot, needed for the convergence check` to Sec. 6 before pytest, and in Sec. 4 after 'Tube gaps converge to 10⁻⁶ eV between the pilot and the full mesh' add '(this check needs the pilot run; without it pytest reports 59 passed, 3 skipped)'.
- limit: Sec. 5 says 'These limits are taken from LIMITS.md' but leaves out several of them. Add: '- **Effective mass and large t′.** m* describes the band only for ħv_F q ≪ Δ; beyond that the dispersion is hyperbolic. At t′ = 0.2t the conduction-band maximum moves off Γ to a ring, which causes the DOS spike near +5 eV in Figure 2; this is a model artefact. Third-neighbour hopping is omitted.' and '- **Idealised tubes and limited convergence.** Tubes are infinite and defect-free, with no environment, doping or temperature. Only the k-mesh was varied; the DOS bin width and minimiser tolerances were not.'
- note: Sec. 3: 'band-edge mass m* = Δ/v_F²' holds only for t′ = 0. Replace with 'band-edge mass given by ħ²/m* = (ħv_F)²/Δ − 2t′(1.5a_cc)² (m* = Δ/v_F² for t′ = 0), valid for ħv_F q ≪ Δ'.
- note: Abstract: '$t'$ shifts the cone rigidly' contradicts the e–h asymmetry stated in the same sentence. Replace with '$t'$ moves the Dirac point to $3t'$ without changing $v_F$, breaks electron–hole symmetry and opens no gap'. In Sec. 3, 'The Dirac point moves rigidly to' should read 'The Dirac point moves to'.
- note: Sec. 3: 'because both bands depend only on |f|' is not enough for t′-invariance of the gaps. Use the LIMITS.md wording: 'because both bands are monotone in |f| (exact within zone folding only; it would not survive curvature)'.
- note: Table 1 caption: 'The cutting-line and real-space routes agree to about 10⁻¹⁴ eV' is true of the eigenvalues, not the gaps. The metallic (5,5) gaps are 2.7e-8 and 1.6e-8 eV, set by the k-refinement. Replace with 'eigenvalues from the cutting-line and real-space routes agree to about 10⁻¹⁴ eV; gaps agree to better than 10⁻⁷ eV'.
- note: Figure 2 caption omits that the DOS panel is clipped at 0.6 states/eV/cell/spin and that the t′ = 0.54 eV curve's spike near +4–5 eV runs off-scale. Append: 'The DOS axis is clipped at 0.6; the t′ = 0.54 eV spike near +5 eV is the off-Γ conduction-band maximum, a model artefact at large t′.'
- note: Sec. 3: 'van Hove peaks at 2.70 eV $= \pm t$' quotes only one peak. Use both placeholders: 'at {{…vhs_peak_energies_eV.0|.2f}} and {{…vhs_peak_energies_eV.1|.2f}} eV (= ∓t)'.
- note: False precision: '873,731 m/s' (6 significant figures) comes from an input t = 2.7 eV given to 2 significant figures. Prefer 'v_F = 8.737×10⁵ m/s, equal to 3ta_cc/2ħ to 6×10⁻⁷' in the abstract and Sec. 3. Hand-typed numbers that duplicate results ('8.74×10⁵' in Sec. 5, '~20 s' in Sec. 6 against 'about 21 s' in Sec. 2) should be placeholders or removed. '9 %' and '2×10⁻⁴' come from notes/review.md; cite it there.
- note: Sec. 4 mentions only the m* suspect from the pilot. The pilot also failed the DOS pointwise check (42 % at nk = 375; mesh-limited, 2.0 % at full size). Add: 'A second pilot failure, the pointwise DOS deviation, was mesh-limited and passes at full size (2.0 %).'

## Sceptic draft (2026-11-24T03:20, after 4 runs)

Re-review of the revised draft: every earlier finding has been fixed correctly. Sec. 4 now gives the right counts: 62 collected, 60 pass (18 unit tests + 42 of 44 checks), 2 skipped by design, and 59 passed / 3 skipped without the pilot. Both match notes/review.md. Sec. 6 now includes the pilot step, and run.py's --size/--out arguments and conftest's results_pilot path confirm it works. Both missing limits have been added to Sec. 5 and match LIMITS.md. The m* formula and its regime, the removal of 'rigid', the monotone-in-|f| reasoning, the Table 1 caption (eigenvalues ~1e-14 eV; gaps <1e-7 eV, with a largest difference of 2.5e-8), the Figure 2 caption, both van Hove peaks, the runtime placeholder (21 s), the removal of 8.74e5, citing notes/review.md for 9 % and 2e-4, and the second pilot failure (DOS) are all present and correct. I checked every placeholder value in paper.md against results/summary.json and they match. No claim goes beyond LIMITS.md. The PDF has 4 pages (/Count 4), within the limit, but I could not render it to check the layout. No suspect findings. Three minor wording notes remain.

- note: Unnecessary precision, carried over from the earlier review: '873,731 m/s' (format ',.0f') appears in the abstract, Sec. 3 and Sec. 5, giving 6 significant figures from t = 2.7 eV. Sec. 3 also reads 'v_F = 873,731 m/s, against the analytic 3ta_cc/2ħ = 873,731 m/s', which looks like a tautology and hides the real 6e-7 agreement. Optional fix in Sec. 3: change to '... = {{...analytic_3t_acc_over_2hbar_m_s|,.0f}} m/s; the two agree to better than 10⁻⁶ (relative).' This is backed by test_acceptance and LIMITS.md ('about 1e-6').
- note: Sec. 3 says the DOS 'matches the analytic DOS to 2.0% per bin away from the singularities'. The metric (run.py:116) covers only 0.3 eV < |E| < 2.9t and more than 0.1 eV from ±t, so it also excludes low energies and the band edges. Replace with: 'matches the analytic DOS to {{results/summary.json:dos.max_rel_dev_from_analytic_away_from_vHs|.1%}} per bin for 0.3 eV < |E| < 2.9t, excluding 0.1 eV around the van Hove peaks.'
- note: Cosmetic: the abstract shows 'isotropic to 2e-07' and 'band splitting is 4e-14 eV' (raw Python exponent format). Consider writing 'isotropic to about 2×10⁻⁷', or keep the placeholder and accept the format.

## Sceptic draft (2026-11-24T03:21, after 4 runs)

Both edits are accurate. The v_F deviation is 6.0e-7 (6.8e-7 at most over directions), below 10^-6. The DOS window wording matches the mask in run.py. paper.md renders the same text. No suspects; the draft is cleared.

- note: v_F wording checked: results/summary.json gives a relative deviation of 6.0e-7 for v_mean and at most 6.8e-7 for v_min/v_max against the analytic 873,730.71 m/s, so 'better than 10^-6 (relative)' is correct.
- note: DOS window wording matches run.py line 116 exactly: |E| > 0.3 eV, |E| < 2.9t, ||E| - t| > 0.1 eV. The 2.0% figure is rendered from dos.max_rel_dev_from_analytic_away_from_vHs.
- note: paper/ is untracked in git, so there is no diff to compare against. Other placeholders that were re-checked (hbar v_F 5.751, v_F 873,731, vHs at -2.70 and 2.70, DOS slope 0.0510) render consistently in paper.md.
