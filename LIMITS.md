# LIMITS: what this work does NOT show

Written by the Sceptic after the full run (`python run.py`, size 1.0). The model and its
parameters are in PLAN.md ("Model", "Nanotubes", "t′ and Δ"). Every number in results/ is a
property of that **model**. None of them is a prediction for real graphene or real nanotubes
beyond what that model can describe.

## 1. The model itself
- **Orthogonal π-only tight binding.** There is one p_z orbital per atom, nearest-neighbour
  t = 2.7 eV, and no overlap integral s. With s ≈ 0.13 the real π bands are electron–hole
  asymmetric (Saito et al. 1998), and the NN-only bands here are exactly symmetric by
  construction. σ bands, and the π bands far from E = 0, are not described.
- **t = 2.7 eV is an input, not a result.** The Fermi velocity follows from it:
  v_F = 3 t a_cc / 2ħ = 8.74×10⁵ m/s. The numerical value agrees with this formula to about
  1e-6. Measured values (about 1.0–1.1×10⁶ m/s, raised further by electron–electron
  renormalisation at low density) are higher. The v_F reported here is the model's value
  and is not a prediction of the measured one. All gaps scale linearly with t.
- **No spin–orbit coupling, no electron–electron interaction, no excitons, no GW
  corrections, no strain, no substrate, no disorder, no temperature.** Optical gaps of
  real semiconducting nanotubes are dominated by exciton binding and self-energy
  corrections. The "band gap" here is the single-particle TB gap only.
- **DOS convention.** The DOS is per spin per 2-atom unit cell and integrates to 2.
  Multiply by 2 to include spin. Castro Neto et al. 2009, Eq. (14), includes spin. The
  analytic curve in results/dos.csv and fig_dos.png has been divided by 2 to match the
  per-spin convention. The DOS is a histogram on a 3000×3000 mesh with 0.01 eV bins. The van
  Hove peak heights therefore depend on the bin width and are not converged values; the
  log-divergence cannot be resolved this way. Only the peak positions (±t) are meaningful.
- **Trigonal warping is present (it is part of the model), but it is only sampled.** The
  cone is isotropic to 2e-7 within |q| ≤ 2e-3 Å⁻¹. Anisotropy appears only at larger q
  (summary.json "trigonal_warping_Ec_eV"). The low-energy DOS slope fitted over
  0.05–0.5 eV is 1.2 % above A_c/(πħ²v_F²) because of this curvature. This is not an
  error in v_F.

## 2. Nanotubes (zone folding)
- **Pure zone folding of flat graphene.** Curvature, σ–π rehybridisation and the change
  of bond lengths and hoppings around the circumference are not modelled. As a result:
  - The (n−m) mod 3 = 0 non-armchair tubes, including the requested **(9,0)**, come out
    exactly metallic. In reality curvature opens a small secondary gap of order
    10 meV × (1/d²) (Hamada et al. 1992; Kane & Mele 1997; Ouyang et al., Science 292, 702
    (2001) measured roughly 80 meV for (9,0)). Only armchair tubes such as (5,5) are
    protected by symmetry.
  - For small diameters (d below about 1 nm: (5,5), (9,0), (10,0), (8,4)), σ–π
    hybridisation lowers the gaps substantially, and for d ≲ 0.5 nm it can even make
    zone-folded semiconductors metallic (Blase et al., PRL 72, 1878 (1994)). The
    zone-folded gaps of (10,0) at 0.948 eV and (8,4) at 0.908 eV are therefore upper-end
    TB numbers and not predictions. PLAN.md "Scope" lists this as out of scope.
- **Gap vs diameter.** E_g·d deviates from 2 a_cc t by up to −4.7 % (0.731 eV nm) at small
  d, and the (n−m) mod 3 = 1 and 2 families separate. Both effects come from trigonal
  warping in the model. They are not a separate result, and the "family" splitting seen
  experimentally also contains curvature effects that are absent here.
- **Metallicity threshold.** A tube is labelled metallic if its numerical gap is below
  1e-5 eV (run.py METAL_TOL). For example, (5,5) has a raw gap of 2.7e-8 eV, which comes
  from the k-refinement tolerance. Over the 0.5–2.5 nm sweep this classification agrees
  with (n−m) mod 3 = 0 and with "K lies on a cutting line" for all 314 tubes. It is not
  tested outside that diameter range.
- **Infinite, defect-free tubes**, with no ends, bundles, environment or doping.

## 3. t′ and Δ (PLAN.md "t′ and Δ": assumed values)
- **The values are illustrative assumptions, not fits.** t′ = 0.27 and 0.54 eV are 0.1t
  and 0.2t from the Castro Neto range. Δ = 0.05, 0.10 and 0.25 eV are chosen values. The
  hBN-induced masses in the literature are about 15–30 meV, smaller than all of the Δ used.
- **The sign of t′ depends on the convention.** Here H contains −t′ g(k), as in Castro Neto
  2009, Eq. (6), so t′ > 0 puts the Dirac point at +3t′. A uniform shift is not observable.
  The physical content is only the electron–hole asymmetry. Other papers, for example
  Reich et al. PRB 66, 035412 (2002), which fit with overlaps, quote t′ with the opposite
  sign convention and alongside s, so their numbers cannot be transferred directly.
- **NNN hopping only.** Third neighbours, which matter for quantitative fits to DFT, are
  omitted. At t′ = 0.2t the conduction-band maximum moves off Γ to a ring around it, which
  produces the spike near +5 eV in the DOS panel of fig_tp_delta.png. This is a property
  of the model at large t′ and not a robust feature of graphene.
- **Δ applied to nanotubes is a formal exercise.** It shows that a staggered potential
  gaps the metallic tubes, giving gap = sqrt(E_g0² + 4Δ²), which is exact for this model.
  No physical mechanism for a uniform A/B mass on a tube is claimed. Within zone folding,
  t′ leaves every tube gap exactly unchanged. This follows from both bands being monotone
  in |f|, and it would not survive curvature.
- **Effective mass with Δ** (m* = 0.0115, 0.0230 and 0.0576 mₑ) is the curvature at K along
  k_x in the limit q → 0. It agrees with ħ²/m* = (ħv)²/Δ − 2t′(1.5 a_cc)² to 2e-4. It describes
  the band only for ħv q ≪ Δ. Beyond that the dispersion is relativistic (hyperbolic) and
  not parabolic.

## 4. Numerics not covered by a convergence study
- The tube gaps converge to 1e-6 eV between the pilot (size 0.125) and the full run. The
  DOS deviation falls with the mesh size. No other parameter (bin width, Nelder–Mead
  tolerances) was varied.
- In fig_tp_delta.png the DOS histogram range is ±10.5 eV, and the t′ = 0.54 eV band
  bottom lies at −11.34 eV, so that curve is truncated at the left edge of the plot. In
  the same figure, curves that coincide (NN-only and Δ = 0.10 over most of the path, and
  t′ = 0.27 with and without Δ) hide each other.
