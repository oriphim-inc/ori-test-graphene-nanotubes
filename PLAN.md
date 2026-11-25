# PLAN — tight-binding graphene and carbon nanotubes

Requested by Oriphim. Planned by Ori (no attachments; standard textbook physics).

## Question
1. Nearest-neighbour (NN) π-band tight-binding model of graphene, t = 2.7 eV:
   band structure over the full Brillouin zone (BZ) and along Γ–K–M–Γ.
2. Locate the Dirac points numerically, expand around them, extract the Fermi velocity v_F.
3. Density of states (DOS).
4. Zone-folded carbon nanotubes (n,m) = (5,5), (9,0), (10,0), (8,4), (13,0):
   metallic or semiconducting, band gap of the semiconducting ones vs tube diameter.
5. Add next-nearest-neighbour hopping t′ and a sublattice mass Δ: effect on the Dirac cones.

Deliverables: band-surface figure, DOS figure, Fermi velocity, gap-vs-diameter plot
(plus a Γ–K–M–Γ figure, nanotube band figures, t′/Δ figure, results/summary.json).

## Model (agreed physics)
- Orthogonal π-band tight binding (no overlap integral s), one p_z orbital per C atom, spinless
  (DOS reported per spin per unit cell; multiply by 2 for spin).
- H(k) = [[+Δ, −t f(k)], [−t f(k)*, −Δ]] − t′ g(k)·1, with f(k) = Σ_δ exp(i k·δ) over the 3 NN
  vectors and g(k) = Σ over the 6 NNN vectors = |f|² − 3.
  E±(k) = −t′ g(k) ± sqrt(Δ² + t²|f(k)|²). Sign convention follows Castro Neto et al.,
  Rev. Mod. Phys. 81, 109 (2009), Eq. (6) (with Δ added).
- Constants:
  - t = 2.7 eV (given by user; also the Castro Neto 2009 value).
  - a_cc = 1.42 Å C–C distance, lattice constant a = √3 a_cc = 2.46 Å
    (Castro Neto 2009; Saito, Dresselhaus & Dresselhaus, *Physical Properties of Carbon
    Nanotubes*, 1998).
  - ħ = 6.582119569e-16 eV s (CODATA 2018).
- Analytic references used only as checks (never as outputs):
  K = (2π/(3a_cc), 2π/(3√3 a_cc)); v_F = 3 t a_cc / (2ħ) ≈ 8.7e5 m/s;
  NN DOS closed form (Hobson & Nierenberg 1953; Castro Neto 2009 Eq. 14);
  van Hove singularities at |E| = t; band edges at ±3t.

## Nanotubes (zone folding)
- Chiral vector C_h = n a1 + m a2, d = |C_h|/π, translation vector T = t1 a1 + t2 a2
  with t1 = (2m+n)/d_R, t2 = −(2n+m)/d_R, d_R = gcd(2m+n, 2n+m), N = 2(n²+nm+m²)/d_R
  hexagons per cell (Saito et al. 1998, ch. 3–4).
- Cutting lines k = μ K1 + k K2/|K2|, μ = 0…N−1, k ∈ (−π/|T|, π/|T|].
- Metallic ⇔ a cutting line passes through K, determined numerically (expect (n−m) mod 3 = 0).
- Gap = min over all lines and k of (E_c − E_v), found numerically (dense k sampling then
  bounded refinement).
- Independent cross-check: build the explicit real-space unit cell of the tube (2N atoms
  on the unrolled sheet with periodic wrapping along C_h), diagonalize the 2N×2N Bloch
  Hamiltonian H(k), compare bands/gaps to the cutting-line result.
- Context sweep: all (n,m) with 0 ≤ m ≤ n and 0.5 nm ≤ d ≤ 2.5 nm, gap vs d, with the
  reference curve E_g ≈ 2 a_cc t / d (Saito 1998; White & Mintmire 1998) overplotted.
  Requested tubes highlighted.

## t′ and Δ (assumed representative values, labelled assumptions)
- t′: literature range 0.02t ≤ t′ ≤ 0.2t (Castro Neto 2009); use t′ = 0.1t = 0.27 eV as the
  representative value, also show 0.2t. Expected: Dirac point shifts to E = 3t′, cone
  survives, electron–hole symmetry broken, v_F unchanged to first order.
- Δ: hBN-substrate-induced masses are tens of meV (Giovannetti et al., PRB 76, 073103 (2007):
  ~53 meV gap; Hunt et al., Science 340, 1427 (2013): ~30 meV). Use illustrative
  Δ = 0.05, 0.1, 0.25 eV. Expected: gap 2Δ at K, massive Dirac dispersion, DOS gap.
- Also apply t′ and Δ to the five nanotubes (Δ breaks the metallic tubes' degeneracy).

## Numerics (Ori may tune; recorded per run)
- Full-BZ surface grid; DOS from a uniform k-mesh over the rhombic BZ (both bands),
  Gaussian or histogram binning; analytic DOS used as a convergence check.
- v_F from linear fits of E(q) for small |q| along many directions around each Dirac point.
- run.py has a --size option (mesh scaling) for pilot vs full runs.

## What counts as the answer
- Two inequivalent Dirac points found numerically at the K/K′ corners, with |E| < 1e-6 eV.
- v_F (m/s) with its isotropy spread, compared to 3ta_cc/2ħ.
- DOS that matches the analytic NN DOS (vanishing linearly at E = 0, log van Hove peaks at ±t).
- A table: (n,m), d (nm), θ, metallic/semiconducting, gap (eV), from cutting lines and from the
  real-space cell, plus the gap-vs-diameter plot.
- Quantified effect of t′ (Dirac energy shift, asymmetry) and Δ (gap 2Δ).

## Scope / scale
Everything runs on a laptop in seconds to a minute; no scaling down is needed.
Out of scope: curvature/σ–π hybridization, spin–orbit, overlap integral s, GW/DFT corrections
(these are limitations, notably for the small-diameter (5,5), (9,0) and (10,0) tubes).
