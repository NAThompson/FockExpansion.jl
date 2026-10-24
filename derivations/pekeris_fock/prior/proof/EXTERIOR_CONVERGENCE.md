# Two convergent series and one determinant: a convergence theorem for the finite-channel model

Claude, 2026-10-20. Revision 3: incorporates Codex REVIEW_5 and REVIEW_6 (codex_review/).
This covers the frozen finite-channel cusp model only, not full helium.

## Setting

Fix real symmetric m×m matrices U, D, W and fixed nonnegative integer channel labels K_i,
with Λ = diag(K_i(K_i+4)). For complex E write T = W + 2E·I. The channel ODE is

    Φ'' + (5/R + 2U) Φ' + (−Λ/R² + D/R + T) Φ = 0,                      (1)

Use the first-order form X = (Φ, Φ') ∈ ℂ^{2m}, X' = A(R)X, where A = A0 + A1/R + A2/R²,

    A0 = [0 I; −T −2U],   A1 = [0 0; −D −5I],   A2 = [0 0; Λ 0].

Its trace is tr A = −2 tr U − 5m/R. Let C_j = A_j^{[m]} be the m-th additive compounds, so
C(R) = C0 + C1/R + C2/R². The wedge Ξ = X_1 ∧ … ∧ X_m of any m solutions solves Ξ' = CΞ.

**Hypothesis (G), the spectral gap at E.** The 2m eigenvalues of A0(E), counted with
algebraic multiplicity, split into a slow group S and a fast group F, with m eigenvalues
each, such that

    γ(E) := min_{μ∈F} Re μ − max_{λ∈S} Re λ > 0.

(Here γ denotes the spectral gap, not the coupling constant g of the cusp model.) This is
the split `cusptail` makes. In the d = 6 runs at E*, the slow group lies in
[−0.281, 0.575] and the fast group in [3.378, 4.932], so γ = 2.80.

**Normalisation near a base energy.** Fix E0 with γ(E0) > 0. On a small complex disk 𝔻
around E0:
- The Riesz projector P_S(E) of A0 for the group S is holomorphic.
- σ(E) := tr(A0 P_S) = Σ_{λ∈S} λ is holomorphic, even where individual slow roots cross.
- For real E, σ and β are real, because the split is conjugation invariant.

**Exterior subspace.** Fix 0 < δ < γ. 𝒮_E is the set of solutions of (1) with
|X(R)| = O(e^{(max_S Re λ + δ)R}). By the exponential dichotomy for separated spectral
groups of A0 + O(1/R) (Coppel, *Dichotomies in Stability Theory*, LNM 629, 1978; Ch. 4–5),
𝒮_E is m-dimensional and independent of δ. This is the subspace the panels and `cusptail`
intend. We do not identify it with the weighted L² condition of the full problem; that
identification is part of the model's definition.

**Interior frame.** F_E(R) is the m×m matrix of solutions that are regular at R = 0. It has
the logarithmic seed normalisation of `cuspcoeff`, which is holomorphic in E. Here
"regular" permits the forced logarithms, so F_E need not be holomorphic at R = 0.

**Matching function.** For a holomorphic exterior frame G_E (any choice; see (d)), let
M(R) = [F G; F' G']. Then

    𝒟(E) := exp(−∫_{R0}^{R} tr A) · det M(R) = e^{2 tr U (R−R0)} (R/R0)^{5m} det M(R)

is independent of R. Replacing G by another holomorphic frame multiplies 𝒟 by a
holomorphic nonvanishing function of E, so the zero set of 𝒟 is intrinsic. E is an
eigenvalue of the finite model exactly when 𝒟(E) = 0.

## Theorem

Fix E0 with γ(E0) > 0 and fix ω > (log 2)/γ(E0). Choose a complex disk 𝔻 around E0 on which
γ(E) ≥ γ0 > (log 2)/ω.

**(a) Simplicity.** For E ∈ 𝔻, σ is an algebraically simple eigenvalue of C0, and every
other eigenvalue σ' satisfies Re(σ' − σ) ≥ γ0. No diagonalisability is needed, so repeated
or defective slow blocks are allowed, provided the fast group stays separated.

**(b) Formal solution and Borel germ.** Normalise ξ_0(E) by a fixed linear functional,
ℓ(ξ_0) = 1, with ℓ nonzero on the σ-eigenvector at E0. Then there are unique holomorphic
β(E) and coefficients ξ_n(E) such that e^{σR} R^β Σ ξ_n R^{−n} formally solves Ξ' = CΞ.
The Borel germ h(t) = Σ_{n≥0} ξ_{n+1} t^n / n! extends to a single-valued holomorphic
function on the half-plane {Re t > −γ0}, jointly holomorphic in (t, E) ∈ {Re t > −γ0} × 𝔻.

**(c) Factorial series.** With x = R/ω,

    ξ(R) := ξ_0 + ω^{−1} Σ_{k≥0} c_k B(k+1, x),
    c_k = [u^k] h(−log(1−u)/ω) = Σ_{j≤k} |s(k,j)| ω^{−j} ξ_{j+1} / k!,
    B(k+1, x) = k! Γ(x)/Γ(x+k+1).

This converges absolutely and locally uniformly on {Re R > 0} × 𝔻′, for any smaller disk
𝔻′ ⋐ 𝔻. For each ε > 0 its K-term tail is O(K^{−Re R/ω + ε}), locally uniformly. So the
abscissa of convergence is at most 0. (It can be smaller: the correction series can
terminate.)

**(d) Identification.** Let R^β = exp(β Log R), using the principal Log on Re R > 0. Then
Ξ_rec := e^{σR} R^β ξ(R) solves Ξ' = CΞ. It is the unique solution with
Ξ = e^{σR}R^β(ξ_0 + o(1)) as R → +∞. The wedge of any frame of 𝒮_E is a nonzero constant
multiple of Ξ_rec. In particular Ξ_rec is decomposable and gives a holomorphic exterior
section.

**(e) Interior.** The generalized Frobenius (Fock) series of F_E converges for all R > 0 on
the chosen log sheet, together with its R-derivative. The convergence is locally uniform in
R ∈ (0, ∞) and in E.

**(f) Matching.** Pair the complementary minors of (F; F') with Ξ_rec at any R* > 0, using
the Laplace expansion of the determinant (signs as in `plucker_matching.jl`). Multiplying
by the Abel factor above gives a matching function 𝒟 that is holomorphic on 𝔻′ and
independent of R*.

Let E* ∈ 𝔻′ be a simple zero of 𝒟. Truncate the Fock series to N_F terms and the factorial
series to K terms. Then there are N0 and K0 such that for all N_F ≥ N0 and K ≥ K0, the
truncated function has exactly one zero in a fixed small disk around E*, and that zero
converges to E* (by Rouché's theorem). If E* is real and the truncations respect real
conjugation symmetry, which they do for real data and real normalisations, the nearby
zero is real.

**Conclusion.** Within the finite-channel model, near every energy with a spectral gap,
the eigenvalue is defined by two convergent series and one determinant, at any matching
radius R* > 0 and for any ω with a strict margin over (log 2)/γ. This uses the linear
Plücker (compound) representation. It does not prove convergence of the direct Y-chart
factorial series; see Remark 1.

## Proof

**(a)** Triangularise A0. The eigenvalues of A^{[m]} are the sums over m-element index
sets, with multiplicity. Any index set other than the all-slow one replaces r ≥ 1 slow
roots by r fast roots. That raises the real part by at least rγ0. There is only one
all-slow index set, even when slow roots repeat. So σ is algebraically simple.
Equivalently, the top exterior power of the generalized slow space is one-dimensional. ∎

**(b) The recurrence.** Write Ξ = e^{σR}R^βξ. This gives
R²ξ' = (B0R² + B1R + C2)ξ, with B0 = C0 − σI and B1 = C1 − βI. Then

    B0 ξ_0 = 0,
    B0 ξ_n = ((β − n + 1)I − C1) ξ_{n−1} − C2 ξ_{n−2}.

Let w be the complex-linear left eigenvector, a row vector with wB0 = 0, so wξ_0 ≠ 0. Then:
- the n = 1 solvability condition gives β = w C1 ξ_0/(w ξ_0);
- normalise only ℓ(ξ_0) = 1. For n ≥ 1 write ξ_n = a_n ξ_0 + v_n with v_n ∈ ker ℓ. Solve
  for v_n using B0 restricted to ker ℓ (an isomorphism onto range B0). The scalar a_n is
  fixed by the solvability condition at step n + 1, whose coefficient is n(w ξ_0) ≠ 0;
- in general ℓ(ξ_n) ≠ 0 for n ≥ 1, e.g. ℓ(ξ_1) = 15/8 in Codex's m = 1 example
  (REVIEW_6);
- the resulting ξ_n are unique and holomorphic on 𝔻 (shrink 𝔻 so that wξ_0 ≠ 0 and ℓ is
  transverse to range B0).

**(b) The Borel equation.** The formal Borel series h satisfies

    (tI + B0) h'' + (2I + B1) h' + C2 h = 0                               (2)

as an identity of formal power series, by termwise transformation. The low-order
identities

    B0 h(0) + B1 ξ_0 = 0,
    B0 h'(0) + (I + B1) h(0) + C2 ξ_0 = 0

are the n = 1, 2 recurrences. They are retained, and they are needed in (d).

**(b) Convergence at t = 0.** Let P = ξ_0 w/(w ξ_0). Then
(tI + B0)^{−1} = P/t + holomorphic. With z = (h, h'), write t z' = M(t, E) z, where

    M(0) = Q = [0 0; −P C2  −P(2I + B1)].

Since P B1 P = 0 by the definition of β, we get Q² = −2Q, so Q has eigenvalues 0 and −2.
For n ≥ 1 the coefficient equations involve nI − Q. These matrices are invertible, with
‖(nI − Q)^{−1}‖ = O(1/n) uniformly on 𝔻. The Fuchsian majorant argument then gives
convergence of the formal series near t = 0, uniformly in E.

**(b) Continuation.** Equation (2) is nonsingular away from Σ ∪ {0}, where
Σ = {σ − σ'} ⊂ {Re t ≤ −γ0}. Continue the germ along paths in {Re t > −γ0} \ {0}.
- The punctured half-plane is not simply connected, but a loop around 0 inside the
  convergence disk returns the same germ. So the continuation is single-valued on the
  punctured half-plane.
- The singularity at 0 is removable.
- Hence h is holomorphic on {Re t > −γ0}.
- Joint holomorphy in E follows from holomorphic dependence of solutions of linear ODEs on
  parameters. ∎

**Growth lemma.** For θ > 0 and Ω_θ = {|Im t| ≤ θ, Re t ≥ 0}, there are C and c, uniform
on 𝔻′, such that for t ∈ Ω_θ with |t| ≥ T0:

    ‖h(t)‖ + ‖h'(t)‖ ≤ C exp(c √|t|) |t|^c.

*Proof.* Take T0 > max(2 sup‖B0‖, θ, 1) and the principal square root.
1. Shear with w = (h, t^{1/2}h'). This gives ‖w'‖ ≤ (a|t|^{−1/2} + b|t|^{−1})‖w‖, with
   a, b uniform on 𝔻′.
2. Integrate horizontally from the compact base set {|t| = T0} ∩ Ω_θ. There h and h' are
   bounded, uniformly in E, because they are analytic.
3. Grönwall's inequality bounds w, hence h and t^{1/2}h'.
4. Equation (2) then bounds h''. All boundary terms in the integrations by parts below
   vanish. ∎

**(c)** Put t(u) = −log(1 − u)/ω and φ(u) = h(t(u)).

- **The Δ-domain avoids the singularities.** Choose η > 0 with log(2 + η)/ω < γ0, and
  0 < η' < π/2. On the Δ-domain {|u| < 1 + η, |arg(u − 1)| > π/2 − η', u ≠ 1} we have
  |1 − u| < 2 + η and |arg(1 − u)| < π/2 + η'. So Re t > −log(2 + η)/ω > −γ0, and φ is
  holomorphic on the whole Δ-domain.
- **The growth at u = 1 is sub-polynomial.** As u → 1, Re t → +∞ with
  |Im t| < (π/2 + η')/ω. The growth lemma (with θ = (π/2 + η')/ω) gives
  ‖φ(u)‖ ≤ C exp(c √(log(1/|1 − u|)/ω + O(1))) · (log)^c. That is O(|1 − u|^{−ε}) for every
  ε ∈ (0, 1), uniformly on 𝔻′.
- **Coefficient bound.** The O-transfer theorem (Flajolet & Odlyzko, "Singularity analysis
  of generating functions", SIAM J. Discrete Math. 3 (1990) 216–240, Thm 1), applied
  componentwise, gives ‖c_k‖ = O(k^{ε−1}), uniformly on 𝔻′. It uses the whole Δ-domain.
- **Laplace integral.** For Re R > 0, the integral ∫_0^∞ e^{−Rt}h dt converges by the
  growth lemma, and t = t(u) turns it into ω^{−1}∫_0^1 (1 − u)^{x−1} φ(u) du.
- **Termwise integration and absolute convergence.** On compact subsets of Re x > 0,
  |B(k+1, x)| ≤ C k^{−Re x} uniformly. Take ε < inf Re x. Then
  ‖c_k‖ |B(k+1, x)| ≤ C k^{−1−Re x+ε}, which is summable. Termwise integration is justified
  by the positive majorant Σ‖c_k‖u^k = O((1 − u)^{−ε}) against (1 − u)^{Re x − 1}. This
  gives absolute and locally uniform convergence.
- **Tail.** The same bounds give a tail of O(K^{−Re x+ε}).
- **The Stirling formula** follows from (−log(1 − u))^j/j! = Σ_k |s(k,j)| u^k/k!. ∎

**(d)**
- **Ξ_rec solves the ODE.** For Re R > 0, integrate (2) against e^{−Rt}. Boundary terms at
  ∞ vanish by the growth lemma. The terms at t = 0 are matched by the retained low-order
  identities. This gives R²ξ' = (B0R² + B1R + C2)ξ.
- **Uniform asymptotics.** Integration by parts gives
  ξ(R) = ξ_0 + h(0)/R + R^{−1}∫e^{−Rt}h'(t)dt. So ξ = ξ_0 + O_δ(1/|R|) uniformly on
  Re R ≥ δ, not just on the real ray.
- **Uniqueness.** Suppose two solutions are both e^{σR}R^β(ξ_0 + o(1)). Their difference is
  o(e^{(Re σ)R}R^{Re β}). That is below the dichotomy threshold e^{(Re σ + γ0/2)R} of the
  compound system C0 + O(1/R). By (a), solutions below that threshold form the
  one-dimensional σ-line, so the difference is a multiple of Ξ_rec. Its leading
  coefficient is o(1) against ξ_0, so the multiple is 0.
- **The slow wedge** (repairs the earlier submultiplicativity step):
  1. Block-diagonalise A0 + O(1/R) asymptotically along the S/F spectral splitting
     (Coppel Ch. 5; Wasow §§11–12 for separated blocks). In coordinates converging to the
     fixed splitting, 𝒮_E is the solution space of the slow block X_S' = (A_S + o(1))X_S.
  2. By Liouville's formula, det H_S(R) = c·exp(∫ tr(A_S + o(1))) = exp(σR + o(R)), with
     c ≠ 0.
  3. The coordinate changes and their exterior powers are bounded with bounded inverses.
     So the wedge of a frame of 𝒮_E is nonzero and O(e^{(Re σ + ε)R}) for every ε > 0.
  4. Hence it lies in the σ-line and equals c·Ξ_rec.
  This works for defective A_S and nonreal spectra. ∎

**(e)** In Z = (Φ, RΦ') the system is R Z' = (H + RJ + R²L) Z, with

    H = [0 I; Λ −4I],   J = [0 0; −D −2U],   L = [0 0; −T 0].

The exponents are K_i and −K_i − 4. The point R = 0 is regular singular, and the
coefficients are polynomial in R, so there is no other finite singular point.
- **Convergence.** The Frobenius solutions with exponents K_i, including the forced
  logarithms (of bounded degree), converge for all |R| < ∞ on the chosen log sheet. The
  proof is the same Fuchsian majorant argument as in (b): the coefficient equations
  involve (n + K_i)I − H, which is invertible off the resonant set with inverse O(1/n),
  and the bounds are uniform for E in compact sets.
- **Analyticity in E.** The seed normalisation is holomorphic in E, so the series is too.
  Termwise differentiation in R is valid locally uniformly on (0, ∞) × 𝔻′. ∎

**(f)**
- **Holomorphy.** 𝒟 is holomorphic on 𝔻′ by (c)–(e).
- **Uniform convergence of truncations.** The truncation errors tend to 0 uniformly on a
  closed disk around E*, as min(N_F, K) → ∞, by the uniform tail bounds of (c) and (e).
- **One nearby zero.** Choose the disk so that it contains no other zero of 𝒟 and none on
  its boundary. Rouché's theorem gives exactly one zero for all N_F ≥ N0 and K ≥ K0.
- **Independence of R\*.** Abel's identity, with the corrected factor
  exp(−∫tr A) = e^{2 tr U (R−R0)}(R/R0)^{5m}, makes 𝒟 independent of R*. ∎

## Remarks

1. **The Y-chart.** Assume, in addition, transversality: the projection of the slow
   spectral space of A0 to the Φ-coordinates is invertible. (G) alone does not give this.
   Codex's counterexample is m = 2, U = diag(2, −2), T = 3I.
   - **Y as ratios.** Under transversality, ξ_0,top ≠ 0 and Y_0 exists. Then
     Y_ij = ± Ξ_{I(i,j)}/Ξ_top is a ratio of convergent Plücker series.
   - **Finitely many poles in each half-plane.** By the uniform estimate in (d),
     ξ_top ≠ 0 for |R| large in Re R ≥ δ. So Ξ_top has finitely many zeros in Re R ≥ δ.
     Some of them may be removable for individual entries.
   - **A pole bounds the abscissa.** A genuine pole p of Y forces the abscissa of the
     direct Y factorial series to be at least Re p.
   - **Not proved.** Equality with a pole boundary is not proved. Neither is the Borel
     singularity set of the nonlinear Riccati series. Plausibly it consists of the
     fast-minus-slow graph-linearisation differences and their sums, but no derivation is
     given.
2. **Rate.** The proof gives only an upper bound, O(K^{−Re R/ω + ε}). Codex's heuristic
   for the observed slowness is that an R^{−q} origin singularity of Ξ_rec contributes
   (log K)^{q−2}, with q ≈ 31 at d = 2.
3. **Finite checks at E\*.** Two hypotheses are numerical checks: the gap γ(E*) > 0 (with
   margin over (log 2)/ω), and simplicity of the zero (𝒟'(E*) ≠ 0). For the Y-chart
   computation, transversality is also needed. All three are finite checks, and interval
   arithmetic could turn them into certificates.
4. **Cost.** There are C(2m, m) Plücker components. The theorem certifies the two-series
   representation for every d, but the practical computation uses the Y-chart.
5. **Citations still to pin down.** I have not verified exact theorem numbers in Coppel
   or Wasow. The Fuchsian convergence steps are argued directly above, so they do not
   rest on those citations.
