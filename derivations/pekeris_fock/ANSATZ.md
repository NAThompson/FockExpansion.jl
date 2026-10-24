# ANSATZ — Pekeris + Fock logarithms (Claude draft, agreed with Codex at sync point 1, 2026-10-21)

Labels: **PROVED** (argument given here), **VERIFIED** (machine/numerical check, file named),
**CITED** (exact reference), **CONJ** (conjectured), **OPEN**.

## 1. What is actually known about the triple point

| Statement | Status |
|---|---|
| FHOS 2005 (CMP 255, 183; math-ph/0312060, Thm 1.1, Eqs. 1.8–1.12): Ψ = e^{F2+F3} φ, φ ∈ C^{1,1}_loc. F2 = Kato cusp terms; F3 = C0 Σ Z (x_i·x_j) ln(\|x_i\|²+\|x_j\|²), C0=(2−π)/(12π) in units H=−Δ+V. In atomic units this is −Z(π−2)/(3π)(r1²+r2²−r12²) ln R + smooth: exactly Fock's ψ21 R² ln R. Only the first logarithm; remainder is C^{1,1}, not analytic. | CITED (Claude and Codex independently; unit conversion x=2r by both) |
| FHOS 2009 (CMP 289, 291; arXiv:0806.1004, Thms 1.1, 1.4): Ψ = Ψ¹ + \|x_k\|Ψ² (or \|x_k−x_l\|), Ψ^i real analytic, **locally near isolated two-particle collisions**; the triple point is excluded. | CITED |
| The hint "Ψ = F1 + R² ln R·F2 with F1, F2 analytic" is not a theorem in either paper, and is **false** for ordinary analyticity: ψ42 ≠ 0 (repo `psi42`, closed form), and analytic F1, F2 restricted to a ray allow at most one power of ln R. | CITED (repo) + PROVED (Codex, codex/codex_literature.md) |
| Fock (1954)/Morgan (1986): formal expansion Σ_n R^n Σ_{p≤⌊n/2⌋} ln^p R ψ_{n,p}; Morgan proves convergence of Fock-type solutions. Codex notes Morgan's abstract leaves open whether the *physical decaying eigenfunction* is representable by it. | CITED (abstract); full text not inspected |
| ln² R occurs (ψ42 ≠ 0, at R⁴), so at least two log sectors are needed. That ln^p R is nonzero for every p (first at R^{2p}) is CONJ: the Fock hierarchy *allows* p ≤ ⌊n/2⌋, but no all-order nonvanishing result was checked. The ansatz allows an infinite tower; necessary infinite log content is not proved. | CITED (repo) for p=2 / CONJ for all p |

## 2. ln R versus ln s, and what happens to the coefficient functions

s = r1+r2, σ = u+v+w = 2εs. On the physical cone R/s ∈ [1/√2, 1] is homogeneous of degree 0 and real-analytic away from R=0.
Then R^n ln^k R ψ_{nk} = Σ_j C(k,j) ln^j σ · R^n ψ_{nk} ln^{k−j}(R/σ). So

Ψ = Σ_j ln^j σ · G_j,  G_j = Σ_{n≥2j} R^n g_{nj}(Ω),  g_{nj} = Σ_{k≥j} C(k,j) ψ_{nk} ln^{k−j}(R/σ)|_{R=1}. **PROVED** (algebra).

The maximal log power and the vanishing orders G_j = O(R^{2j}) are unchanged; the leading coefficient of ln σ is still
ψ21 R² = −Zκ(r1²+r2²−r12²), κ=(π−2)/(3π), because ln(R/σ) only moves angular terms into lower sectors.

**But the G_j are not analytic in (u,v,w) at R=0.** The degree-2 part of G_0, ψ20 + ψ21 ln(R/s), is not a quadratic polynomial in
(r1,r2,r12): least-squares residual 6.1e-3 rms against target rms 2.3 (also not for ln(s+r12): 3.0e-3, or ln R: 1.0e-2).
**VERIFIED** numerically, `julia/check_psi20_polynomial.jl` (uses the repo's ψ20 with its Clausen terms). Frankowski–Pekeris-type
bases add powers of R = ((s²+t²)/2)^{1/2} for this reason (Korobov 2000 Eq. 2, CITED by Codex).

## 3. The chosen global ansatz

Ψ = e^{−σ/2} Σ_{k≥0} ln^k σ · F_k(u,v,w),  F_k = Σ A^{(k)}_{lmn} L_l(u) L_m(v) L_n(w),  ε² = −E (or fixed ε),

with the class restriction (agreed with Codex): each F_k is **log-free**, i.e. has a convergent homogeneous expansion
Σ_n R^n f_{nk}(Ω) at the corner with no hidden ln R, and F_k = O(R^{2k}). Without that restriction the split is not unique
(gauge F0 → F0 − ln σ·G, F1 → F1 + G). **PROVED** (trivial).

Constraints at R=0: F_0(0) = Ψ(0); F_1 = −Zκ Ψ(0)(r1²+r2²−r12²) + O(R³); F_k = O(R^{2k}). In Pekeris variables
r1²+r2²−r12² = (w²+2uw+2vw−4uv)/(8ε²). (Fock/FHOS, CITED.)

Decay: e^{−σ/2} is a *gauge*, not the boundary condition: an unrestricted Laguerre sum can undo it, so decay must be imposed as a coefficient-space admissibility condition (Section 5). Along r1=r2, e^{−εs} = e^{−√2 εR} matches the free hyperradial exponential when ε²=−E; that is motivation, not the interacting angular asymptotics. In the escape channel the true decay is
e^{−√(2(E_He⁺−E)) r2} = e^{−1.344 r2}, slower than the basis' e^{−1.704 r2}, so the series must build the difference; in a 1D model
the Laguerre coefficients of such a function still decay geometrically (ratio |2β−1|/(2β+1) ≈ 0.12, β = 1.344/(2ε)), while a growing
exponential gives geometric *growth*. CONJ for the 3D problem (1D generating-function computation only).

Finite-K truncation: if F_K e^{−σ/2} lies in the physical domain at the simple eigenvalue, σ²P F_K = 0 forces F_K ∝ Ψ, contradicting
F_K = O(R^{2K}); so an exact finite tower has no log content. **PROVED under that domain hypothesis** (Codex's qualification accepted).
Hence any finite-K computation is an approximation; the exact object is the infinite tower.

Why it should be complete: Laguerre completeness gives density of polynomials × e^{−σ/2} in L²; density in the Coulomb *form domain* (a form core) is needed and is NOT proved here (OPEN; Codex point 6c). If it holds, k=0 alone is already complete and the log sectors do not add completeness, they change the *rate*: they remove the R² ln R (and, with k≥2,
R^{2k} ln^k R) obstructions. Remaining rate-limiters: non-polynomial angular parts of G_j (Section 2) and the uncaptured higher logs.

## 4. The recurrence (stage 2; details in COLLAB.md and sym/)

P := −(2v+w)(2u+w)(u+v) e^{σ/2}(H−E)e^{−σ/2} (polynomial operator, cubic coefficients, no u–v mixed derivative).
Sector j (multiplied by σ²):  σ²P F_j + (j+1) C1 F_{j+1} + (j+1)(j+2) C2 F_{j+2} = 0,
C1 = 2σ Σ_j (Σ_i M_ij) ∂_j + σ Σ_i b_i − Σ_ij M_ij,  C2 = Σ_ij M_ij = 4ε²(u+v)(2u+w)(2v+w) ∝ r1 r2 r12 (the physical weight).

| operator on A^{(·)} | offsets (general E) | offsets (ε²=−E) | shell range ΔN |
|---|---|---|---|
| P (k=0, Pekeris) | 57 | **33** | ±3 / ±2 |
| σ²P (diagonal in log level) | 225 | 177 | ±5 / ±4 |
| C1 (couples level j+1 → j) | 117 | 117 | ±4 |
| C2 (couples level j+2 → j) | 57 | 57 | ±3 |

Finite: every derivative carries its own coordinate after the bare-derivative terms cancel (checked symbolically). The coupling is
upper-triangular in the log level (j ← j, j+1, j+2). For a finite K this is solved top-down; the exact tower has no top sector, so
'top-down' is not an algorithm for it — a finite-K closure is an approximation (PROVED qualification, Codex).
VERIFIED: 9/9 spot coefficients equal Koutschan–Zeilberger's machine copy of Pekeris Eq. (22) (theirs = ours/(2ε));
direct-application check (sym/verify_recurrence.py, 3 random seeds, all sectors, exact); Codex's independent derivation agrees
on all 441 coefficients (codex/codex_crosscheck.md).

## 5. Boundary condition as admissibility / minimality (stage 3; owned by Codex, codex/codex_minimality.md)

- PROVED (Codex, finite matrices): at fixed ε, grouping 3 shells (k=0, general E) or 5 shells (log tower) gives block-tridiagonal
  Hermitian H, G with growing rectangular blocks; for Im z ≠ 0 every tail H_L − zG_L is invertible and the backward Schur
  recursion S_i = D_i − U_i S_{i+1}^{-1} V_i is a finite matrix continued fraction.
- PROVED (Codex, conditional): if the nested polynomial spaces form a form core, Galerkin resolvents converge strongly for Im z ≠ 0;
  the eigenvector is recovered as a residue at the real pole, not by real-axis tail inversion.
- The correct "minimality" statement is an admissibility statement: the physical solution is the (one-dimensional, at a simple
  eigenvalue) set of coefficient sequences whose reconstruction lies in the Coulomb form domain and satisfies the weak equation.
  PROVED as a definition + representation theorem; its identification with a Pincherle-type subordinate/minimal subspace of the raw
  recurrence is CONJ. No growing-block Pincherle theorem exists in the sources checked (Gautschi 1981 scalar; Levrie–Bultheel 1996 and
  block-Jacobi deficiency theory: fixed block size).
- Laguerre ℓ² ≠ physical normalizability (PROVED, Codex: weights e^{−σ}du dv dw vs e^{−σ}Q du dv dw, counterexamples both ways);
  sufficient: Σ(1+N)³|A|² < ∞. In 1D, decaying ⇔ geometrically decaying coefficients and growing ⇔ geometric growth (PROVED, generating
  function); "physical ⇔ algebraic decay" is false as a general criterion.
- Fock constants: Ψ(0) = Σ A^{(0)} needs absolute summability (not implied by ℓ²); a21 is best read through the global finite-part
  identity a21 = −(1/8) FP∫ h(R)/R dR − b/8 (a21_global/FINITE_PART.md, re-derived by Codex), which converged ~11× better than the local
  projection on our states (0.47636 vs 0.47245 at N=13; reference 0.4767479).

Rate claim (Codex 6f): log sectors improve the observed convergence (NUMERICAL, julia/convergence.png); no all-N rate theorem.
