# Phase 2 — Fock angular "atoms" in the Pekeris/Laguerre recurrence (Claude; Codex results credited)

Labels: PROVED, VERIFIED, CITED (repo), CONJ, OPEN.

## A. Where the non-polynomial content is (through order 4)
In Hylleraas distances (r1, r2, r12):
- Polynomial: R ψ10, R²ψ21, R³ψ31, R⁴ψ42, and the Z¹, Z³ and Y₄ₗ parts of R⁴ψ41. VERIFIED: least-squares residual ~1e-16
  (atoms_order2_check.jl); exact homogenizations by Codex.
- R²ψ20 = polynomial + Z·R²χ(α, cos θ). VERIFIED: residual 1.2e-15 after removing Zχ. The phase-1 residual of 6.1e-3 is
  Zχ (1.05e-2 alone), partly offset by ψ21·log(R/s).
- Non-polynomial at orders 3–4: ψ30 (Clausen L/T plus four elliptic moments), ψ41^(2) (the Z² part; Green's-function quadrature) and
  ψ40 (numerical only). CITED from repo src; their closure under differentiation is OPEN.

## Order-2 atoms in closed form
P = |x1+x2| = sqrt(2R² − r12²), d = r1 − r2, s = r1 + r2, L5 = r12·P + s·d.
- α = 2 atan(r2/r1).
- β = asin((R² − r12²)/R²).
- Ω = logs − DL is pure weight 2 with constant coefficients:
  - DL = Cl2(α−β) + Cl2(π+β−α) − Cl2(π−α−β) − Cl2(α+β);
  - logs = β·log((2r1r2)²/L5²) + (β+α)·log((s²−r12²)/(2r1r2)) + (β−α)·log((r12²−d²)/(2r1r2)).
- χ = cΩ/(6π) − (v/3)·log((s+r12)/R) + (c/6)·log((d+r12)²/L5) + r12·Pβ/(3πR²) + algebraic terms, where c = s·d/R² and
  v = (R² − r12²)/R².
- This is one analytic expression, valid for r1 > r2 and r1 < r2 (VERIFIED to 1e-61 against the repo, chi_letters.py).

## B. Closure under differentiation
- PROVED (explicit formula; checked against finite differences to 6e-33): dΩ = α·dlog[(s²−r12²)/(r12²−d²)] + β·(P × rational).
  Every log letter cancels.
- Closed but not minimal set (Claude, 8 atoms): {Ω, α, Pβ, log(R/s), log((s+r12)/R), log((d+r12)²/L5), 1, P}.
- Minimal set (Codex, rank 5, PROVED): {1, L = log((s+r12)/R), α, V = (P/r12)(β+π/2), U}. For the physical log σ ansatz,
  L → L̃ = L − ((π−2)/π)·log(s/R). Adopted as the reference.
- Recurrence: finite (PROVED by construction, every block a polynomial-coefficient operator).
  - Each target row is multiplied by a polynomial clearing factor: uvw for most rows; uvw·(2u²+2uw+2v²+2vw+w²)/2 for the "1" row;
    uvw·P² for the P-type rows.
  - The uvw factor is needed for Laguerre finiteness: Γ has face denominators 1/u, 1/v, 1/w, which otherwise leave bare-derivative
    infinite tails.
  - Mine: 16/15/8 blocks for σ²P/C1/C2 (count corrected by Codex), 25 955 offsets in total, shell half-bandwidth 10.
  - Codex (minimal): 11/10/5 blocks, half-bandwidth 10.
  - Coupling is triangular in transcendental weight (2 → 1 → 0) and in log level (j ← j, j+1, j+2).
- Holonomic picture: for orders 3–4 (elliptic ψ30, Green-integral ψ41^(2) and ψ40), closure is expected from D-finiteness but is not
  derived (CONJ).

## C. Numerics (fixed ε = 1.704; (s, x, y) weak form; tanh–sinh angular quadrature; 320-bit)
- Variants: polynomial; +log; +log²; +order-2 atoms (χ and log(R/s), each times polynomials); +ψ30 (Float64).

| N | polynomial | +log | +log² | +order-2 atoms (size) |
|---|---|---|---|---|
| 8 | 9.9e-7 | 1.6e-9 | 6.5e-10 | 4.6e-14 (289) |
| 9 | 5.0e-7 | 3.0e-10 | 7.5e-11 | 3.9e-15 (404) |
| 10 | 2.7e-7 | 6.2e-11 | 1.2e-11 | 2.3e-16 (545) |
| 11 | 1.5e-7 | 1.4e-11 | 1.6e-12 | 2.4e-17 (713) |

- Order-2 atoms: the log–log slope against size is about −8, against about −3 for polynomial. The semilog curve stays convex, so this is
  faster algebraic convergence, not exponential (NUMERICAL; Codex fits slopes −7.7 to −8.5).
- Quadrature check: rerunning at HINV = 20 and 384 bits leaves the energies and the a21 values (to 12 digits) unchanged.
- ψ30 atom: better at N ≤ 7 (2.7e-14 with 209 functions). From N = 8 the energy falls below E0, so those runs are invalid; the cause is
  Float64 ψ30 values and gradients amplified by cond(S) ~ 1e32.
- a21 (Codex's global finite-part estimator, generalized to atom states; reproduces Codex's earlier values exactly), order-2-atom states,
  N = 8..11: 0.4767906, 0.4767919, 0.4767502, 0.4767484. Reference 0.47674788. Supported: about 6 digits (uncertainty ~1e-6).
  No new digits.
- Plots: convergence_atoms.png (log–log and semilog).
