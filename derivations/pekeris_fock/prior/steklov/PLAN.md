# Steklov–Poincaré matching for the free Fock constants of helium (plan)

E0 = -2.903724377034119598311 is given. Find the free Fock constants by matching at
the hypersphere R = R* (R^2 = r1^2 + r2^2). The interior is a Fock series. The exterior
is a variational (Ritz) solve in Hylleraas-type functions, not hyperspherical channels.

## Coordinates

Only S states. r1 = R cos(phi), r2 = R sin(phi), r12 = R u, where phi = alpha/2 is in
[0, pi/2] and u is in [|cos phi - sin phi|, cos phi + sin phi]. Then

    dtau = (8/pi) R^5 cos(phi) sin(phi) u dR dphi du      (S-state volume, rescaled)
    dOmega = (8/pi) cos(phi) sin(phi) u dphi du            (normalized, integral 1)

This is Codex's uniform disk measure: q = cos(alpha) = cos(2 phi), s = sin(alpha) cos(theta) = 1 - u^2.
The S-state weight r1 r2 r12 cancels every Coulomb singularity: 1/r1 times cos(phi) is smooth,
and 1/r12 times u is smooth. It also cancels the direction factors in |grad|^2. The cusp
r12 = 0 sits at the endpoint u = |cos phi - sin phi| -> 0, i.e. at phi = pi/4. Splitting the
phi range at pi/4 makes every integrand piecewise analytic, so tensor Gauss–Legendre
in (phi, u) converges spectrally. Exchange symmetry is phi -> pi/2 - phi with u fixed,
so we integrate over phi in [0, pi/4] and double.

## Interior (R < R*): channel-Galerkin Fock series in the cusp gauge (fallback)

Psi = exp(S) Phi, with S = R U(Omega) = -Z(r1+r2) + (g/2) r12. Then

    Phi'' + (5/R + 2U) Phi' - Lambda^2 Phi / R^2 + (2/R) grad_Omega U . grad_Omega Phi + (W + 2E) Phi = 0,
    W = |grad S|^2 = 2Z^2 + g^2/2 - Z g (rhat1 - rhat2) . rhat12.

Galerkin in the exchange-even S harmonics Y_i (labels (m, j), K = 2(m+j), degree l = m+j <= d_int)
gives Codex's channel ODE and Fock recurrence (cuspmodel / cuspcoeff). I rebuild it here,
using the (phi, u) quadrature for U, C (charge), and W, and D_ij = (lam_j - lam_i + 5) U_ij - 2 C_ij.
Regular solutions: Phi = sum_n sum_k R^n (ln R)^k c_{nk} a, with one free seed a_i per channel.
a_i is the coefficient of R^{K_i} Y_i (no log) in Phi, i.e. <Phi_{K_i,0}, Y_i> = a_i.

Why not exact-angle Fock functions: the repo gives psi_{n,k} in closed or quadrature form only
up to n = 4 (psi40 takes about 1 s per point through a Green's function). At R* ~ 1 the
series needs N ~ 50–100 orders for 1e-10, so exact angles are unusable for the interior map.
The repo psi20/psi21/psi31/psi30 are used only as checks: the projections of the channel-model
Phi_{n,k} at low n should match exact-angle projections as d_int grows.

Recurrence memory: I accumulate F(R*), dF/dR(R*) on the fly. I store c_{nk} only for n <= 6,
for the gauge restoration below. That lets d_int go beyond Codex's 24 (m = 169).

## Trace space and Steklov operator

T_L = span{ tau_i = exp(R* U) Y_i : l_i <= L }, with L <= d_int (cusp-gauge trace, as Codex had).
For an interior solution whose trace has coefficients [t; 0] (zero beyond L), a = F(R*)^{-1} [t; 0].

    B_int(s, t) = oint g_s d_R Psi_t dS = R*^5 s' [ G_U t + G_{L,m} (F' F^{-1})_{m,L} t ],
    G_ij = <exp(2R*U) Y_i Y_j>,  G_U,ij = <exp(2R*U) U Y_i Y_j>.

In the exact problem B_int is the interior energy form, so it is symmetric. In the channel model
it is only approximately symmetric. I symmetrize it and report the asymmetry.

## Exterior (R > R*): Ritz

Form: a_ext(u, v) = int_{R>R*} [grad u . grad v + 2(V - E) u v] dtau. It is coercive at E0 for
Dirichlet data. The Hylleraas gradient form is

    grad u.grad v = u1 v1 + u2 v2 + 2 u12 v12 + c1 (u1 v12 + u12 v1) + c2 (u2 v12 + u12 v2),
    c1 = (r1^2 + r12^2 - r2^2) / (2 r1 r12),  c2 likewise.

Trial function: sum_i t_i chi_i + sum_j b_j phi_j.
- Lift chi_i = exp(S) Y_i(q, s); on the sphere this equals tau_i. A second lift,
  exp(R* U(Omega)) Y_i(Omega) exp(-gamma (R - R*)), is used to show that the lift does not matter
  (and to make validation (i) nontrivial).
- Correction phi_j = (R - R*) * Sym[ r1^i r2^j r12^k exp(-za r1 - zb r2) ], i + j + k <= Omega.
  Exponent sets: (Z, kappa) symmetrized, where kappa = sqrt(-2(E+2)) = 1.3444 is the He+(1s) + e
  channel; plus a compact (zeta, zeta) set. Optional exp(S)-factor set.
- Quadrature: R = R* + t/beta (Gauss–Laguerre), phi Gauss–Legendre on [0, pi/4] (doubled),
  u Gauss–Legendre on [cos - sin, cos + sin]. All integrands are analytic on these cells.
- Exterior Steklov form: A_ext = A_tt - A_tb A_bb^+ A_bt (Schur complement). Basis near-dependence
  is handled by eigen-truncating A_bb. That is still a subspace, so the result is still variational.
  If Float64 limits it, I switch to Double64.

## Matching

S_L = sym(B_int) + A_ext. I solve S_L t = lambda M t, with M = R*^5 G_LL (the trace L2 Gram).
For the exact operators S >= 0 (it equals 2<Psi,(H-E0)Psi> for the glued, continuous Psi), and
S has a null vector at E0. With an exact interior and a Ritz exterior, S_h >= S, so
lambda_min >= 0 and lambda_min -> 0 from above. lambda_min(E) is decreasing in E. The E* where
it crosses zero is a by-product energy, and it is a check: E* - E0 should go to 0.
The null vector t0 gives a = F^{-1}[t0; 0], normalized by a_1 = 1. Equivalently, minimize
the glued quadratic form over (a, b). The Schur complement is just b eliminated.

## Reading off the constants

- Cusp-gauge seeds a_i (as in Codex's CSV).
- Physical a21 (total projection, as in the repo psi20): a21 = <psi20, v>/<v, v> with v = s.
  Since Phi_10 = 0, psi20 = Phi20 + U^2/2, so a21 = 2 a_2 + 2 <U^2/2, Y_2> (Y_2 = 2s). Analytically
  2<U^2/2, Y_2> = Z g (pi+4)/(6 pi) - g^2/8. I check this number by quadrature.
- Order 4 (two constants, Y_(2,0) and Y_(0,2)): physical <psi40, Y> = sum_k <U^k/k! Phi_{4-k,0}, Y>.
  Phi_{n,k} comes from the stored channel coefficients. I report both gauges.
- Orders 6 and 8 cusp-gauge seeds as "next constants".

## Validation

(i) g = 0, Z = 2, E = -4: Psi = exp(-2(r1+r2)), Phi = 1. All cusp seeds beyond a_1 are 0, and
    physical a21 = 0. The physical order-4 projections <(-2 sigma)^4/24, Y_(2,0)> are nonzero
    and known in closed form, which tests the gauge restoration. Use the plain lift so that the
    exterior Ritz is actually exercised. Expect lambda_min ~ 1e-13 and seeds ~ 1e-12.
(ii) Nonzero known constant: the same exact state with the plain gauge (U = 0, so Phi = Psi).
    Its seeds are total projections, e.g. <psi40, Y_(2,0)> != 0. Here the channel model has a
    V-truncation error, so this tests convergence in d_int to a known nonzero value.
    Not machine precision.

## Helium study and expectations

Scans: exterior Omega and exponents; L; d_int (beyond Codex's 24 if feasible); N; R* in {0.5, 0.75, 1, 1.5, 2}.
E0 +- 1e-6 for lambda_min. Figure: |a21 - ref| and the successive differences against Codex's channel sequence.

Expectation: the exterior Ritz should converge quickly, since there are no Fock logs outside R*.
The interior channel truncation (V, U, W coupling in finite harmonics, algebraic in d_int) may
well be the true bottleneck. Decoupling d_int from L and from the exterior settles which error
dominates in Codex's sequence. If d_int dominates, a21 will look like Codex's sequence at equal d_int.
