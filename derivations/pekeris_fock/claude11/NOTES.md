# claude11 — phase 13 working notes [C]

## Coordinates
Spectral grid: X = alpha - beta, Y = alpha + beta, u = cos alpha = x1, v = cos beta = x2, tx = 2X/pi, ty = 2(Y - pi)/pi.
Physical square D0: |X| <= pi/2, pi/2 <= Y <= 3pi/2. All three Coulomb radicals are entire in (alpha, beta).
The weighted operator is SEPARABLE: w (Lambda^2 - c) ~ M_X - M_Y, M = d(cos d) + (c/16) cos (claude9/fock_spectral_copy.jl).
Degenerate lines cos X = 0, cos Y = 0 = the collinear circle w^2 = -cos X cos Y = 0 (all its sheets). Each degenerate line has
exponents (0,0) for M: analytic or log.

## Continuation method (extend.jl)
Y stays physical (L_Y eigenbasis of sfock); X is continued across the physical edge X1 = side*pi/2 to the next degenerate
line X0 = side*3pi/2. Per Y-mode the PDE is a second-order ODE in X. Element 2 = [Xmid, X1] (Chebyshev-Lobatto, analytic at
X1, value at X1 = physical spectral value, 60 digits); element 1 = [X0, Xmid] with psi = sum_m A_m(X) l^m, l = log(-cos X),
A_m analytic at X0 (conditions A_m(X0) = 0 for m >= 2 forced by analyticity of the lower equations; the free homogeneous parts
of A_0, A_1 fixed by C^1 matching at Xmid). The whole hierarchy (k <= 4, all p, E = 0, a21 = a40 = a42 = 0, Y4 harmonic
additions copied from the physical solve) is continued order by order; sources on element 1 carry their own l-powers.
Output: A_m(X0, Y) = log coefficients of the continued branch on the far edge, per Z-sector (Vandermonde over Z = 0..4).
Self-consistency: no analyticity across interior lines is ASSUMED beyond what element convergence shows; a hidden singular
line inside an element (alpha = pi, beta = 0 lines cross element interiors) would destroy n = 64 vs 80 agreement.

## Sheets
side = -1: path crosses alpha = 0 (u = 1) and beta = pi (v = -1); far edge = arc (e1, -e2). RESULT: no logs through psi40
  (all |A_1| < 1e-29 at 128 bits, Z = 2), i.e. this branch is analytic at that arc.
side = +1: path crosses alpha = pi (u = -1) and beta = 0 (v = 1, the e-e Coulomb line); far edge = arc (-e1, e2) on the
  sheet with r1 and r12 flipped; alpha = (3pi/2 + Y)/2 in (pi, 3pi/2), u = cos alpha in (-1, 0), v = -sin alpha.
  RESULT: log degree exactly 1 at every order k = 2, 3, 4 (N^2 = 0 numerically): psi20 already has a log there.

## Fingerprint (co-signed by Codex x25 with isogeny correction)
Type II limit of T5 at the circle; Im N = H1(E_q), E_q: a^2 = uv q^4 - 2(1-v^2) q^2 - uv, 2-isogenous to Legendre
lambda = (u-1)/(u+1); normalized periods p = [u(u+1)]^(-1/2) F(lambda) satisfy
  [4u^2(u^2-1) D_u^2 + 4u(3u^2-1) D_u + (4u^2+1)] p = 0.
C3 log coefficients at a generic circle point are mixed Tate. So an elliptic (E_q) constituent in the far-edge log coefficient
of the psi40 sector germ proves psi40 is not in C3. Test: guess the ODE of B(alpha) = A_1(X0, Y) and look for this factor.
