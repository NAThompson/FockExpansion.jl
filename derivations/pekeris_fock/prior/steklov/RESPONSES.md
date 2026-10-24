# Responses to Codex reviews

## Review 1 (REVIEW_1_codex.txt, plan stage)

1. Measure. Agreed. My harmonics are orthonormal on my (phi, u) rule to 7e-15.
   My U, C, W matrices agree with Codex's cuspmodel to 1e-14 (U), 2e-14 (C) and 4e-13 (W).
   Refining nquad from 40 to 120 moves them by at most 2e-14 (test_model.jl, check_a21conv.jl).
   Integrands are evaluated only at interior Gauss nodes. The direction cosines are multiplied by the
   r1 r2 r12 weight before summation.

2. Cusp PDE and recurrence. Agreed. D is built from the identity with D[:,1] = 0.
   I store c_{n,k} only for n <= 6 and accumulate F, F', F'' on the fly. N defaults to
   max(44, 2 d_int + 28). An N check is in the stage-2 runs.

3. Trace compression. Agreed: B_LL = R*^5 [G_U,LL + G_{L,m} (F' F^{-1} J)].
   That is what interior_setup computes (DtN = Fp * (F \ Ein), full m x nL). Seeds are F \ [t; 0]
   (a linear solve, not an inverse). The suggested L study was done.
   The interior-only 2s x 2s test (validate_2s2s_split.jl) shows the forced-zero high trace costs
   little once L >= 12. The error there is dominated by d_int. I did not implement the
   total-form Schur elimination of high traces.

4. Positivity. Accepted. The PSD and monotonicity statements hold only for exact operators.
   The channel interior is not an exact solution, so the computed flux form is neither symmetric
   nor guaranteed nonnegative. Changes made:
   - I do not symmetrize. I solve the nonsymmetric problem S t = lambda M t and take the eigenvalue
     of smallest modulus. Symmetrizing visibly destroyed the exact null vector in validation (i):
     lambda went from -3.9e-7 to 1e-14 when I removed it.
   - lambda_min is reported only as a diagnostic, never as a bound. At finite L, S_L has no exact
     null vector even with exact operators.
   - I did not implement the physical interior energy matrix, and record it as an open item.

5. a21 and order-4 constants. Agreed. Checked by quadrature: 2<U^2/2, Y_2> = 0.6327465149117193,
   and the formula gives 0.632746514911721. Order-4 projections use powers of the full angular U
   evaluated at quadrature points (gauge_moments), not powers of the truncated matrix.

6. Exterior numerics. Taken further than suggested. Float64 Gram matrices of the Hylleraas set
   lost needed directions: a21 moved by 1e-5 as the eigen cutoff went from 1e-12 to 1e-16
   (ext_test2.jl). The exterior is now assembled in Double64:
   - Hylleraas-Hylleraas block: analytic in r12 (all integrands are polynomials in r12 after
     the S-state weight), checked against 3D quadrature to 7e-16.
   - Lift and channel blocks: 3D Gauss rules.
   - Radial rule: Legendre panels at R*+{0.25, 0.75, 2, 5}, then Laguerre.
   Quadrature is converged to ~1e-15 in a21 (test_quad2.jl). Columns are Jacobi-scaled. Cholesky
   is used, with an eigen-truncation fallback at 1e-28. Retained ranks are reported.

7. Validation. Added the 2s x 2s nonzero cusp-gauge seed test (interior only, since E = -1 is in
   the continuum): a_(2,0) -> -1/32 with error 1.4e-8 at d = 24 and 6.9e-9 at d = 28, decaying
   like d^-4.8. a_2 and a_(0,2) are zero to 1e-15. The plain-gauge model uses D = -2C, U = W = 0,
   which is probe.jl's ungauged equation, as the review requires.

8. Selection rule. Done (codex_selection.jl -> codex_selection.csv). With the same channel
   interior and exterior, the constrained least-squares a21 and the mismatch-eigenvector a21 differ by
   3.6e-8 at d = 4, 1.6e-10 at d = 8, 2e-12 at d = 12, and <1e-12 at d = 16. The selection rule is
   negligible next to every effect studied here.

## Review 2 (REVIEW_2_codex.txt, results stage)

1. Angular quadrature at high L. This was the priority finding. I reran the production exterior with
   the refined rules (nphi, nu) = (40, 40) for L = 8, 10, 12, 14 and (56, 52) for L = 12
   (quad_hiL.jl -> quad_hiL.log). The refined runs use separate cache files. Change in a21 against
   the (24, 20) rule, B10P6, R* = 1:
   L=8: < 1e-13; L=10: < 1e-13; L=12: 1e-13 at (40, 40) and at (56, 52); L=14: 2e-12.
   lambda changes by <= 7e-12, and the a3/a4/a5 seeds by < 1e-12. Codex's Gram-matrix
   discrepancy is real, but it does not reach the matching result at the reported digits.
   I did not refine the radial rule at L >= 12. It was checked only at L = 4, to 1e-15.
   This is listed as an open item.
2. d tail. Agreed on the sign: the tail is +2.0e-9, added rather than subtracted, and it is a
   model-conditional extrapolation, not a bound. I also accept the rewording. Interior truncation
   drove the d = 24 -> 48 change (+5.6e-8). At d = 48 the remaining d-tail (2e-9) is smaller than
   the last L step (-6.2e-9) and the R* spread.
3. Final value. I now compare the L steps within one basis (B10P6): 12 -> 14 is -6.16e-9 at d = 48.
   R* spread is used as one combined consistency envelope, not added in quadrature.
   Reported value: a21 ~ 0.47674788, with a provisional uncertainty of about 1e-8 (not a bound,
   and no claim of eight correctly rounded decimals).
4. Higher seeds. The R* = 0.75 offset in a3 (+1.1e-7) is not explained. Neither the d nor the L
   steps at that radius account for it, so I report a3 to about +-1.2e-7 and say so.
5. lambda. Wording accepted. lambda is a combined matching residual. It falls with L, and the
   exterior basis moves it by about 2.5e-9.
6. Channel comparison. Accepted. +2.8e-8 is an observed equal-d difference between configurations
   that differ in R* and trace restriction, not a measurement of the old exterior error. The
   supported claim is the practical one: Steklov lets d go to 48 cheaply.
7. Code. match() now errors on a complex selected eigenvalue or eigenvector, instead of warning and
   continuing. Not done: a higher-precision recurrence at d = 48, and analysis of the 15 directions
   that B12P6 drops. B12P6 and B10P6 differ in a21 by only 7e-10 at L = 12, d = 48.
