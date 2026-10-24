# TASK_codex_16 (Codex, phase 8 continued) — 2026-10-23
Read COLLAB.md PHASE 8 (new Claude entries). Write codex/phase8/x16_*, append-only COLLAB.md / THEORY.md Log. Small serial jobs.
1. Your (H,k) 55-state rational family covers all four panels via rational maps (panel 1: H = r1^2/r2^2, k = (r12^2-Q)/r2^2; panel 2 = exchange;
   panel 3: H = r12^2/e, k = 2(r1^2-r2^2)/e; panel 4 = exchange). Export the DISTANCE-coordinate connections (d/dr1, d/dr2, d/dr12 matrices,
   exact rational entries) and the exact output rows (each panel's contribution to R^3 psi30, i.e. R^3 (c0_j N0_j + c2_j N2_j) as a combination
   of atoms with coefficients in Q(pi)(r1,r2,r12) or with an explicit overall algebraic factor such as R or sqrt2 — say exactly) for ALL FOUR
   panels, in a simple text format: "panel|axis|row|col|numerator_sympy|denominator_sympy" and "panel|output|atom|coefficient_sympy".
   Also state reachable atoms from the outputs (drop unreachable/duplicate states) and sharp row clearing factors.
2. Cross-check Claude's independent panel presentation claude8/panel1_connection.json (format in claude8/panel1_solve.py; states in
   claude8/panel1_fp.py docstring: J, K, S0..S3, J5..J8 (poles i,-i,ih,-ih), Js, e, f0; index = state*3 + jet; entries re/im parts over Q(h,Y2))
   against your family at fresh points: e.g. compare the derivatives of the physical N0, N2 (both presentations) with respect to h and Y2.
3. Admissibility: continue toward the minimal weight / whether the exact psi30 panel output lies in the closure of polynomial x weighted atoms.
Append dated entries; co-sign/dissent.
