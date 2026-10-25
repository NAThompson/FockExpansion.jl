# TASK_codex_21 (Codex, phase 11) — 2026-10-23
Read COLLAB.md PHASE 10/11 and codex/phase10/x20_report.md. Write codex/phase11/x21_*, append-only COLLAB.md / THEORY.md Log. Serial, < 2 GB.
Goal: decide Conjecture A (rank-5 T5 constituent of the physical psi40, in sectors Z^1, Z^2 (e-e x panels 1,2) or Z^3 (nuclear x panels 3,4)).
Exact data available: Claude's order-5 operator on x2 = 1/3 (claude10/op_e2_Q.sobj, printed in COLLAB), 60-digit psi40 sectors
(claude9/spec40.txt; evaluator claude9/psi40_big.jl), old exact connections (phases 2, 7, 8).
1. Route 1 (you lead): mixed-Hodge / Gauss–Manin. For the actual relative period (with the atan(q/h(y))^2 density as two extra parameters
   and the chain boundaries), determine Gr^W pieces and whether H^2_tr(S_e2) occurs. Look for a cheap *invariant* that decides nonvanishing,
   e.g. a residue / leading-coefficient pairing <density class, transcendental cycle> that is an explicit computable number or function;
   or a symmetry (electron exchange, the s -> h/s antisymmetry of the panel integrand, charge sums) that forces the pairing to vanish.
2. Audit Claude's PHASE 11 framing (Q, P-harmonic solutions, quotient-vs-sub tests) and the local-exponent shortcut (analyticity across
   the physical edges vs T5 logs there): does analyticity force the T5 projection to vanish, or only constrain it? Give a proof either way.
3. Propose the cheapest decisive computation if neither settles it, with a cost estimate. Mark PROVED/VERIFIED/CITED/CONJ/OPEN.
