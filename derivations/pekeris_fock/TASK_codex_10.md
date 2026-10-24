# TASK_codex_10 (Codex, phase 5, Goal 2 support) — 2026-10-22
Read COLLAB.md (PHASE 5) and THEORY.md. Write codex/phase5/goal2_*, NOTES_codex.md, append-only COLLAB.md and THEORY.md Log. Keep memory small.
1. Compute the rank of the smallest differential module (over the algebraic coefficient field you need; also report the rational rank
   after restriction of scalars if feasible) containing {1, chi, R^3 psi30}, using your explicit connection (codex/phase3/panel_connection.py,
   classical module in structure.md). Suggested method: cyclic rows r_0 = l (the linear form giving psi30), r_{alpha+e_k} = d_k r_alpha +
   r_alpha C_k, evaluated at generic points with exact Taylor coefficients of C (e.g. 2D Cauchy-integral/FFT on small complex polydiscs,
   or truncated power-series arithmetic), rank of {r_alpha(p)} until it stabilises at two consecutive orders, at >=2 generic points.
   Do it per piece (classical; each elliptic panel; N0 vs N2) and combined. Then prune: give an explicit minimal generator set.
2. Comment on Claude's single-valuedness position in COLLAB.md PHASE 5 (agree/dissent, with reasons), and on which radicals survive.
3. Later (when Claude's claude5/goal2_* recurrence exists): independently re-derive the closure and spot-check coefficients.
