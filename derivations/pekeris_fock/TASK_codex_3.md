# TASK_codex_3 (Codex, phase 2) — 2026-10-21
Read COLLAB.md (PHASE 2 section) first. Write only in codex/phase2/, NOTES_codex.md, REVIEW_3_codex.txt and append-only
COLLAB.md entries. Do not read claude2/ until your A and B are written (sync points); julia/ may be read for harness reuse only after B.
Read-only: /Users/nathompson7/FockExpansion.jl (src/, derivations/k4/README.md), earlier workspace outputs/.
Python: venv/bin/python (sympy, mpmath). Julia: `julia --project=/Users/nathompson7/FockExpansion.jl` gives the package (BigFloat works for psi20).

A (independent). For psi10, psi20, psi21, psi31, psi30, psi40, psi41, psi42 (repo closed forms / numerical evaluators):
  write R^n psi_{n,k} in (r1,r2,r12) (or u,v,w) and separate (i) polynomial part, (ii) non-polynomial remainder in closed form.
  Give a minimal list of "transcendental atoms" (degree-0 homogeneous functions or R) through order 4, and say which psi's are only
  available numerically. Explain Claude's earlier finding: the degree-2 log-free part is not a quadratic polynomial (rms residual 6e-3).
  Output codex/phase2/codex_atoms.md (+ scripts). Mark proved / verified / cited / conjectured.
B (independent). Order-2 extension: find the smallest atom set containing the order-2 non-polynomial content that is CLOSED under
  differentiation in (u,v,w) modulo multiplication by rational functions (holonomic / D-finite closure; radicals such as
  R = sqrt(r1^2+r2^2) count as algebraic atoms). Then derive the recurrence for Psi = e^{-sigma/2} sum_k log^k sigma sum_i g_i F_{k,i}:
  clearing factor, offsets/term counts per (sector -> sector) block, coupling structure (triangular in transcendental weight?).
  Machine-verify by direct operator application on random finite combinations (same style as sym/verify_recurrence.py, which you may read).
  Say plainly if finiteness fails and why. Output codex/phase2/codex_order2_recurrence.md (+ json in the format of sym/recurrence.json,
  keyed by "(target_atom,source_atom)").
D. Does the atom-extended ansatz change the admissibility <=> minimality picture of codex/codex_minimality.md? (short note,
  codex/phase2/codex_minimality2.md)
At the end, append a dated COLLAB.md entry summarising A, B, D and listing files. A later task will ask for the cross-check and review.
