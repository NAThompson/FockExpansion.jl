# TASK_codex_5 (Codex, phase 3) — 2026-10-22
Read COLLAB.md (PHASE 3) first. Write only in codex/phase3/, NOTES_codex.md, REVIEW_5_codex.txt, append-only COLLAB.md.
Read-only: /Users/nathompson7/FockExpansion.jl (src/FockExpansion.jl psi30, classical(), elliptic moments N0/N2, panel code;
derivations/, paper/ if useful), the earlier workspace outputs/ (e.g. compact_elliptic_moments/, classical_polylog_reduction/), claude2/.

1. psi30 structure. Express R^3 psi30 (Z, E, a21 general) as: polynomial + Clausen/log/angle part + elliptic part, in distances.
   Give either a closed form or a holonomic (D-finite) characterization: a finite set of functions g_i, closed under d/dr_k modulo
   rational (or algebraic, with explicit radicals) coefficients, containing psi30 — i.e. a rational connection d g = C_k g as in your
   order-2 work — with initial/normalization data fixing psi30 among the solutions. The elliptic moments N0, N2 are integrals in q of
   algebraic*arctan-type integrands; derive the differential system they satisfy in the outer variables (b, y) or distances.
2. Minimal closed atom set (order 2 + order 3 together) and whether the Laguerre recurrence stays finite (clearing factors, counts,
   shell bandwidth). Machine-verify the connection numerically against the repo's psi30 (Float64 / BigFloat evaluator).
3. Reachability: can psi41^(2) (Z^2 part of psi41) and psi40 be treated the same way within this step? Give a reasoned yes/no.
4. Append a dated COLLAB.md entry with results and files. Claude will then cross-check.
