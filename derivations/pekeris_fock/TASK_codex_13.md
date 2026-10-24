# TASK_codex_13 (Codex, phase 7) — 2026-10-22
Read COLLAB.md (PHASES 5–7) and THEORY.md. Write codex/phase7/*, NOTES_codex.md, append-only COLLAB.md and THEORY.md Log. Small serial jobs.
1. Per-block connections. Blocks: (a) the classical part C (your 18-generator module; rank 12); (b) each elliptic panel j=1..4 (N0_j, N2_j with
   their 13/39/55-state Pfaffian system, pulled back to (t,q) or distances). For each block give an atom list (prefer a basis whose connection
   coefficients are RATIONAL in (t,q), or adjoin the minimal radicals explicitly, e.g. sqrt(D), quartic roots, and restrict scalars), the
   connection matrices, and the maximal numerator/denominator degree of the entries (exact symbolic if possible; else exact finite-field line
   degrees as in claude6/ratrec.py). Target: tens, not hundreds. Also report the per-block rank after removing redundancy.
2. Global single-valuedness: the panel endpoint u = min(1,h) switches at b_j = 0 (b_1, b_2 = +-(r1^2-r2^2)/R^2 vanish on the exchange line;
   b_3, b_4 = -v vanish on v = 0). Either prove the per-block atoms (in a rational-coefficient basis) extend analytically across these interior
   lines (e.g. by showing the switching is only in the presentation, or by re-basing the panel integrals with a fixed endpoint), or give the
   minimal piecewise structure and the matching relations.
3. Admissibility: prove (or give the needed vanishing weights for) every block atom x polynomial x e^{-sigma/2} in the Coulomb form domain.
4. When Claude's claude7/ recurrence appears, independently derive at least one block's recurrence and cross-check counts/coefficients.
Append dated entries; co-sign/dissent on Claude's phase-7 claims.
