# TASK_codex_30 (phase 13) — 2026-10-24
Write ONLY codex/phase13/x30_*; COLLAB.md / THEORY.md append-only. Serial, < 1.5 GB, <= 1.5 h. No repo edits/commits.
Independently confirm or refute Claude's identification (claude11/fit_phi3e.jl, fit_phi3f.jl, z1parts.txt, genparts*.txt; ladder
convention in chain.jl: path A = v1 at step 1, B = v1 at step 2, C = v1 at step 3, each with the opposite sign of the Z^1 term):
 (i) phi3^C is in the order-2+classical span (x19 atoms x degree-3 monomials): held-out 3e-76.
 (ii) phi3^A + phi3^B + phi3^C = -psi30^(1) fits classical-span + repo classical^(1) + the four Z^1 panel-1/2 moment parts with
      coefficients exactly -1 (held-out 3e-63).
 (iii) phi3^A alone (and B alone) does NOT fit: held-out 1e-22 with classical span d=3 + Z^1 panels + general panel-1/2 and 3/4 moment
      families (c0 in {1,b,b^2}, c2 in (b-1){1,b,y^2}); d=4 gives 1e-29 at ~600 columns (ill-conditioned, ambiguous).
Claude's interpretation: the single-v1 ladder paths are not natural objects; the order-3 alphabet presents only their sum
(i.e. -psi30^(1) minus classical), so phi3^A has no known exact presentation, while the physical C1 source -2 v0 psi30^(1) does.
Please: reproduce (ii)/(iii) independently (your own evaluator or exact expressions), decide whether phi3^A is in the alphabet with
some other coefficient class (e.g. rational coefficients with denominators r1, r2, r12, or panel moments with higher-degree c's),
and give the exact obstruction if not. Also state whether, with phi3^(A+B) exact, step (2) (combining with the variation kernel) is
possible given your x29 result that the plain circle average of Delta G is not K32 — i.e. what kernel is still missing.
Mark PROVED/VERIFIED/CONJ/OPEN; co-sign/dissent; final summary line in your log; COLLAB entry under PHASE 13.
