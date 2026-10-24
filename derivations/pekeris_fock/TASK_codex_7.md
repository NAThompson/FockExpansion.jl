# TASK_codex_7 (Codex, phase 3 review + co-sign) — 2026-10-22
Read COLLAB.md (PHASE 3 entries). You may read claude3/. Write only REVIEW_7_codex.txt, codex/phase3/, NOTES_codex.md, append-only COLLAB.md.
1. Review claude3/atoms3.jl (min5 atoms: Lt, alp = alpha - pi/2, V, U from global letter forms; parity handling with odd x-powers for
   exchange-odd atoms; BigFloat psi30 atom via the repo, rtol 1e-50 at 384 bits), run3.jl, check3.out, param_check.out and all
   run3_*.out, a21_*.out. Independently spot-check at least one energy (e.g. min5 N=8 or 9) and one a21 value with your own code.
   NOTE: the machine has 16 GB; Claude's min5psi N=10,11 job is running — keep your jobs small (single BigFloat job, modest N).
2. Judge the explanation of the old sub-E0 energies (Float64 contamination; parameters irrelevant to the upper bound) and the
   statement that the psi30-atom space is independent of E0/a21 (parameter parts are cubic polynomials).
3. Slope: is there evidence that the log-log slope steepened beyond -8 with the minimal atoms (claude3/run3_min5.out; local slopes about
   -12.5, -13.8, -15.1, then -10.5 at N=10->11)? Is the semilog curve still convex? Honest a21 digits?
4. Co-sign or dissent, per claim, in a dated COLLAB.md entry "Codex position on phase-3 claims". Also restate your own phase-3 status
   (psi30 closure: local finite, minimal rank and global interface open; psi41^(2), psi40 reachability).
