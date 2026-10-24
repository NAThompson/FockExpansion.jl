# TASK_codex_2 (Codex) — 2026-10-21

Read COLLAB.md first (new Claude entries at the end). Write only in codex/, codex_*, NOTES_codex.md, REVIEW_2_codex.txt,
and append-only entries in COLLAB.md. You may now read all of Claude's files (sym/, julia/, ANSATZ.md when present).
Julia: `julia` (no extra packages needed for plain BigFloat/LinearAlgebra); Python: venv/bin/python.

## 1. Stage 3 write-up (you own it): codex/codex_minimality.md
Turn NOTES_codex.md into a precise statement:
 - the block-banded shell recurrence (sizes, bandwidth, for k=0 and for the log tower with K sectors);
 - the cleanest *correct* formulation of the decay/boundary condition as a minimality / admissible-subspace condition
   (spectral/form-domain version, matrix-continued-fraction / Schur-complement version for nonreal spectral parameter),
   with each claim marked PROVED (give the proof or exact citation), CITED, NUMERICAL, or CONJECTURED;
 - a cheap numerical illustration if you can: e.g. shell norms ||a_N|| of the truncated 33-term eigenvector (julia/pekeris33.jl
   builds the pencil) vs N, and what a non-physical (wrong-eps or forward-iterated) solution does; does physical ⇔ algebraic decay,
   nonphysical ⇔ geometric growth hold empirically?
 - Fock constants as functionals: which sums converge, which need finite parts.

## 2. a21 as a global functional of the computed states (optional but valuable)
Claude's local read-off converges slowly (0.4724 at N=13 vs 0.476748; julia/a21_readoff.out). Apply the global finite-part identity of
/Users/nathompson7/Documents/Codex/2026-09-05/how-much-compu/outputs/a21_global/FINITE_PART.md
(a21 = -1/8 FP∫ h(R)/R dR - b/8, h = 2 P_v[C psi]/R - 2E f) to the exported states julia/states/fixedlog_N{6..13}.txt
(format in julia/dump_states.jl; Psi(0)=1; eps=1.704; use E = E_variational from the header or E0 — say which) and report a21(N).
Check your implementation on a manufactured case first. Report in codex/codex_a21_functional.md.

## 3. Review of Claude's stage-4 numerics: REVIEW_2_codex.txt
Critique julia/pekeris33.jl, julia/galerkin_log.jl (moments, operator application, variants), julia/a21_readoff.jl, the plot,
and the claims in COLLAB.md. Point out errors, unjustified claims, or better comparisons (e.g. equal-basis-size comparison).

## 4. Co-sign / dissent
At the end append to COLLAB.md a dated "Codex position on the final claims" entry: CO-SIGN or DISSENT, per claim, for the claims in
the Claude entries and ANSATZ.md.
