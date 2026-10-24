# TASK_codex_4 (Codex, phase 2: sync 4 cross-check, review, co-sign) — 2026-10-21
Read COLLAB.md (all PHASE 2 entries) first. You may now read everything in claude2/. Write only in codex/phase2/,
REVIEW_4_codex.txt, NOTES_codex.md and append-only COLLAB.md. Do not edit Claude's files.

1. Sync 4 — recurrence cross-check. Claude's order-2 system uses a different (non-minimal, 8-atom) basis:
   {Om = logs - DL, al, beP = P*beta, lR = log(R/s), la = log((s+r12)/R), lb = log((d+r12)^2/L5), one, P}
   (claude2/order2_atoms.py, chi_letters.py; Gamma table; claude2/order2_recurrence.json with keys "('target','source')" in
   sigma2P/C1/C2, offsets as before, row factors in 'D', E=-eps^2, Z=2, lambda=1).
   (a) Independently check Claude's recurrence: e.g. run your direct-application verifier against claude2/order2_recurrence.json
       (Claude's own claude2/verify_order2.py is running; do not rely on it). (b) Relate the two atom bases explicitly
       (e.g. U, V, L in terms of Om, beP, P, la, lb) and check consistency of the connections. (c) Compare counts block by block
       where comparable; explain differences (row factors). Write codex/phase2/codex_crosscheck2.md.
2. Review Claude's numerics (claude2/atoms_ritz.jl, run_atoms.jl, run_atoms.out, run_atoms4.out, a21_global_atoms.jl,
   a21_atoms.out, the (s,x,y) weak form, moment-table assembly, quadrature, conditioning, the psi30 atom handling and the
   N>=8 variational violation). Run your own a21 functional on claude2/states/atoms2_N*.txt (format: "n p a b atom coef",
   Psi = e^{-eps s} sum coef s^n log^p s x^a y^b atom(x,y), atoms chi, lR=log(sqrt((1+x^2)/2)) as in atoms_ritz.jl) and say
   whether the a21 digits can be trusted and to what uncertainty. Is convergence "faster than algebraic"? Write REVIEW_4_codex.txt.
3. Minimality (D): confirm or update codex/phase2/codex_minimality2.md in light of the numerics.
4. Co-sign or dissent, per claim, in a dated COLLAB.md entry "Codex position on phase-2 final claims".
