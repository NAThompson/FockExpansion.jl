# TASK_codex_8 (Codex, phase 4) — 2026-10-22
Read COLLAB.md (PHASE 4). Write only codex/phase4/, REVIEW_8_codex.txt, NOTES_codex.md, append-only COLLAB.md. Machine has 16 GB: keep jobs small.
1. Independently (no Claude source) implement the Fock-constrained trial space described in PHASE 4: seed f0 with the exact order<=2 Fock
   jet (Z=2, E=E0, psi00=1), a free D*e^{-eps s} term, plus degree>=3 terms of your choice (e.g. polynomials, log s, log^2 s, chi and
   log(R/s) atoms, or your minimal atoms). Solve the Ritz problem at fixed eps=1.704 for small/moderate N and read a21 = a/c0.
   Verify independently that the constructed seed has exactly the Fock jet (e.g. symbolic/numerical R-expansion along rays).
2. Cross-check against Claude's claude4 results when they appear (claude4/run4*.out, states), at matching N if possible.
3. Review: is the a21 readout from a constrained Ritz vector a consistent estimator (does it converge to a21 as the space becomes dense in
   the form domain? what is its error order relative to the energy error?). If cheap, propose/implement a residual-corrected estimator.
4. Append a dated COLLAB.md entry; co-sign/dissent on Claude's phase-4 claims once they are written.
