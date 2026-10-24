# TASK_codex_14 (Codex, phase 7 cross-check + remaining items) — 2026-10-22
Read COLLAB.md PHASE 7 entries. Write codex/phase7/x14_*, append-only COLLAB.md and THEORY.md Log. Small serial jobs.
1. Cross-check Claude's classical-block recurrence: claude7/classical_ops.json holds the cleared polynomial operators per block
   ((name, target, source) -> {deriv: polynomial in u,v,w,epsilon,PI}), atom names g0..g17 = your classical atom order, Pg13/Pg14/Pg17 = P*a,
   P*p, P*1; row factors in "D". Independently compare with your own derivation: at least the invariant subblock {1, a, a^2} (your
   angle_recurrence at eps=1) and one block involving L/T atoms, up to the row factors. When claude7/classical_stencils_*.json exist, spot-check
   stencil coefficients at a few equation indices.
2. Panels: estimate the size of the explicit global recurrence for the elliptic blocks from your 55-state (h,z) connections pulled back to
   distances: reachable sector count after restriction of scalars (K(A,z1) for panels 1/2, K(P,A,z3) for 3/4), row-clearing degrees and
   rough stencil sizes. Say whether a full expansion is feasible, and whether a smaller rational presentation (e.g. the 13-atom cyclic per-panel
   basis with line degree ~48, or symmetric functions of the pole integrals) would be better.
3. Step 6: confirm (or not) that the enriched Goal-1 backward theorem extends to polynomial x {classical atoms, panel atoms with the needed
   weights} spaces — admissibility + a bounded positive normalization suffice?
Append dated entries; co-sign/dissent on Claude's phase-7 claims.
