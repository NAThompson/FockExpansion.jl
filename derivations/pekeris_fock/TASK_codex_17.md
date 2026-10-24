# TASK_codex_17 (Codex, phase 8 cross-check) — 2026-10-23
Read COLLAB.md PHASE 8. Write codex/phase8/x17_*, append-only COLLAB.md / THEORY.md Log. Small serial jobs (machine has 16 GB; Claude's checks
are running — keep memory < 2 GB).
1. Independently derive the cleared Hamiltonian operator for at least two blocks of one panel (e.g. panel 3, one diagonal and one off-diagonal
   sigma2P block, and one C1 block) from your own distance connection and the distance Hamiltonian, and compare with
   claude8/panel3_terms.txt (format: block_key|du dv dw|i j k a b (u^i v^j w^k eps^a PI^b)|num|den) up to the row factor in
   claude8/panel3_rowfactors.json. Optionally spot-check a few stencil coefficients by your own Laguerre expansion.
2. State the final admissibility position for the order-3 atoms (classical unweighted; panels with delta^4 or better) and whether the
   Goal-1 backward theorem applies to poly x {order-2, classical, weighted panel} spaces.
3. Co-sign or dissent on Claude's phase-8 claims in a dated COLLAB.md entry.
