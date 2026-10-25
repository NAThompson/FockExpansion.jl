# TASK_codex_18 (Codex, phase 9 / order 4) — 2026-10-23
Read COLLAB.md PHASE 9. Write codex/phase9/x18_*, append-only COLLAB.md / THEORY.md Log. Small serial jobs, memory < 2 GB.
Repo (read-only): /Users/nathompson7/FockExpansion.jl — src/fourth_order.jl (psi40 numerics), src/green_k4.jl (S^3 Green's function G,
iterated G2), src/feynman_k4.jl (psi41^(2) 1-D integral), derivations/k4/README.md (reductions: a21 sector, E sector, Z^0/Z^4 closed forms).
You own psi40 (Claude owns psi41^(2); we cross-check each other afterwards).
1. psi40 particular solution: (Lambda^2-32) psi40 = 12 psi41 + 2 psi42 - 2 V psi30 + 2 E psi20. Split by powers of Z and by sector (the
   README already reduces the a21 and E sectors and gives Z^0, Z^4 closed forms). Derive a Green's-function (or lower-dimensional) integral
   representation for the hard part, i.e. the Z^2/Z^3 terms V x (Clausen and elliptic-panel parts of psi30), V x chi, and Z^2 psi41^(2).
2. Identify the function class: does it reuse the order-3 elliptic curves (iterated integrals / elliptic polylog of one weight higher) plus
   polylogs, or does a genuinely new class appear (higher-genus periods, K3-type, new elliptic curves)? Find the singular loci (at order 3
   the panels came from coincidences of the quartic's roots); list any new branch loci / curves with their j-invariants or genus.
3. If feasible: measure the holonomic rank (structured blocks, not cyclic vectors) of the hardest piece at exact finite-field points.
4. Attack Claude's (H3) holonomy-by-induction argument in PHASE 9: is it a proof? What hypotheses are missing (integration domains, the
   pure-inverse/harmonic projection, log R sectors, Coulomb singular sources)? Does anything bound the rank growth with k?
Mark claims PROVED / VERIFIED / CITED / CONJ / OPEN. Append a dated COLLAB.md entry; co-sign/dissent on Claude's hypotheses.
