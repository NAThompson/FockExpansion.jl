# TASK_codex_19 (Codex, phase 9) — 2026-10-23
Read the newest PHASE 9 entries of COLLAB.md. Write codex/phase9/x19_*, append-only COLLAB.md / THEORY.md Log. Serial, memory < 2 GB.
1. Independently check Claude's psi41^(2) closed form (claude9/closed41.py -> closed41.pkl, coefficients claude9/coef5_ident.txt) with
   your own order-2 atoms/connection: exact or high-precision PDE check and purity; co-sign or dissent on "psi41^(2) in rank-5 module".
2. psi40 membership test (you own psi40). Data: claude9/spec40.txt (alpha theta C0..C4, 60 digits; pure psi40 at E=a21=0, Z^j sectors).
   Test whether C1, C2, C3 lie in polynomials(r1,r2,r12) x {order-2 module, classical 18 atoms, 4x23 panel states}, using your own
   panel evaluators (x16 family). Least squares with held-out points, as in claude9/fitbig.jl. Report residuals vs polynomial degree.
   A failure must be distinguished from too-low degree / missing denominators (try multiplying by row-clearing or discriminant factors).
3. Theory with the line-charge observation (COLLAB): T4[v0 x panel] and T4[v1 x panel] as integrals over (q, Feynman t) of the elementary
   zonal U-functions (repo feynman_k4.jl): identify the curve/surface of the remaining integration variables; is it fibred by the old
   E_(H,k) curves (iterated elliptic integrals), or does it require a new curve / K3 period? Use reasoning effort; concrete eliminations.
Mark PROVED / VERIFIED / CITED / CONJ / OPEN; co-sign/dissent.
