# TASK_codex_26 (phase 13, joint) — 2026-10-24
Write ONLY codex/phase13/x26_*; COLLAB.md / THEORY.md append-only. Serial, < 2 GB, 1-2 h. No repo edits/commits.
Read claude11/NOTES.md, claude11/extend.jl, claude11/chain.jl, claude11/claim_chain.md (Claude's new numerical results).

Claude's results (please audit and co-sign/dissent):
 R1. side = -1 far edge (X0 = -3pi/2; path crosses u = 1, v = -1): continued psi20..psi40 have NO log there (analytic branch).
 R2. side = +1 far edge Xb = 3pi/2 (path crosses u = -1 and v = 1; arc (-e1, e2) on the sheet with r1, r12 flipped): log degree
     exactly 1 for psi20, psi30, psi40 (sectors Z^1..Z^3), N_b psi40 converged n64/n80 to ~1e-37 (n80 accurate ~1e-50).
 R3. Iterated variations on the strip J x Y_phys (J = [pi/2, 3pi/2]): N_a N_b psi20 = N_a N_b psi30 = 0, N_a N_b psi41 = 0;
     N_a N_b psi40 is nonzero ONLY in the Z^2 sector and on arc a (physical edge X = pi/2, Y physical; u = cos alpha, v = sin alpha)
     equals 1/9 - (7/45) sin 2 alpha = 1/9 - (14/45) u v to 52 digits (N = d/dl, l = log(-cos X)); N_b N_a N_b psi40 = 0.
     The guessed ODE of N_aN_b psi40 (Z^2) on arc a was exactly [(1+t^2)^2 D^3 + 6t(1+t^2) D^2 + (66+6t^2) D], t = tan(alpha/4).

Your workstream:
 A. Audit extend.jl/chain.jl (math of the separable reduction w(Lambda^2 - c) ~ M_X - M_Y, the l-ansatz equations, the analyticity
    conditions A_m(X0) = 0 for m >= 2, matching, the claim that N commutes with the hierarchy, and that only V_0 carries the Y4
    additions). Reproduce one number independently if cheap (e.g. psi20's N_b at one Y from the closed-form psi20 of the repo
    src/ or an independent ODE integration of one Y-mode).
 B. KEY QUESTION: are the loops (a, b) T5-blind? Using the T5 local system (exact L5 on x2 = 1/3 and/or the 2-D connection; the
    claude10 50-digit monodromy data; reflections at u = +-1 and the v = 1 local monodromy), decide whether, for T5, the type-II
    images satisfy W1(a) = W1(b) (then N_a N_b = 0 on T5 and this chain cannot see T5) or not (then a T5 constituent in a generic
    position would typically give a non-terminating chain N_b N_a N_b ... != 0). Make the dependence on the position of T5 in
    D.psi40 (sub / middle / quotient) precise.
 C. Explain (or derive) the rational double variation 1/9 - (14/45)uv: which explicit C3/order-4 terms produce it? Is it consistent
    with psi40 in C3? Propose the cheapest further loop/branch on which T5 could NOT be blind (e.g. far edges in Y, infinity, a
    second step across a Y edge), given the separable machinery (Y fixed physical, X continued; or the symmetric X <-> Y swap).
Mark PROVED/VERIFIED/CONJ/OPEN; co-sign or dissent; COLLAB.md entry under PHASE 13; final summary line in your log.
