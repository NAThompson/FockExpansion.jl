# TASK_codex_25 (phase 13, joint) — 2026-10-24
Read COLLAB.md PHASE 11-13, THEORY.md Phases 10-11, your x20-x24 reports. Write ONLY codex/phase13/x25_*; COLLAB.md / THEORY.md
append-only. Serial, < 2 GB, budget 1-3 hours. No repo edits, no commits.

## Context: Claude's new route (claude11/), please co-sign or dissent
User instruction: "keep pushing until a concrete obstruction or an explicit formula". Claude proposes a 1-D fingerprint for T5 that avoids
the full relative pairing: the LIMIT at the collinear circle w^2 = 1-u^2-v^2 = 0 (u = x1, v = x2).
 (1) Exact: for S_e2, Y^2 = X(X^2+2BX+AC), B^2 - AC = -P1^2 w^2 (P1 = q^4-1). So on w = 0 every elliptic fibre is nodal: Y^2 = X(X+B)^2.
     The vanishing circle (around the roots X = -B +- i P1 w) swept over 1-cycles in q gives periods -> 2 pi i \int dq/sqrt(-B(q)),
     -B = uv q^4 - 2(1-v^2) q^2 - uv. With r = q^2 and u^2+v^2 = 1 the branch points are r = 0, inf, (u+1)/v, (u-1)/v, i.e.
     E_lim = Legendre curve with lambda = (u-1)/(u+1) = -tan^2(phi/2) (u = cos phi). Claude's claim: T5 has a type II degeneration
     at the circle (exponents {0,0,1,1,2} on x2 = 1/3: N^2 = 0, rank N = 2), Gr^W_1 = H^1(E_lim), and Im N on T5 restricted to the
     circle is spanned by the two E_lim periods.
 (2) C3's pure constituents (Tate, E_old for panels 1/2: eta^2 = (1-u)s^2 - 2vs + (1+u) times s, and the analogous panel 3/4 curve
     (1+v)s^2 + 2us + (1-v)) degenerate MULTIPLICATIVELY on the circle (two roots collide), so for any C3-class function the variation
     around the circle, restricted to the circle, is mixed Tate (no elliptic constituent).
 (3) Hence a fingerprint: if some branch b of psi40 has, at some arc of the (complexified) circle, a log coefficient N b whose restriction
     to the circle has an H^1(E_lim) constituent, then psi40 has a constituent outside C3 (Conjecture A, positive).
 Claude is computing such log coefficients numerically: in (alpha,beta) with u = cos alpha, v = cos beta all three Coulomb radicals
 are resolved; the physical branch is analytic on the open square |tx|,|ty| < 3 (tx = 2(alpha-beta)/pi, ty = 2(alpha+beta-pi)/pi;
 argument: analytic at the four D0 vertices + sliding meridians along lines with only normal crossings), and the spectral Chebyshev
 decay rate 10^-0.773/degree ~ (3+2 sqrt 2)^-1 shows singularities on the far grid lines tx, ty = +-3 (non-physical sheets of the
 circle arcs), already for psi20/psi30. The test is whether the log coefficient of psi40 (C1/C2 sectors) there contains E_lim periods.

## Your workstream (independent of Claude's numerics)
 A. Prove or refute (1): limiting MHS of T5 at the circle; Jordan type; identify Im N with H^1(E_lim) periods; check against Claude's
    exact operator claude10/op_e2_Q.sobj on x2 = 1/3 at u0 = sqrt(8/9) (analytic local solutions vs 2 pi i \int dq/sqrt(-B) over the
    two E_lim cycles, numerically to 30+ digits). Watch the lambda formula and any extra ramification (sqrt factors, which arc).
 B. Prove or refute (2) for the FULL C3 list (order-2 rank-5 module, the 18 classical atoms, the 4x23 panel atoms, and the
    algebraic/radical factors): list all pure constituents and their degenerations on the circle. If any C3 constituent can produce an
    elliptic Gr at the circle (in particular one isogenous to E_lim), say so: it would kill the fingerprint.
 C. Prove or refute (3) and its converse: if T5 is a constituent of D.psi40 in a quotient or middle position, must some branch of psi40
    have N b with an E_lim constituent at some arc (which arc / which sheet)? State exactly what is and is not implied.
 D. Audit Claude's "analytic on the open square |tx|,|ty| < 3" argument (claude11/NOTES.md once written) and later any numerical
    certificate Claude writes in claude11/ (claim files will be named claude11/claim_*.md).
Mark PROVED/VERIFIED/CITED/CONJ/OPEN; co-sign/dissent. Finish with a one-paragraph summary line in your log and a COLLAB.md entry under
PHASE 13.
