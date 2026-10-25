# TASK_codex_29 (phase 13, joint) — 2026-10-24
Write ONLY codex/phase13/x29_*; COLLAB.md / THEORY.md append-only. Serial, < 1.5 GB (a 2.5 GB Claude job is running), <= 2.5 h.
Claude is running ladder A at n=240/704 bits (claude11/ladA_hi_240.log, checkpoints in claude11/ck/) to guess the annihilator of
B_A(u) (far-edge variation on arc b) and test it for the E_q factor E = 4u^2(u^2-1)D^2+4u(3u^2-1)D+(4u^2+1) up to isogeny, the
sqrt(u) twist and gauge. Do TWO things, in this priority order:
 1. PIPELINE CALIBRATION (most useful): build a known test function on arc b (u in (-1,0)) of comparable shape that MUST
    contain E_q, e.g. f = R1(u) p1(u) + R2(u) p1'(u) + (a Tate part: products of logs/dilogs of the radicals sqrt(1+-u),
    sqrt(1-u^2), u with rational coefficients, several terms) where p1 is an E_q period (your x25 normalization, any branch
    analytic on the arc), at 180+ digits. Run Claude's guesser (claude11/guess_hi.jl logic: derivatives at many points, QR gap)
    or your own equivalent; report the smallest (r,d) at which the true operator is detected and whether ore_algebra's
    factor() / a gauge-equivalence test recovers E (or its isogenous/twisted form). Give a reusable script x29_ftest.sage
    that, given an exact operator L (rational coefficients in u), decides whether some right/left factor or LCLM component is
    gauge-equivalent to E or to its sqrt(u)-twist [4(u^2-1)D^2+8uD+1] (e.g. compare local exponents and monodromy/ symmetric
    square invariants, or use ore_algebra's factor plus a Hom/gauge solve). Validate it on the test function's operator.
 2. If time remains: attempt to prove or refute that the rotation average of Delta G (x28) equals K32.
Mark PROVED/VERIFIED/CONJ/OPEN; co-sign/dissent; final summary line in your log; COLLAB entry under PHASE 13.
Sage: /Applications/SageMath-10-3.app/Contents/Frameworks/Sage.framework/Versions/10.3/venv/bin/sage (ore_algebra installed).
