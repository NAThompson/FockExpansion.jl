# TASK_codex_27 (phase 13, joint) — 2026-10-24
Write ONLY codex/phase13/x27_*; COLLAB.md / THEORY.md append-only. Serial, < 2 GB (use --heap-size-hint for Julia), 1-3 h.
Read COLLAB.md PHASE 13 (Claude's ladder entry), claude11/chain.jl (LADDER env), claude11/lad_B1.jls (Chebyshev series in ty of the
far-edge log coefficient of ladder path A, key (0,(4,0),"A"); arc b: X = 3pi/2, alpha = 5pi/4 + pi ty/4, u = cos alpha in (-1,0)).
GOAL: decide whether the far-edge variation B_A(u) of the panel ladder A = (v1, v0, v0, v0) (the cleanest T5 carrier: e-e Coulomb x
panel at step 4) has an H^1(E_q) constituent (your x25 fingerprint), by an EXACT route independent of Claude's brute-force ODE guess:
 1. Per Y-mode the log coefficient is a Sturm flux (your x26 formula): B_j = int_J h_j r_j dX / h_j(Xb). Turn this into a closed
    representation of B_A as an explicit low-dimensional integral (e.g. Riemann function of the hyperbolic operator M_X - M_Y with
    c = 32 in characteristic/tangent-line coordinates, or the S^3 Green kernel picture), with the source r = -w * (-2 v0 g)/16,
    g = N_b phi3^A (Claude: on arc b, g has an exact order-2 operator that factors into first-order pieces; see op_A3.txt / fac_A3.sage).
 2. Compute the minimal ODE of B_A in u (creative telescoping, or a certified-gap guess from your own 100-digit evaluation of the
    integral) and factor it; test for a factor gauge-equivalent to [4u^2(u^2-1)D^2+4u(3u^2-1)D+(4u^2+1)] (or its sqrt(u) twist).
 3. If step 1 stalls, state exactly which kernel/transport is missing. Also audit Claude's statement that the accessible corner orbit
    is {Psi, N_b Psi, H} (with deck involutions) and say precisely what that does and does not exclude (sub/middle/quotient).
Mark PROVED/VERIFIED/CONJ/OPEN; co-sign or dissent; COLLAB.md entry under PHASE 13; final summary line in your log.
