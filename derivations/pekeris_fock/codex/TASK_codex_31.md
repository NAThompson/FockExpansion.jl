# TASK_codex_31 (phase 13) — 2026-10-24
Write ONLY codex/phase13/x31_*; COLLAB.md / THEORY.md append-only. Serial, < 1.5 GB, <= 1.5 h. No repo edits/commits.
Claude found (numerically) a CLOSED FORM for the flux kernel K32 (c = 32, b = 3pi/2, I = J = (pi/2, 3pi/2)):
   K32(Y; X, eta) = Theta(z0 - 1) * T3(z0) / (2 pi sqrt(z0^2 - 1)),   T3(z) = 4z^3 - 3z,
   z0 = cos((b+Y)/2) cos((X+eta)/2) + cos((Y-b)/2) cos((eta-X)/2)   (= u_obs u' + v_obs v', observation on the circle).
Evidence: int_I K32(Y;X,eta) h(eta) d eta = h(Y) h(X)/h(b) holds for the exact modes h0 = sin, h_pm = 1 +- sqrt2 cos (ratio 2pi*
(1/2pi) = 1 to 25 digits, claude11/riemann/test_kernel.py) and for four generic numerical c=32 modes (ratio 1 to 1e-6..1e-8,
test_modes.jl). So B(Y) = int_J int_I K32 R dX d eta (R = modal right side, start analytic at X = pi/2) is an explicit
characteristic-region (z0 > 1) integral with an algebraic kernel: the U(1)-reduced Riemann function is NOT an elliptic average here
because the observation point lies on the fixed circle, where the rotation orbit degenerates.
Tasks: (A) PROVE or refute the formula: show K satisfies M_Y K = M_X K = M_eta K on z0 > 1 (exact symbolic check), the correct jump/
boundary behaviour on z0 = 1 (characteristic) and the endpoint normalization K(b;X,eta) = delta(X-eta), and that no extra boundary
terms arise at X = pi/2 (regular start) or from log sources at X = b; state the general-c version (T_{(k-1)}?) if easy.
(B) Given B_P(Y) = int int K32 (w v0/8) P (P = Z^1 panel-1/2 moment part of psi30, x16 23-state module), assess concretely whether
creative telescoping of this (X, eta, q[, atan-representation]) integral is feasible in hours, and propose the smallest exact route.
Mark PROVED/VERIFIED/CONJ/OPEN; co-sign/dissent; final summary line in your log; COLLAB entry under PHASE 13.
