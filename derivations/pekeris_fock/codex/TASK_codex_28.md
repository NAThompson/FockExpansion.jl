# TASK_codex_28 (phase 13, joint) — 2026-10-24
Write ONLY codex/phase13/x28_*; COLLAB.md / THEORY.md append-only. Serial, < 2 GB, budget <= 1.5 h. No repo edits/commits.
Coordinator asks: evaluate K32 in closed form and derive B_A(u)'s annihilator by creative telescoping (ore_algebra), else fall back.
Claude's feasibility analysis (please confirm or refute, with a concrete check):
 (1) K32(Y;X,eta) = sum_j h_j(Y)h_j(X)h_j(eta)/(n_j h_j(b)) satisfies M_Y K = M_X K = M_eta K with K(b;X,eta) = delta(X-eta):
     for fixed eta it is the rotation-REDUCED Riemann function of the hyperbolic continuation of Lambda^2-32. Reduced over the
     (x3,x4) circle it is a complete elliptic integral of the x18 azimuthal kernel curve (EPD-type with half-integer index),
     not elementary. Unreduced, the variation kernel Delta G(z) of the S^3 resonant pure inverse around z = x.y = 1 is algebraic,
     Q(z)/sqrt(1-z^2), but then B_A is a relative integral over a moving 3-D cap times v0*phi3^A.
 (2) Either way creative telescoping needs a holonomic ideal for phi3^A in 2-3 variables (never computed; its arc-b
     restriction alone needs > 150 unknowns) or a 4-parameter cap integral with log density (x24 chain): beyond reach.
 Tasks: A. Decide whether K32 has a closed form (elementary / hypergeometric / elliptic); give it, or show the elliptic
 reduction explicitly (e.g. check sum_j numerically against your candidate at one (Y,X,eta) using the exact modes
 h0 = sin X, h_pm = 1 +- sqrt2 cos X and Claude's SOps eigenvectors). B. Give Delta G(z) exactly. C. State whether ANY
 creative-telescoping route for B_A is feasible within a few hours, and if so specify it concretely. Co-sign or dissent.
Final summary line in your log; COLLAB entry under PHASE 13.
