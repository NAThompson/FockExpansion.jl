# TASK_codex_6 — focused answers, 2026-10-22

## 1. Minimal rank: unresolved

No numerical minimal rank is reported. The available number **216 is a chartwise algebraic-field generator bound**, not a cyclic rank, rational rank, or evidence of stabilization. No two-point derivative-order stabilization experiment has been completed. This is the task's explicitly permitted “rank out of reach” outcome.

The concrete missing computation is a combined connection and output row for the complete physical combination (classical functions, the four pulled-back panels, and chi), followed by differential row reduction. If X obeys d_i X=C_i X and f=v X, the next coefficient row is d_i v+v C_i. Iterating only multiplication by C_i is wrong. The supplied executable connection handles each panel separately; it does not implement this combined output-row differentiation, identify cross-panel identities, or quotient the old module. High-order finite differences of actual scalar values would be a badly conditioned replacement. At a single point all scalar derivatives are numbers and do not by themselves have the proposed vector rank.

Even a stable numerical row rank in the padded representation would initially be a rank in that representation: actual-function relations among its generators can lower it. The existing endpoint and duplicated jet relations must first be removed (and cross-panel relations investigated). Numerical witnesses can establish a nonzero minor in a specified formal representation; they do not alone prove actual-function independence or an upper bound at all derivative orders. No rank is inferred from the earlier 65x65 curvature check. Fixing Z,E,a21 does not fix these missing steps.

## 2. Radicals and finite recurrence

Put Q=r1²+r2², P²=2Q-c², c=r12 and

    D=(Q-c²)²-4r1²r2²
     =((r1+r2)²-c²)((r1-r2)²-c²).

A sufficient coefficient field for the present connection is

    K=C(r1,r2,c), F0=K(R=sqrt(Q), P=sqrt(2Q-c²), A=sqrt(D)).

Constants include i, pi, sqrt(2) and the fixed physical constants. This is a sufficient field, **not a proved necessary/minimal radical list**. D is negative in the strict triangle interior; A can be imaginary, with compatible complex branches. The physical functions remain real.

**Exact new result:** for all four maps in structure.md,

    Delta_j²=4D/Q²,  Delta_j=2A/Q (choose signs consistently).

The panel h values are r1/r2, r2/r1, c/P, c/P; hence they require no additional square roots. For each panel adjoin just

    z_j=sqrt((2-4y_j²+Delta_j)/(2(1-b_j))).

Its four quartic roots are ±z_j and ±h_j/z_j. Thus four further quadratic extensions suffice for every pole in the saved connection. The total degree over K is at most 2^7=128, often smaller; this is not an independence assertion. y_j are -c/(sqrt(2)R), -c/(sqrt(2)R), -r1/R, -r2/R. Endpoint square roots also add nothing: sqrt(T(u))=-2y u for y<0 and u=1 or h. Likewise sqrt(1+b) is sqrt(2)r1/R, sqrt(2)r2/R, c/R, c/R.

The classical coefficient field does not require separate quarter-angle roots: cot(w), tan(w), and their pi/4 shifts are rational in exp(2iw), and
exp(i alpha/2)=(r1+i r2)/R, exp(i p0/2)=(P+ic)/(sqrt(2)R).
For w=(alpha-pi/2 ±(p0-pi))/4 these expressions already give exp(2iw). Arguments/angles are atoms, rather than coefficient radicals. The old order-two connection uses rational functions and P/c (its h), so is contained in the same field.

**Finite-factor argument (proved locally):** in characteristic zero, derivatives of algebraic field elements stay in the same finite extension. Restrict scalars using a K-basis of F; each algebraic-coefficient connection becomes a finite rational matrix. Each denominator has a rational norm, so clearing rational denominators gives polynomial coefficient matrices of finite degree. Polynomial multipliers and finitely many distance derivatives become a finite Laguerre stencil after the linear perimetric substitution and the usual additional coordinate row clearing. Therefore these radicals cause a finite factor, not an infinite atom cascade. A loose sufficient dimension bound from this construction is 216*128=27648 over K, expressly not a minimal rank.

Degrees established here: Q and P² have degree 2, D degree 4, and the distance coordinate factors have degree 1. These are **radicand degrees, not complete row-clearing degrees**. Pole collisions z_i-z_j, endpoints u-z_i, root derivatives T_q(z_i), and algebraic-basis denominators must also be included; their rational norms have not been expanded. No total clearing degree or bandwidth is claimed. The precise local denominator prescription remains connection.md. This is enough to prove finiteness locally, not a single globally regular recurrence.

Reproduction: `PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase3/verify_radicals6.py`. Its saved output checks all four discriminant pullbacks and the factorization exactly, and prints h² and angle-pole collision polynomials. In particular those collisions reduce to P² or electron-coordinate factors in the physical maps.

## 3. Obstructions and scope

An elliptic period is not an obstruction to finite differential closure: the supplied explicit finite connection includes these integrals with their endpoint data. Nothing in the displayed formula requires a distance transcendental outside the locally D-finite class. This statement is local and fixed order; it does not assert an all-order Fock alphabet.

The actual unresolved obstruction to using these *particular generators* globally is the moving endpoint u=min(1,h), electron folding, and continuation at discriminants/collision faces. Differentiating a piecewise panel atom repeatedly without interface matching can produce distributional terms. Multiplying an interior equation by a polynomial vanishing on an interface does not supply the missing matching conditions. Nontrivial period monodromy requires branch/initial data, not infinitely many local generators. A global physical atom basis, rational minimality, and a globally valid recurrence/domain remain unproved. No failure of D-finiteness has been found.

## Status

- **Proved by explicit algebra / exact checked identities:** common discriminant; sufficient finite radical tower and degree bound; local restriction-of-scalars/denominator-clearing implication.
- **Previously numerically verified, not rerun for this note:** panel connection and reconstructed psi30 values, as recorded in phase-three logs.
- **Unresolved:** minimal rank, two-point order stabilization, minimal radical necessity, total clearing degrees, and global interface realization. These are not described as conjecturally computed numbers.
