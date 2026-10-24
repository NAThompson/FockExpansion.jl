# Triple-coalescence literature audit — 2026-10-21

Scope: task 2 of TASK_codex_1.md only. COLLAB.md and TASK_codex_1.md were read first; neither sym/ nor julia/ was inspected. Live web access was available. Labels below distinguish **CITED** (primary text actually inspected, with locator), **DERIVED** (our inference/calculation), and **UNVERIFIED** (full text or exact assertion not checked). Links identify the inspected versions; theorem numbering is that of the linked preprint when specified.

## Main answer

**DERIVED from the explicit logarithmic coefficients below:** A representation

    Psi = A + R^2 log(R) B

with A and B analytic in ordinary distances, Cartesian coordinates, or in R and regular angular coordinates, cannot describe normalized singlet helium near the triple point. On a generic fixed ray such analytic functions have ordinary convergent power series in R, so this expression permits log degree at most one. The nonzero R^4 log(R)^2 coefficient contradicts it. This argument also excludes any finite algebraic/analytic change of variables whose restriction to a ray is an ordinary power series. It does not exclude analyticity in an *additional independent variable* t = R^2 log R, or a logarithmic/nonanalytic coordinate change; that is a different and weaker assertion.

**CITED:** Neither FHOS theorem examined below states this additive analytic representation at triple coalescence. Their 2005 factorization has a C^{1,1} remainder; their 2009 analytic theorem concerns an isolated pair collision and explicitly removes multiple collisions.

**Caution:** The existence of a formal Fock hierarchy, convergence of certain Fock-type solutions, and representation of the physical decaying eigenfunction by a globally convergent Fock series are separate propositions. Morgan's own abstract distinguishes them.

## FHOS 2005: actual regularity theorem

**CITED:** S. Fournais, M. Hoffmann-Ostenhof, T. Hoffmann-Ostenhof, T. Ø. Sørensen, *Sharp regularity results for Coulombic many-electron wave functions*, Commun. Math. Phys. **255**, 183–227 (2005), DOI [10.1007/s00220-004-1257-6](https://doi.org/10.1007/s00220-004-1257-6). Inspected [arXiv:math-ph/0312060v1](https://arxiv.org/pdf/math-ph/0312060), Theorem 1.1, Eqs. (1.8)–(1.12):

    Psi = exp(F2 + F3) phi3,     phi3 in C^{1,1}_loc.

In their H = -Delta + V convention,

    F2 = -(1/2) sum_(ell,i) Z_ell |x_i-X_ell|
         +(1/4) sum_(i<j) |x_i-x_j|,
    F3 = ((2-pi)/(12 pi)) sum_(ell,i<j)
         Z_ell ((x_i-X_ell) dot (x_j-X_ell))
         log(|x_i-X_ell|^2+|x_j-X_ell|^2).

The theorem includes optimality among universal factors independent of the individual solution and eigenvalue. It gives no analytic remainder. “Optimal” here must not be interpreted as ruling out a solution-dependent analytic desingularization.

**DERIVED, normalization warning:** Their unit kinetic coefficient and pair interaction coefficient one differ from the present atomic-unit H. Setting x = 2 r transforms their atomic operator into one-half the present operator with the same nuclear Z. Thus locally, in our physical distances,

    F2 = -Z(r1+r2) + r12/2,
    F3 = -Z (pi-2)/(3 pi) (r1_vector dot r2_vector) log(R^2)
         + a smooth homogeneous quadratic term,
       = -Z (pi-2)/(3 pi) (r1^2+r2^2-r12^2) log R
         + a smooth homogeneous quadratic term.

The extra term comes from log(4R^2) and may be absorbed into the remainder. This reconciles the leading coefficient with the Fock convention below; copying FHOS's numerical prefactor without rescaling would be wrong. Exponentiating this universal log factor creates higher log powers, but C^{1,1} regularity of phi3 does not determine their exact coefficients: R^4 log^2 R already has bounded second derivatives.

## FHOS 2009: excludes the triple point

**CITED:** Same authors, *Analytic structure of many-body Coulombic wave functions*, Commun. Math. Phys. **289**, 291–310 (2009), DOI [10.1007/s00220-008-0664-5](https://doi.org/10.1007/s00220-008-0664-5). Inspected [arXiv:0806.1004v1](https://arxiv.org/pdf/0806.1004), Theorem 1.4 and definitions (1.7)–(1.10). Equations (1.12)–(1.13) give

    Psi = Psi_k^(1) + |x_k| Psi_k^(2),
    Psi = Psi_kl^(1) + |x_k-x_l| Psi_kl^(2),

with real analytic coefficient functions in Cartesian variables, locally near the corresponding isolated two-particle collision sets. Other collisions are removed in defining these sets. At x1=x2=0 the hypotheses do not hold. Theorem 1.1 is the corresponding one-particle Coulomb statement, not a triple-collision logarithmic theorem. The proof uses the Kustaanheimo–Stiefel transform and analytic hypoellipticity (Remark 1.5).

## Fock 1954 and explicit logarithmic orders

**UNVERIFIED original full text:** V. A. Fock, Izv. Akad. Nauk SSSR, Ser. Fiz. **18**, 161–172 (1954), commonly titled *On the Schrödinger equation of the helium atom*. English translation in K. Norske Vidensk. Selsk. Forh. **31**, 138–152 (1958). Original equation numbers and proof hypotheses were not verified here. Do not attribute a precise convergence theorem to that paper on this audit alone.

**CITED independent primary calculation:** E. Z. Liverts and N. Barnea, [arXiv:1505.02351v1](https://arxiv.org/pdf/1505.02351), published after revision as Phys. Rev. A **92**, 042512 (2015). Eq. (1) states

    Psi = sum_(n>=0) R^n sum_(0<=p<=floor(n/2)) psi_np(alpha,theta) (log R)^p.

Eqs. (17) and (19), for psi00=1, give

    psi21 = -Z (pi-2)/(3 pi) sin(alpha) cos(theta),
    psi42 = Z^2 (pi-2)(5 pi-14)/(180 pi^2)
            [1 - 2 sin(alpha)^2 sin(theta)^2].

Their alpha=2 atan(r2/r1); cos(theta)=(r1^2+r2^2-r12^2)/(2 r1 r2). In particular psi42 is not identically zero. These are formulas checked in the linked preprint, not recalled equation numbers from the final journal version.

**DERIVED:** The allowed log powers are:

| radial degree n | log degrees allowed by the Fock hierarchy |
|---|---|
| 0,1 | 0 |
| 2,3 | 0,1 |
| 4,5 | 0,1,2 |
| 6,7 | 0,1,2,3 |
| general n | 0 through floor(n/2) |

This is an allowed-order statement, not an assertion that every coefficient at every angle is nonzero. First log: R^2 log R; first squared log: R^4 log^2 R. The bound implies a formal regrouped sector F_p begins at homogeneous degree 2p. “F_p=O(R^{2p})” requires a suitable bounded angular coefficient class and actual summability; it is not by itself an analyticity theorem.

**DERIVED:** Because s=r1+r2=R q(Omega), with 1<=q<=sqrt(2),

    log R = log(s) - log q(Omega).

Changing from log R to log s (or log sigma = log(2 eps s)) mixes coefficients triangularly at each fixed homogeneous order. It preserves the maximum logarithmic degree and the first possible radial order, while generally inserting non-polynomial angular functions into the coefficients. A log(s) basis is therefore structurally plausible but does not make its log-free sectors analytic in perimetric distances.

## Morgan 1986: strong convergence result, important qualification

**CITED, publisher abstract only:** J. D. Morgan III, *Convergence properties of Fock's expansion for S-state eigenfunctions of the helium atom*, Theor. Chim. Acta **69**, 181–223 (1986), [DOI 10.1007/BF00526420](https://link.springer.com/article/10.1007/BF00526420), abstract. It describes pointwise convergence for all hyperradii of Fock expansions constructed using functional analysis, improving Macek's mean-convergence result. Its radial symbol is the squared physical hyperradius: its generators R^{1/2} and R log R translate into our R and R^2 log R (up to a constant). The abstract also states that infinitely many Fock-type solutions exist for every complex energy without a boundary condition at infinity. It explicitly leaves open whether the physical exponentially decaying eigenfunction at a point-spectrum energy is representable by that expansion.

**UNVERIFIED:** Full text was behind subscription access; theorem numbers, angular function spaces, and the exact class of admissible free angular data were not checked. Consequently this report does **not** certify either a globally convergent physical-ground-state Fock expansion or joint analyticity in (R,R^2 log R). The abstract's wording cannot establish those stronger propositions. Some later papers cite Morgan as proving convergence without repeating this qualification; that abbreviated citation is inadequate for our global ansatz proof.

## Abbott–Gottschalk–Maslen 1987: three papers, not one

**CITED publisher abstract:** P. C. Abbott and E. N. Maslen, *Coordinate systems and analytic expansions for three-body atomic wavefunctions. I. Partial summation for the Fock expansion in hyperspherical coordinates*, J. Phys. A **20**, 2043–2075 (1987), [DOI 10.1088/0305-4470/20/8/023](https://doi.org/10.1088/0305-4470/20/8/023). It develops the infinite coupled angular differential equations and compares partial summation approaches of Pluvinage and Hylleraas. No all-order analytic two-term triple-point theorem was verified.

**CITED publisher abstract:** J. E. Gottschalk, P. C. Abbott and E. N. Maslen, *II. Closed form wavefunction to second order in r*, J. Phys. A **20**, 2077–2104 (1987), [DOI 10.1088/0305-4470/20/8/024](https://doi.org/10.1088/0305-4470/20/8/024). It obtains a closed form through second hyperradial order for singlet S states and discusses relations among coordinate systems. This is an explicit finite-order result, not a closed-form solution to all orders. (The institutional database lists final page 2103; the traditional bibliography lists 2104.)

**UNVERIFIED full text:** J. E. Gottschalk and E. N. Maslen, *III. Derivative continuity via solutions to Laplace's equation*, J. Phys. A **20**, 2781–2803 (1987), [DOI 10.1088/0305-4470/20/10/022](https://doi.org/10.1088/0305-4470/20/10/022). Exact derivative-continuity results and equation numbers not audited. No conclusion in this report relies on them. Part I scan was found in search but failed full opening; Parts II/III fulltexts were not retrieved.

## Frankowski–Pekeris 1966: what was actually used

**CITED metadata:** K. Frankowski and C. L. Pekeris, *Logarithmic terms in the wave functions of the ground state of two-electron atoms*, **Physical Review** **146**, 46–49 (1966), [DOI 10.1103/PhysRev.146.46](https://journals.aps.org/pr/abstract/10.1103/PhysRev.146.46); erratum **150**, 366 (1966), [DOI 10.1103/PhysRev.150.366](https://journals.aps.org/pr/abstract/10.1103/PhysRev.150.366). This predates the separate journal Physical Review A: “PRA 146,46” is bibliographically incorrect. APS full text required credentials; the original basis restrictions and the erratum's content remain unverified.

**CITED later primary work, exact locator:** V. I. Korobov, *Coulomb three-body bound-state problem: Variational calculations of nonrelativistic energies*, Phys. Rev. A **61**, 064503 (2000), [author-hosted PDF](https://theor.jinr.ru/~korobov/papers/variational_PRA.pdf), Eq. (2), attributes the compact representation

    exp(-k s) sum c_lmnij s^l u^m t^(2n) (s^2+t^2)^(i/2) (log s)^j

(up to dummy index names) to Frankowski–Pekeris, with s=r1+r2, t=r1-r2 and u=r12. Thus powers of log(s), together with additional radial powers, were in use; it was not simply an analytic A plus one R^2 log R term. This author is a primary source for his printed comparison, but this is **indirect evidence** of the 1966 basis, not inspection of the original equation. Exact admitted values/maxima of i,j and whether every factor was already present in the original versus later refinements must remain UNVERIFIED. Do not silently infer that their finite variational basis incorporated the entire exact infinite log tower.

## Consequences for the proposed program

**DERIVED/recommendation:** Use the polyhomogeneous Fock index set n>=2p as a local admissibility constraint, and specify a log-free homogeneous/angular coefficient class to give the sector decomposition meaning. Do not impose ordinary analyticity of each perimetric sector based on either FHOS paper. Keep global completeness/convergence and boundary/minimality selection as new obligations. Local equations leave free harmonic data; the literature audited here does not turn these into a unique physical solution without global information.

No files outside codex/ were modified by this literature workstream. No paper authors were contacted, subscriptions bypassed, or unverified fulltexts represented as read.
