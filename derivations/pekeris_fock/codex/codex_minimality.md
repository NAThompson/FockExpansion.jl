# Shell recurrence, physical admissibility, and Fock functionals

Codex, 2026-10-21. Normalization is P=-Q exp(sigma/2)(H-E)exp(-sigma/2), Q=(u+v)(2u+w)(2v+w), sigma=u+v+w. Fix eps>0 before varying E in a spectral pencil. Statements marked PROVED below include their argument; NUMERICAL means finite computation only; CONJECTURED includes missing analytic hypotheses rather than an asserted theorem. Literature citations and their access limits are in codex_literature.md. No growing-block Pincherle theorem is invoked.

## 1. Exact algebra and dimensions

**PROVED (counting).** Shell N consists of l+m+n=N, dimension d_N=(N+1)(N+2)/2. Exchange u<->v has f_N=floor(N/2)+1 fixed triples. Thus the singlet and triplet dimensions are d_N^+=(d_N+f_N)/2=floor((N+2)^2/4), and d_N^-=(d_N-f_N)/2. A coefficient attached to l<m occurs twice in the full expansion: its contribution to the unweighted coefficient norm is 2|A_lmn|^2. The normalized symmetric basis has coefficient sqrt(2) A_lmn there.

**PROVED algebraically; independently machine-verified.** For k=0,

    sum_h B_N,h(E,eps) a_(N+h)=0,

where B_N,h is d_N by d_(N+h) (replace dimensions by d^+ in the singlet sector). Negative-index shells are absent. P has 57 offsets and h=-3,...,3 at fixed eps and general E; after E=-eps^2 it has 33 offsets and h=-2,...,2. The saved independent symbolic derivation and 441 cross-comparisons are in codex_crosscheck.md. Bandwidth means allowed coupling range, not invertibility of any outer block.

**PROVED (product rule plus the same polynomial expansion).** Put A_N=(a_N^(0),...,a_N^(K)), with K the largest retained log power, hence K+1 sectors and dimension (K+1)d_N^+. Then

    (sigma^2 P) F_j +(j+1) C1 F_(j+1)+(j+1)(j+2) C2 F_(j+2)=0.

The sector matrix is upper triangular with two upper diagonals. Shell half-bandwidths are 5 for sigma^2 P, 4 for C1, 3 for C2 at general E; 4,4,3 at E=-eps^2. Offset counts respectively are 225/177,117/117,57/57 (general/fixed). Consequently the whole finite-sector recurrence has half-bandwidth 5 or 4. Grouping b successive shells for half-bandwidth b makes the system block tridiagonal, but block sizes still grow and off-diagonal blocks are rectangular. For the exact infinite tower each shell has infinitely many sector coordinates. F_j=O(R^(2j)) constrains local radial jets, not Laguerre shell indices. It does not justify j<=N/2.

**PROVED qualification.** An infinite upper-triangular hierarchy has no highest sector; 'solve from the top down' is not an algorithm for it. Finite-K closure needs an approximation prescription and error analysis. If exp(-sigma/2)F_K is itself in the physical operator domain at the simple ground eigenvalue, PF_K=0 implies it is a ground eigenfunction; its positive-order corner vanishing excludes it. Induction excludes all positive sectors under that extra hypothesis. The sector PDE alone gives no such implication.

## 2. Admissible-subspace formulation

**PROVED (physical form construction).** On full six-dimensional configuration space the kinetic form has domain H^1. The three-dimensional Hardy inequality, applied in each electron coordinate or the relative electron coordinate, gives ||u/r||_2<=C||gradient u||_2. Cauchy-Schwarz and Young then bound each Coulomb form by delta||gradient u||_2^2+C_delta||u||_2^2 for every delta>0. Adding a sufficiently large norm shift makes the Coulomb form norm equivalent to H^1; completeness proves closedness and a lower bound. Restriction to the closed singlet rotationally invariant subspace preserves this construction. Its associated self-adjoint operator is what H means below; simplicity of the specified ground eigenvalue is an explicit spectral input.

**PROVED (definition and conditional equivalence).** Let H be the physical self-adjoint Coulomb Hamiltonian in the singlet S subspace of L^2(dmu), with its physical closed quadratic form q. Define the admissible solution subspace at E to be

    { A : the reconstructed Psi is in D(q),
          q(phi,Psi)=E <phi,Psi> for every phi in D(q) }.

For a tower also impose a specified log-free homogeneous corner class, summability topology, and normalization of any flat-sector ambiguity. Only the sum Psi must lie in D(q); individual sectors need not do so. The representation theorem for a closed form identifies these weak solutions with ker(H-E). Indeed the displayed identity is precisely the definition of the associated operator with H Psi=E Psi. At a simple ground eigenvalue the image of this coefficient subspace in physical functions is one dimensional. Uniqueness of coefficients additionally needs injective reconstruction. This is an admissibility condition, not a scalar ratio-minimality theorem. 'Minimality' may be used as shorthand only with this definition stated.

**CONJECTURED / analytic work still required.** All formal recurrence solutions satisfying chosen coefficient bounds reconstruct physical weak solutions, and all physical solutions admit the globally convergent, gauge-fixed sector representation. The polynomial multiplier Q vanishes on boundary strata, so an interior identity is not a substitute for the domain and boundary conditions. Global convergence of the log tower is not supplied by local Fock algebra. Ground-state simplicity and the physical Coulomb form realization are standard spectral inputs here, not consequences of the recurrence.

**PROVED (Parseval and coordinate change).** For an unweighted Laguerre L^2 expansion,

    sum |A_lmn|^2 = integral exp(-sigma)|F|^2 du dv dw,
    ||Psi||_phys^2 = pi^2/(32 eps^6) integral exp(-sigma) Q |F|^2 du dv dw.

Here dmu=8 pi^2 r1 r2 r12 dr1 dr2 dr12, |Jac|=1/(8eps^3), r1 r2 r12=Q/(32eps^3). These norms are not equivalent. On an angular interior cutoff, F~rho^-2 near the origin is physical-L^2 but not unweighted-L^2; at infinity F~exp(rho/2)rho^-2 has the opposite property. These examples concern norms, not eigenfunctions. Since Q<=C sigma^3 and its Laguerre matrix has bandwidth 3 and entries O((1+N)^3), a weighted Schur bound on finite sequences gives the sufficient condition sum (1+N)^3|A|^2<infinity for physical norm convergence by completion. It does not by itself ensure finite kinetic energy or the weak equation. Log-sector cross terms require a joint norm estimate.

## 3. Schur complements: the rigorous finite statement and a conditional limit

**PROVED (finite matrices).** At fixed eps, let H_L and G_L be the physical energy and Gram matrices of a finite linearly independent polynomial basis. Apart from a common positive constant they are -<L,P(E=0)L> and <L,Q L>. They are Hermitian and G_L>0 when integrals and boundary forms are valid. At Im z !=0, H_L-zG_L and every principal tail are invertible: a null vector would give Im(z) x*G_L x=0. After grouping three polynomial shells, both matrices are block tridiagonal. Elimination from the last block gives

    S_i(z)=D_i(z)-U_i(z) S_(i+1)(z)^(-1) V_i(z),
    D_i=H_ii-zG_ii, U_i=H_i,i+1-zG_i,i+1,
    V_i=H_i+1,i-zG_i+1,i.

All dimensions conform even when blocks grow. V_i(z) is not U_i(z)* for nonreal z: the conjugated spectral parameter would differ. This is a finite matrix continued fraction, not a forward transfer system. In log-augmented trial bases log moments generally destroy polynomial shell bandedness; use an ordinary finite Schur complement rather than assuming the formal sector bandwidth is the bandwidth of a physical Hermitian matrix.

**PROVED under explicit core hypothesis.** Suppose nested finite trial spaces V_L have union dense in D(q) in form norm. Choose a real c below the form lower bound. The shifted variational problem is coercive; orthogonal projection in the q-c inner product proves Galerkin convergence of (H_L-cG_L)^(-1) loads to (H-c)^(-1) in form norm. Finite resolvents, embedded using physical orthogonal projections, have norm <=1/|Im z|. The resolvent identity, or bounded analytic continuation of their scalar matrix elements from the subthreshold resolvent set, then gives strong resolvent convergence for nonreal z. Thus physically tested resolvent matrix elements and effective Schur data defined through those elements converge. At a real isolated eigenvalue the resolvent has a pole: recover the eigenprojection by a contour integral/residue, not by assuming tail inverses exist on the real axis.

**CONJECTURED for an unqualified raw-shell limit.** Establishing that these particular polynomial spaces form a Coulomb form core is a missing proof in this collaboration; L^2 Laguerre completeness alone is insufficient. Nor does strong physical resolvent convergence automatically imply convergence of every raw coefficient Schur matrix. Raw Laguerre coefficient extraction is not obviously continuous in the physical norm (the constant dual involves 1/Q near the triple corner). Physical orthogonalization gives a sound head/tail splitting but may remove the simple band pattern. These distinctions prevent a claimed infinite matrix continued fraction theorem for the raw recurrence. Rectangular shell couplings cannot simply be inverted to obtain a fixed-size transfer map.

## 4. Cheap illustration; what it does not establish

**NUMERICAL.** codex_shells.jl includes the read-only pencil builder and computes Float64 generalized eigenvectors at shell cutoffs 12,16,20. It uses multiplicity-two weights off the exchange diagonal and normalizes total coefficient norm to one. See codex_shells.out. At cutoff 20 the shell norms for N=4,8,12,16,20 are 2.0348e-3,1.9310e-4,4.2994e-5,1.0811e-5,8.0556e-7. Cutoff sensitivity is substantial near the artificial boundary: shell 12 is 1.0011e-5 at cutoff 12 versus 4.2994e-5 at cutoff 20. Small finite-pencil residuals certify the finite solve only. These data do not distinguish an asymptotic power law from geometric decay, and do not establish a physical/nonphysical dichotomy. No unique wrong-eps forward solution exists without choosing additional free boundary data; rectangular blocks make that comparison ambiguous.

**PROVED 1D counterexample to 'physical iff algebraic'.** For beta>-1/2, the Laguerre coefficients of exp((1/2-beta)x) are

    a_n=(beta+1/2)^(-1) [(beta-1/2)/(beta+1/2)]^n.

Proof: integrate the Laguerre generating function against exp(-(beta+1/2)x), or integrate each polynomial term. The physical function exp(-beta x) has geometrically decreasing coefficients for beta>0 (beta=1/2 even has finite support); for -1/2<beta<0 it grows and the coefficients grow geometrically. For beta<=-1/2 even these projection integrals fail to exist. Corner logarithms can lead to algebraic tails, but neither decay class is a universal physicality criterion. No 3D nonphysical growth theorem is established.

## 5. Fock constants as functionals

**PROVED sufficient summability.** Corner evaluation F_0(0)=sum A^(0) requires absolute summability or another justified convergence theorem; bare l^2 does not suffice. For a derivative multi-index q, L_n^(q)(0)=(-1)^q binomial(n,q), so absolute weighted sums sum |A| product_i (1+n_i)^(q_i)<infinity permit termwise corner derivatives. By Cauchy-Schwarz, sum (1+N)^p |A|^2<infinity with p>2|q|+3 is sufficient (there are O(N^2) triples). These are sufficient, not necessary, conditions. They are often too restrictive for the actual angular singular structure. Abel damping t^(l+m+n) regularizes formal sums for t<1, but a t->1 limit must be justified separately.

**PROVED conditional local functional.** With psi00=1, physical R=sqrt(r1^2+r2^2), v=sin(alpha)cos(theta), and normalized angular projection P_v[v]=1, a21 is the R^2 log-free projected coefficient after subtracting lower orders and b R^2 log R, b=-Z lambda(pi-2)/(3pi). Conversion from log sigma includes log(sigma/R), and the exponential contributes lower-order products. If the polyhomogeneous remainder is uniformly o(R^2) in projection norm, this finite-part limit converges. Ordinary Taylor jets of non-polynomial angular terms are not a replacement.

**PROVED conditional global identity.** Let f=P_v[Psi], C=RV and h=2P_v[C Psi]/R-2E f for an exact regular decaying solution. The polynomial 2 x_1 dot x_2 is harmonic homogeneous of degree two in six dimensions, so v=2 x_1 dot x_2/R^2 has angular eigenvalue 2(2+4)=12. The projected radial ODE is therefore f''+5f'/R-12f/R^2=h. Its fundamental powers R^2,R^-6 have Wronskian -8R^-5. Variation of constants gives

    f=-1/8 [R^-6 integral_0^R t^7 h(t)dt + R^2 integral_R^infinity h(t)/t dt].

For f=A R+R^2(b log R+a21)+..., h=-7A/R+8b+..., expansion gives

    a21=-1/8 FP integral_0^infinity h(R)/R dR-b/8.

Subtract (-7A)/delta-8b log(delta) from the cutoff integral. The -b/8 comes from the first Green integral. Decay controls the far end; the corner requires this finite part, not an ordinary integral. For trial states this is an estimator, not an identity for their local coefficient. Its defect is a weighted residual finite part; small energy error does not bound it. codex_a21_functional.md supplies the independent Mellin implementation and manufactured check.
