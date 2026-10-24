# Enriched physical zero-tail selector

Codex, 2026-10-22, TASK_codex_11. PROVED below means a supplied mathematical argument; finite numerical checks are separately VERIFIED. Assumptions and the plain-space core proof are in ../phase5/goal1_theorem.md. Fix eps>0 and the exact isolated simple positive helium ground state (E0,psi0). Rounded E0 is not covered by the exact-eigenvalue conclusion.

## 1. Concrete finite spaces

Write s=r1+r2, x=(r1-r2)/s, y=r12/s, R²=r1²+r2², w=exp(-eps*s), sigma=2eps*s. The full angular triangle is |x|<=y<=1. Consider the nested spaces U_N spanned by:

* all plain singlet polynomial functions w*p(r1,r2,r12) of total degree <=N;
* w*s^n*x^a*y^b*log(sigma)^p for p=1,2, nonnegative n,a,b, a+b<=n, a even, with any nested finite shell cutoffs exhausting the indicated families;
* w*s^n*x^a*y^b*g for g in {Ltilde,alpha-pi/2,V,U}, n>=2, a+b<=n, a even for Ltilde,V, a odd for alpha-pi/2,U; for V require b>=1 exactly.

For the precise build3 cutoffs use p=1: 3<=n<=N-1; p=2: 4<=n<=N-2; atoms: 2<=n<=Na=N-2. Changing log s to log sigma only mixes lower log sectors. The V products are products of yV with the remaining y^(b-1); this is the dimensionless r12*V/s regularization actually used in build3. Literal r12*V products are admissible too. No bare V is admitted.

Build3 replaces the constant plain function by a fixed Fock seed. Thus its literal finite span need not contain our full V_N. Section 3 proves density even for that literal variant; no silent equality of these two finite spans is used. The finite Fock seed consists of admissible polynomial/log terms. No order-three derivative atom is included in this theorem.

## 2. PROVED: all these products are in the physical form domain

It suffices to prove that each regularized angular atom g is bounded with Cartesian weak gradient <=C/R away from R=0. This is stronger than integrability on selected approach paths. Distances and the Cartesian norms |x1+x2| and |x1-x2| are Lipschitz. Since R<=s<=sqrt(2)R, all their ratios by s have gradients O(1/R).

Ltilde is log((s+r12)/R)+(pi-2)/pi*log(R/s), as in atoms3.jl. Both log arguments stay in positive compact intervals. Consequently Ltilde is bounded with gradient O(1/R). The same statement holds for alpha-pi/2, alpha=2atan2(r2,r1): its distance derivatives have magnitude <=2/R, including one-electron collisions by continuous extension. Exchange makes this atom odd.

Put H=|x1+x2|=sqrt(2R²-r12²). On the entire physical domain,

    beta = 2atan2(H,r12)-pi/2,
    yV = (H/s)*(beta+pi/2).

The angle map has gradient O(1/R), since H²+r12²=2R² and both numerator norms are Lipschitz. This proves the bound for yV, including H=0 (equal opposite vectors), r12=0, nuclear collisions and collinear faces. By contrast V itself has a 1/r12 pole and fails H1 near electron coalescence.

Here is a cancellation proof for U, avoiding unsafe separate logarithms on faces. Fold r1>=r2, so 0<=alpha<=pi/2 and |beta|<=alpha. Set a=alpha-beta, b=alpha+beta, hence a,b>=0, a+b<=pi. From codex/phase2/codex_atoms.md, A=log| (cos beta-cos alpha)/(cos beta+cos alpha) |, B=log| (sin alpha+sin beta)/(sin alpha-sin beta) |, and dU=beta*dA+alpha*dB+pi*dJ. Half-angle identities give

    A=log tan(a/2)+log tan(b/2),
    B=log tan(b/2)-log tan(a/2),
    J=log(1+sec(a/2))+log tan(pi/4-b/4).

Therefore the differential of the actual combined U is

    dU = Fa(a) da + Fb(b) db,
    Fa(z)=-z/sin z + pi*sin(z/2)/(2*cos(z/2)*(1+cos(z/2))),
    Fb(z)= z/sin z - pi/(2*cos(z/2)).

The apparent endpoint poles cancel. Their endpoint limits are

    Fa(0)=-1, Fa(pi)=1-pi/2,
    Fb(0)=1-pi/2, Fb(pi)=-1.

Both functions extend continuously to [0,pi], hence are bounded. Thus U has bounded derivatives with respect to alpha,beta on the closed folded triangle, is bounded there by integrating its differential from any interior value, and has Cartesian gradient O(1/R). The trace U=0 on alpha=pi/2 follows from Fa(pi-z)=Fb(z) and U(pi/2,0)=0. Odd exchange extension is consequently continuous and Lipschitz across r1=r2, with no jump distribution. This is the actual U branch of atoms3.jl, not independent principal branches for its separately singular logs. At a nuclear collision the folded angles both approach zero; the proof covers their joint limit. At electron collision (a,b) approaches (0,pi); at the opposite collinear equal-radius configuration it approaches (pi,0). All endpoints have just been included.

A monomial s^n*x^a*y^b times any bounded angular Lipschitz atom has value O(s^n), gradient O(s^(n-1)). For V use yV and b-1 instead of b. Powers of x,y are bounded with gradients O(1/R); the radial degree n, rather than polynomial divisibility, controls this estimate. For a log sector, value is O(s^n*(1+|log s|^p)) and gradient O(s^(n-1)*(1+|log s|^p)). The six-dimensional radial volume is O(R^5 dR), so the worst n=0 kinetic integral is bounded by a constant times

    integral_0^1 R^3*(1+|log R|^4) dR < infinity.

At infinity w beats every finite power and logarithm; the same derivative bounds apply. On an annulus all angular atoms are Lipschitz in the actual Cartesian variables, so there is no extra face/line singularity or matching jump. Cutting out a ball R<h and letting h->0 shows these local gradients are the global weak gradients: the cutoff error is O(h^(2n+4)*(1+|log h|^4)). Thus every stated product belongs to H1(R6), in the required symmetry subspace. Coulomb potential terms are finite by the usual three-dimensional pair Hardy inequalities and Cauchy--Schwarz. Physical weak forms, without collision cutoffs or omitted boundary distributions, are well defined.

## 3. PROVED: core, including the literal fixed-seed build3 variant

For spaces containing V_N, the phase-5 polynomial core theorem immediately implies form density. For literal build3, all exchange-even positive-degree distance monomials times w occur eventually (n>=1), although the constant is replaced. These positive-degree polynomials alone are a form core in six dimensions.

Proof of the latter strengthening: let F be a continuous functional on invariant H1 annihilating them, extended by rotational/exchange averaging. The H1-holomorphic function F(w*exp(t dot X)) on |Re t|<eps has all positive-degree Taylor coefficients zero: the average of a positive-degree Cartesian homogeneous polynomial is either zero or a positive-degree distance polynomial. It is therefore the constant C=F(w). Fourier inversion gives F(w*phi)=C*phi(0) for every Schwartz phi. Choose smooth compactly supported phi_h(X)=eta(X/h), eta(0)=1. Then ||w*phi_h||H1=O(h²), forcing C=0. Fourier uniqueness and the compact-support argument of the previous core proof now give F=0. This uses the zero H1 capacity of a point in six dimensions, not continuity of a corner value. In particular the fixed-seed replacement does not destroy density, though it changes every finite pencil.

## 4. Physical pencil, normalization, zero tails and convergence

Let q be the physical Coulomb weak form on Q=H1 in the singlet S subspace, c=-5, a=q-c<.,.>. Take the bounded, representation-independent positive head

    ell(u)=<w,u>_physical / <w,w>_physical.

It is L2-bounded and hence form-bounded; ell(psi0)>0. It is not the enriched seed coefficient and is not a triple-corner evaluation. For either family choose b0=s*w/ell(s*w), which is present from N>=1 and has ell(b0)=1. Define W_N=U_N intersect ker ell. Given any finite spanning list, replace each generator phi by phi-ell(phi)*b0 and retain a nested independent subset of the resulting physical functions, quotienting exact dependencies. Together with b0 these supply coordinates. Centering changes only the fixed low shell occupied by b0. No infinite log-sector gauge choice is needed for a finite physical span.

For this basis {b0,b1,...}, ell(bj)=delta_j0. The test pencil is explicitly

    K_ij(z)=q(bi,bj)-z<bi,bj>_physical.

Use the unmultiplied weak Hamiltonian, not coefficient extraction from a denominator-cleared atom PDE. Set all generator coordinates outside cutoff N to zero in every sector, set the head coordinate to 1, and solve

    u_N=b0+v_N, v_N in W_N,
    q(t,u_N)-E0<t,u_N>=0 for every t in W_N.

This is the enriched physical zero-tail backward problem. In these centered coordinates it deletes just the head row. Its matrix is generally dense. Shellwise backward Schur elimination is ordinary finite matrix elimination with fill-in; this theorem does not identify it with the sparse formal differential-sector recurrence.

Let psi=psi0/ell(psi0), E*>E0 a spectral lower bound on psi0-perp, and gamma=E*-E0. Normalize psi0 in L2 for the following constant. The phase-5 gap argument, applying to any bounded ell, gives coercivity of b=q-E0<.,.> on W=ker ell with

    beta=[(E*-E0)/(E*-c)] /
          [1+(E0-c)*||ell||_(a*)²/|ell(psi0)|²] >0.

Indeed v=t*psi0+v_perp in W implies |t|<=||ell||_(a*)*||v_perp||a/|ell(psi0)|. On psi0-perp, b>=((E*-E0)/(E*-c))*a. This proves the displayed bound. b<=a and is nonnegative globally. Hence all finite tail matrices are positive definite and the solutions exist uniquely. Since b0 is fixed, density of union U_N implies density of union W_N in W by subtracting ell(.)b0. Energy orthogonality then proves

    ||u_N-psi||a <= beta^(-1/2)
       * inf_{v in U_N, ell(v)=1} ||v-psi||a ->0.

Writing e_N=u_N-psi, d_N=b(u_N,u_N)=b(e_N,e_N), one also has

    beta*||e_N||a² <= d_N
       <= inf_{v in U_N,ell(v)=1} ||v-psi||a²,
    d_N decreases to zero.

The omitted row g_N=b(b0,u_N) equals d_N, since b(v_N,u_N)=0. Thus the missing source vanishes in the limit and density recovers every physical test, giving the unique normalized physical ground state. For prescribed exterior functions tau_N with ||tau_N||a->0 and ell(tau_N)=0, the induced interior correction satisfies ||delta v_N||a<=beta^-1*||tau_N||a. Uncentered exterior data can first be centered by subtracting ell(tau_N)b0 with the corresponding bounded-factor estimate. No raw coefficient-size tail bound is substituted for this hypothesis.

For Im z!=0, on W the a-representative of b_z is I-(z-c)T, where a(t,Tv)=<t,v>, T positive self-adjoint with spectrum in [0,1/(E0-c)]. The same interval bounds every finite compression. Thus beta(z)=min_mu |1-(z-c)mu|>0 gives uniform finite inverse bounds, affine quasioptimality with constant 1+(1+|z-c|/(E0-c))/beta(z), and form convergence of normalized tails. Their limiting source g(z)*ell is nonzero off the real spectrum. In centered finite coordinates wherever the full K_N(z) is invertible,

    u_N coefficient vector = K_N(z)^(-1)e0 / (K_N(z)^(-1))00,
    g_N(z)=1/(K_N(z)^(-1))00.

The tail equations still make sense at E0 when the full matrix is singular. Physical full resolvents converge for nonreal z, ground Ritz energies decrease to E0, and phase-aligned Ritz ground vectors converge in form norm by the same core/min--max proof. Near E0 the ground Ritz vector has nonzero ell and tail coercivity persists; the nearby zeros of g_N are exactly the ground Ritz energies. Any redundant raw representation must be quotiented before using matrix inverses.

## 5. Rate and scope

PROVED rate transfer: if the affine best-approximation error is <=C*N^(-r), then ||u_N-psi||a<=C*beta^(-1/2)*N^(-r) and 0<=g_N<=C²*N^(-2r). More generally insert any proved approximation function eta(N) in these bounds. The enriched infimum is no worse than the plain one for the version retaining V_N. No explicit r is established: local Fock expansions, finite energy tables and form density alone do not give a global weighted approximation theorem. Literal fixed-seed build3 has the separate core proof above, not an assumed comparison at each finite N.

This closes the enriched *physical weak-pencil* extension requested here, including all order-two products and exact zero tails. It leaves open a faithful infinite enriched coefficient reconstruction/gauge theorem, equality to the cleared formal atom recurrence, classification of every formal solution, an explicit shell exponent, and continuity of a21. The seed's corner coefficient cannot be adopted as ell without a stronger topology.
