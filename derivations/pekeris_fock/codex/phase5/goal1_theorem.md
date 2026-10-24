# Goal 1: a form-core and a convergent normalized backward criterion

Codex, 2026-10-22, TASK_codex_9. Labels PROVED mean the argument is supplied here, not a claim of external peer review. CITED inputs/access limits are in goal1_sources.md. This concerns fixed eps>0, Z=2, lambda=1 and the singlet rotationally invariant physical Coulomb form. No writes to Claude code. All limits below use exact arithmetic/forms; rounded E0 and quadrature add separate errors.

## 1. Setup and the result

Let X=(x1,x2) in R^6, r_i=|x_i|, t=|x1-x2|, s=r1+r2, w=exp(-eps s)=exp(-sigma/2). Let V_N be w times all exchange-even polynomials in (u,v,w_perimetric) of total degree <=N, equivalently in (r1,r2,t). These are the finite plain Laguerre spaces. Let q be the physical form of H on the singlet S Hilbert space Hcal. Its domain Qcal is the invariant part of H^1(R^6). Hardy plus Young bounds make a=q-c<.,.> an equivalent H^1 inner product for c<-4; c=-5 suffices because the repulsive electron interaction is nonnegative and the two hydrogen Hamiltonians give H>=-Z^2=-4. Write ||.||a for this norm.

CITED spectral input: the physical helium ground eigenvalue E0 is isolated and simple with a strictly positive eigenfunction psi0. Let E*>E0 bound the spectrum on psi0-perp from below, gamma=E*-E0>0. No numerical value of gamma is needed for convergence.

PROVED below:

1. Union V_N is a form core. This removes the earlier explicit density hypothesis for the *plain fixed-eps* spaces.
2. The formal constant Laguerre coefficient ell=A000 extends to a bounded functional on Qcal, although it is not bounded on physical L2. Also ell(psi0)>0. Thus normalize psi=psi0/ell(psi0).
3. At the exact E0 the following computable finite criterion converges to psi in form norm: fix A000=1, set all coefficients above shell N to zero, and solve all recurrence rows through N except the constant row. The principal tail matrix is nonsingular. The dropped constant-row residual tends to zero. This is a genuine normalized zero-tail/backward finite-section criterion for the ground state, without inverting any rectangular outer block.
4. At nonreal z the same normalized tail problems converge in form norm to a unique *inhomogeneous* constrained solution. Merely converging at such z does not imply a homogeneous eigenfunction: the dropped row must also vanish.
5. Nested form-admissible enrichments retaining every V_N inherit core density, resolvent convergence and ground Ritz convergence. The identical raw single-row Miller statement for their gauged/dense coefficient representations requires a separately specified continuous normalization and a conforming physical pencil. Formal log-sector recurrences alone do not supply it.

## 2. PROVED: direct form-core density

Every polynomial-distance times w belongs to H^1: distances are Lipschitz, their a.e. gradients are bounded, polynomial growth is beaten by w. Rotational/exchange averaging P is a bounded projection in H^1.

Suppose a continuous linear functional F on invariant H^1 annihilates union V_N. Extend it to H^1 by F(P .). For any Cartesian polynomial p(X), Pp is a polynomial in the Gram invariants |x1|^2, |x2|^2, x1 dot x2, symmetrized under electron exchange. One way to see this is to average monomial tensors over rotations: invariant contractions pair indices with Kronecker tensors; an epsilon contraction vanishes with only two vector arguments. Hence w Pp is in some V_N, since x1 dot x2=(r1^2+r2^2-t^2)/2. Therefore F(wp)=0 for every Cartesian polynomial p.

The map t -> w exp(t dot X) is H^1-valued holomorphic on the connected tube |Re t|<eps (Euclidean norm in R^6), because s>=|X| and its parameter derivatives have polynomial growth with a uniform exponential margin on compact subtubes. Every Taylor derivative of F(w exp(t dot X)) at zero vanishes. The holomorphic identity theorem therefore gives F(w exp(i k dot X))=0 for every real k.

For a Schwartz function phi, Fourier inversion converges in H^1 after multiplication by w: ||w exp(i k dot X)||H1<=C(1+|k|), and the Fourier transform is rapidly decreasing. Consequently F(w phi)=0. For eta smooth and compactly supported, eta/w is compactly supported H^1 (w^-1 is locally Lipschitz); approximate eta/w in H^1 by smooth compactly supported phi_j. Multiplication by w is bounded H^1->H^1, so F(eta)=lim F(w phi_j)=0. Compactly supported smooth functions are H^1-dense; F=0. Hahn--Banach proves density, and equivalence of H^1 and a norms proves the form-core claim. This argument is independent of the unavailable full theorem hypotheses in Klahn--Bingel.

## 3. PROVED: continuous coefficient head and its kernel

Write Qp=(u+v)(2u+w_perimetric)(2v+w_perimetric), Ceps=pi^2/(32 eps^6), so dmu=Ceps Qp du dv dw_perimetric. For finite Laguerre expansions,

    ell_j(f)=integral f exp(-sigma/2) L_j du dv dw_perimetric
            =Ceps^-1 integral f exp(-sigma/2) L_j/Qp dmu.

For every fixed triple j, Cauchy--Schwarz and the six-dimensional radial Hardy inequality give

    |ell_j(f)| <= Ceps^-1 ||f/R||L2
                 * ||R exp(-sigma/2) L_j/Qp||L2 <= C_j ||f||a,

where R=|X|. The second norm is finite: its square is Ceps integral R^2 exp(-sigma) |L_j|^2/Qp du dv dw. Near the triple corner the radial power is rho^(2-3+2)=rho^1. On each nonzero axis only one factor of Qp vanishes, linearly in the two transverse nonnegative coordinates, whose 1/(a+b) singularity is integrable. The compact angular integral is therefore finite; infinity is exponentially controlled. The formula extends uniquely to Qcal. It need not define a convergent bare-l2 expansion for arbitrary form functions. In the symmetric basis use ell_000 unchanged and multiplicity two for off-exchange-diagonal raw coefficients.

For ell=ell_000 the representing density is positive; positivity of psi0 gives ell(psi0)>0. Let b0=w; ell(b0)=1. Let W=ker ell and W_N=V_N intersect W. These are precisely finite spans with the constant Laguerre coefficient zero. W_N are form-dense in W: given v_N->v in form norm, replace v_N by v_N-ell(v_N)b0. Likewise the affine spaces b0+W_N are form-dense in {ell=1}. Fixed finite heads are therefore continuous under the convergence proved here, improving the earlier L2-only caveat. This does not establish convergence of head *inverses* at poles or convergence of unbounded shell-index norms.

## 4. PROVED: coercivity of the real-eigenvalue tail

Set b_E(u,v)=q(u,v)-E<u,v>. On W, b_E0 is coercive in ||.||a, although it is only semidefinite on the full form space. Indeed decompose v=t psi0+v_perp with ||psi0||L2=1. The constraint gives |t|<=||ell||a* ||v_perp||a/|ell(psi0)|. The spectral gap gives

    b_E0(v,v)=b_E0(v_perp,v_perp)
       >= k ||v_perp||a^2,
    k=(E*-E0)/(E*-c)>0.

Since ||v||a^2=(E0-c)|t|^2+||v_perp||a^2, one valid constant is

    beta0=k/[1+(E0-c)||ell||a*^2/|ell(psi0)|^2] >0.

Solve u_N=b0+v_N, v_N in W_N, with b_E0(t,u_N)=0 for t in W_N. Coercivity gives a unique solution and the finite-dimensional variational characterization: u_N uniquely minimizes b_E0(u,u) on ell(u)=1 within V_N. Orthogonality and coercivity give

    ||u_N-psi||a <= sqrt(M0/beta0)
                        inf_{v in b0+W_N} ||v-psi||a ->0,
    M0=sup_{lambda in spec H} (lambda-E0)/(lambda-c) <=1.

More directly the error is the best approximation in b_E0 on W. The nonnegative defects d_N=b_E0(u_N,u_N)=b_E0(u_N-psi,u_N-psi) decrease to zero. There is no universal N-rate here; this is a best-approximation rate bound. beta0 is positive but not numerically certified in this task.

The plain coefficient pencil is K_N(E)=-[H_N-E G_N] up to the same positive physical normalization, with symmetric raw coefficients conjugated by the multiplicity weights. For a finite polynomial the exact differential recurrence equals the weak form by integration by parts. Thus deleting the constant row, imposing A000=1 and zero tail is exactly this affine Galerkin problem, not an arbitrary boundary rule. If g_N is the raw dropped P row, the weight convention gives d_N=-Ceps*g_N. Only its sign and proportionality matter if a common normalization was suppressed. g_N tends to zero. The other rows vanish exactly through N; no assertion is made about their omitted exterior rows at finite N.

At the infinite limit, density supplies every physical test function. Equivalently, stationarity on W says b_E(u,.) is a scalar multiple of ell; testing u (ell u=1) shows that scalar is b_E(u,u). At E0 its vanishing is equivalent to the full weak eigen-equation. Simplicity and normalization give uniqueness. Thus within form-convergent normalized reconstruction, admissibility and this backward limit agree. This proves the ground-state criterion, not a scalar component-ratio Pincherle theorem for every formal solution or every eigenvalue.

PROVED tail stability: if a prescribed exterior polynomial tau_N (shells N+1 through N+3) has ||tau_N||a->0, its induced change h_N in the normalized interior solution lies in W_N and satisfies b_E0(t,h_N)=-b_E0(t,tau_N). Consequently ||h_N||a<=M0/beta0 ||tau_N||a->0. Raw last-shell amplitudes alone do not imply this hypothesis; high-shell kinetic/Gram weights and the increasing number of coordinates matter. Also ||u_N-psi||a<=sqrt(d_N/beta0). These are norm-based stability/error bounds, not a certified numerical shell exponent.

## 5. PROVED: nonreal tails, resolvents and spectral recovery

On W, represent L2 pairing by the bounded positive self-adjoint operator T in the a inner product: a(t,Tv)=<t,v>. Then b_z corresponds to I-(z-c)T. Its spectrum is bounded away from zero for Im z!=0, since spec T lies in [0,1/(E0-c)]. The same bound holds on every W_N. With

    beta(z)=inf_{0<=mu<=1/(E0-c)} |1-(z-c)mu| >0,
    M(z)<=1+|z-c|/(E0-c),

Galerkin quasioptimality ||u_N(z)-u(z)||a <=[1+M(z)/beta(z)] inf_{v in b0+W_N}||v-u(z)||a proves convergence. The +1 bound follows by estimating u_N-v with the finite inverse. The limit solves b_z(t,u)=0 on W and ell u=1. Its source is a multiple of the form dual ell; ell is not an L2 load in general. For nonreal z that multiple cannot vanish, because H has no nonreal eigenvalues. A convergent normalized tail solution off spectrum is therefore expected, not a failure of the method. At arbitrary real z tail poles can occur: no unconditional convergence claim is made there.

For unrestricted V_N and any f in L2, solve b_z(t,R_N(z)f)=<t,f>. The same argument gives form convergence for each nonreal z, as well as for c. In particular the embedded finite resolvents converge strongly. For plain polynomials group three adjacent shells: K_N(z) is block tridiagonal, with growing rectangular off-diagonal blocks. Eliminate from the last block via S_i=D_i-U_i S_(i+1)^-1 V_i, including the correspondingly eliminated load. All principal tails are invertible nonreal (or below -4), by their physical positive Gram matrix. This computes exactly R_N, requiring no off-diagonal inverse and no growing-block Pincherle theorem.

Ground Rayleigh--Ritz eigenvalues E_N decrease to E0 by the form core. Second min--max values are >=E*, so normalized ground Ritz vectors, phase aligned, converge in form norm: their perpendicular mass is <=(E_N-E0)/gamma and their perpendicular a norm squared <=(E_N-E0)[1+(E0-c)/gamma]. This gives rank-one eigenprojection convergence even in operator norm. For a fixed small contour around E0 below E*, eventual uniform distance from the finite spectra plus strong resolvent convergence justifies the contour/residue recovery of this projection. Strong resolvent convergence *alone* was not enough to assert absence of pollution at an arbitrary real contour; the ground min--max separation supplies that missing argument here. Use -(2 pi i)^-1 integral (H-z)^-1 dz with positive orientation.

A further fully physical algorithm is normalized shifted inverse iteration using the same backward Schur solves at c=-5. In infinite dimension, for f with nonzero overlap with psi0, the normalized (H-c)^-m f converges in L2 with ratio <=[(E0-c)/(E*-c)]^m times its initial overlap factor. Taking N->infinity for each fixed m in the finite solves, then m->infinity, gives psi0. f=w is positive and has nonzero overlap. This supplies a spectral iteration rate, not a shell-cutoff rate or a justified arbitrary interchange of limits. It also avoids evaluating a singular full matrix at E0.

## 6. PROVED scope / OPEN extensions

If V_N subset U_N subset Qcal are nested finite spaces, the core theorem is inherited immediately. Finite log powers times polynomials are H1 at the triple point (in six dimensions even gradients of log R are integrable) and exponentially decaying at infinity; collision faces and angular atoms still require inspection. Include only atom-polynomial products actually in Qcal, remove Gram dependencies, and retain the entire plain V_N sequence. This proves full physical Galerkin/Schur and ground Ritz convergence for those enriched spaces. It does not prove that every individual U,V letter or every multiplier is admissible. Bare V~s/t fails H1 at electron coalescence; exact t regularization removes that pole but other strata must still be checked.

The normalized real-tail theorem extends abstractly to U_N with any fixed bounded ell nonzero on psi0, provided the ell=0 tail spaces are conforming and dense. Enriched coefficient gauges can prevent the raw 'one coefficient equals 1' condition from being a physical ell. Their physical pencils usually have dense log/atom overlaps, so backward elimination is a general Schur algorithm rather than the formal sector's fixed-band recurrence. An infinite log tower requires a specified reconstruction topology and gauge; no completeness theorem for separately restricted log sectors is asserted.

OPEN: a universal N-rate needs a quantitative weighted Sobolev/polyhomogeneous approximation theorem for psi0 (and admissible enriched atoms), plus norm constants. A scalar/block ratio definition that classifies *all formal* growing-shell solutions needs reconstruction and an exponential-dichotomy/Weyl-subspace theorem; our result instead gives an actual convergent ground-state finite-section selector. Extending the exact coefficient-Miller proof to a particular enriched formal recurrence requires identifying its physical test pencil, normalization/gauge and dense admissible tail spaces. No merely interior multiplied PDE establishes collision admissibility. Energy convergence still gives no a21 continuity theorem.

## 7. PROVED interpretation of coefficient limits

Form convergence gives convergence of every fixed ell_j(u_N) to ell_j(psi). It does not by itself prove that the ordinary partial sums of the *limiting* Laguerre coefficient array converge in form norm: the finite-section arrays are cutoff-dependent Galerkin solutions, not literal partial sums of one array. The selector proved here reconstructs the physical state by its form-convergent finite-section sequence. A global raw-series basis theorem or a classification of all formal sequences remains an additional question. The canonical finite-head limits are unambiguous in the plain representation.

PROVED scalar spectral test near E0: because the lowest Ritz vectors converge in form norm and ell(psi0) is nonzero, their constant coefficients are eventually nonzero. The Schur determinant identity then makes the nearby zero of the dropped row g_N(E) exactly the lowest Ritz eigenvalue E_N, as long as the tail pencil is invertible there. Tail coercivity at E0 persists in a fixed neighborhood by the uniform form bound, so this condition holds eventually as E_N->E0. Thus these zeros converge monotonically to E0. This conclusion concerns the ground branch, not every zero or arbitrary real tail pole.
