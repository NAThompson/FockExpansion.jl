# TASK_codex_31 — exact characteristic flux kernel and panel-source route

2026-10-24. Serial work; no agents, commits, repository edits, or changes to
inputs. New deliverables are x31_*; shared logs are append-only.

**PROVED / VERIFIED:** the proposed K32 is the endpoint-normalized regular
Sturm product kernel, with its characteristic support interpreted as a
locally integrable distribution. In particular it gives the first-log flux
functional for an analytic start and integrable polylogarithmic sources.
The endpoint delta is a weak limit, not the pointwise formula evaluated at
Y=b. No elliptic rotation average remains in this kernel. This does NOT
imply that its integral against a panel source is elementary or nonelliptic.

**VERIFIED / OPEN:** the Z-linear panel 1/2 source has two explicit 23-state
presentations, and its exact output rows are exported here. Thus a finite
holonomic input for B_P is available. A two-variable relative direct image
of these modules is a concrete smaller route than an expanded q/atan chain.
Setup is feasible on a small-hours scale; completion of an exact scalar
annihilator and its physical boundary certificate in hours is not established.
No annihilator, elliptic-factor decision, or decision of Conjecture A is
claimed.

## 1. Coordinates and exact operator identities — PROVED / VERIFIED

Put a=pi/2, b=3pi/2, s=b-X, r=b-eta, t=b-Y, all in (0,pi). Set

    L_s = -d_s(sin(s) d_s) - (c/16)sin(s),
    z = cos(s/2)cos(r/2)cos(t/2)
        + sin(s/2)sin(r/2)cos(t/2)
        + sin(s/2)cos(r/2)sin(t/2)
        + cos(s/2)sin(r/2)sin(t/2).

This z is exactly the task's z0 and is symmetric in s,r,t. The operator L_s
is M_X, with identical formulas for r,t. For any one-variable F, the exact
chain-rule identities give

    (L_s-L_t)F(z) = (sin(s)-sin(t))/4 * Q_c F(z),
    Q_c = (z^2-1)D_z^2 + 3z D_z - c/4.

The other two pairs follow by symmetry. More explicitly,

    -sin(s)z_s^2 + sin(t)z_t^2
        = (sin(s)-sin(t))(z^2-1)/4,
    -sin(s)z_ss-cos(s)z_s + sin(t)z_tt+cos(t)z_t
        = 3(sin(s)-sin(t))z/4.

For c=32, F=(4z^3-3z)/sqrt(z^2-1) satisfies Q_32 F=0.
`x31_check.py` verifies all three pairs as exact rational identities, not
floating-point substitutions. It uses p=tan(s/2), q=tan(r/2), h=tan(t/2),

    N=1+pq+ph+qh, D^2=(1+p^2)(1+q^2)(1+h^2), z=N/D.

The radical D cancels from each checked coefficient identity.

## 2. Characteristic front and distributional equation — PROVED

Define

    K(s,r,t) = (1/(2pi)) 1_{z>1} (4z^3-3z)/sqrt(z^2-1).

At e=z-1=0 its expansion is

    K = (1/(2pi sqrt(2))) e_+^(-1/2)
        * (1 + (35/4)e + (315/32)e^2 + ...).

Thus there is an inverse-square-root wavefront, not a finite jump. Its
coefficient, together with endpoint normalization below, fixes the constant.
The conormal principal symbol of L_s-L_t on z=1 vanishes by the first
identity in section 1: this front is characteristic for every operator pair.
There is no additional front at z=-1: the chosen support never reaches it.

To differentiate, use the analytic distribution family e_+^lambda and its
identity D e_+^lambda=lambda e_+^(lambda-1), continued to the nonintegral
half-integers. Apply Q_c term by term to this convergent local expansion.
The same scalar ODE holds distributionally; the coefficients cancel and
there is no delta or delta derivative. Multiplying a formal derivative of
Theta by its divergent boundary value is invalid and is not this calculation.
The pullback is legitimate in the open cube: at a front point with z=1,
the full gradient in (s,r,t) is nonzero. Indeed, for fixed s,t the two r
roots are simple except when a root touches r=0 or pi, which is outside
the open cube. Hence

    (M_X-M_Y)K = (M_eta-M_Y)K = 0

as distributions in the open cube, including across the front. Equality
of the three individual operators means this pairwise distributional
identity, not that each individual M annihilates K.

## 3. Exact cap substitution, endpoint delta and weak boundary flux — PROVED

For fixed interior s,t let

    d = sin(s)sin(t), A=cos((s-t)/2), B=sin((s+t)/2),
    R=sqrt(1+d), r0=2 atan2(B,A), alpha=acos(1/R).

Then z=R cos((r-r0)/2), and the complete support interval is
r0-2alpha < r < r0+2alpha. It lies in [0,pi]. This also follows from
z(0)=cos((s-t)/2)<=1 and z(pi)=sin((s+t)/2)<=1, with the unique interior
maximum R>1. Contacts occur at r=0 when s=t and at r=pi when s+t=pi.

The following change of variable covers the whole cap:

    r(theta)=r0+2 asin(sqrt(d/(1+d))*sin(theta)),
    -pi/2 < theta < pi/2,
    z(theta)=sqrt(1+d*cos(theta)^2).

It gives the particularly simple exact measure

    K(s,r,t) dr = (1+4d*cos(theta)^2)/pi dtheta.             (1)

This both removes the square-root singularities and proves

    int_0^pi K(s,r,t) dr = 1+2sin(s)sin(t) <= 3.           (2)

For a continuous test function f, as t decreases to zero, d tends to zero
and r(theta) tends uniformly to s. Dominated convergence in (1) gives

    int K(s,r,t)f(r) dr -> f(s).

As t increases to pi, r(theta) tends uniformly to pi-s, giving f(pi-s).
Consequently

    K(b;X,eta)=delta(X-eta),
    K(a;X,eta)=delta(eta-(a+b-X)),

in the weak sense, with delta measured in d eta. There is no factor two,
pi, or endpoint half-weight for interior X. By symmetry the same results
hold when X or eta tends to either endpoint.

For f analytic near [0,pi], its cap action U(s,t) is analytic up to each
endpoint in either s or t when the other variable is interior. In (1), r0
is analytic, and powers of sqrt(d) with odd exponent integrate to zero
by theta -> -theta. The arcsin argument is bounded by 1/sqrt(2), so there
is no hidden inverse-sine endpoint singularity. In particular, at t=0,

    U(s,t)=f(s)-t M_s f(s)+O(t^2),
    sin(t)partial_t U(s,t) -> 0.

The analogous reflected statement holds at pi. These give the regular
weak endpoint conditions for K and justify endpoint integration by parts.
One must use weak limits: at the isolated contacts s=t or s+t=pi, a
pointwise r-flux of the displayed density need not have a finite limit.
The distributional boundary flux still vanishes. Equation (1) proves it
after testing against analytic endpoint germs. Thus those contacts do not
supply additional boundary deltas.

## 4. Identification with the Sturm product kernel — PROVED

Let h_j be a nonzero eigenfunction of the regular/log-free Sturm realization
of M on I, with eigenvalue lambda_j, in the x27 normalization. In s coordinates
write f_j(s)=h_j(b-s). Its endpoint germ is analytic and f_j(0)!=0: the
regular recurrence for (L_s-lambda_j)f=0 forces the whole germ to vanish if
the constant term vanishes.

Set U_j(s,t)=int K(s,r,t)f_j(r)dr. The distributional identities and weak
zero-flux endpoint conditions above imply

    (L_s-lambda_j)U_j=0,   U_j(0,t)=f_j(t).

The regular endpoint recurrence uniquely selects the analytic solution, so

    U_j(s,t)=f_j(s)f_j(t)/f_j(0).

Returning to X,Y,eta proves the exact product identity for EVERY regular
mode, not just the three elementary tests:

    int_I K(Y;X,eta)h_j(eta)d eta = h_j(X)h_j(Y)/h_j(b).

By completeness of the same regular self-adjoint Sturm realization used
in x27, this identifies K with

    sum_j h_j(Y)h_j(X)h_j(eta)/(||h_j||^2 h_j(b)),

as a weak spectral kernel. No absolute pointwise convergence of this series
is asserted or needed. Equations (1)-(2) independently define the bounded
cap measure and avoid relying on numerical collocation convergence.

`x31_numeric.py` independently checks the three exact modes at ordinary
points and both contact loci, 65-digit arithmetic: worst error below 1e-64.
A general analytic test also checks the weak operator identity at these
loci, with errors below 1e-63. Endpoint and reflected-endpoint errors decrease
linearly through t=1e-12; weak endpoint flux tends to zero. These are
numerical verifications, not the proof of uniqueness or distributions.

## 5. Flux, analytic start, and logarithmic source endpoints — PROVED with stated hypotheses

For a mode of the continued inhomogeneous solution, the exact identity is

    d_X[cos(X)(h_j f_j' - h_j' f_j)] = h_j R_j.

If f_j is analytic at a, its start flux is zero. At b, put epsilon=b-X.
If R_j=O(1+|log(epsilon)|^m), or more generally h_j R_j is integrable,
the Wronskian flux has a finite endpoint limit. For the logarithmic branch

    f_j = B_j log(-cos X) + less singular terms,

that limit is h_j(b)B_j. Integrating h_j R_j over an epsilon endpoint
interval gives O(epsilon*(1+|log epsilon|^m)); the remaining local terms
are epsilon times powers of log and have zero additional limiting flux.
Thus

    B_j = int_a^b h_j(X)R_j(X)dX / h_j(b),
    B(Y) = int_I int_I K(Y;X,eta)R(X,eta)d eta dX.

The sign here is for +log(-cos X); code using -log(-cos X) reverses it.
There is no integration-by-parts boundary term attached to an integrable
logarithm in R itself. A 1/epsilon source or a nonregular start would
require a different endpoint analysis; neither is silently included.
For a nonzero start flux W_j(a), add W_j(a)/h_j(b) to the result.

The double integral is absolutely defined for the usual polylog endpoint
source class. For example, if |R(s,r)| is bounded by a finite sum of
products g(s)h(r), with g,h finite powers of endpoint logarithms, (2),
symmetry, and Cauchy-Schwarz for the positive measure K ds dr give

    int int K |g(s)h(r)| ds dr <= 3 ||g||_2 ||h||_2.

Endpoint strips tend to zero by the same bound with the appropriate
indicator. Thus source logs and the cap's corner contacts do not hide an
extra source face. This covers finite products of endpoint logs, not
uncontrolled power singularities of individual unrecombined atoms.

This theorem applies to the specified source contribution R. It does not
remove upper-log hierarchy terms if they occur elsewhere in the physical
order-four equation, nor replace an unvaried source by its variation.
B_P is the panel-source contribution; equality to an individual ordered
path B_A is not inferred. Continuation/framing of the source must still
be the intended one.

## 6. General c — PROVED

Let nu^2=1+c/4, and T_nu(z)=cosh(nu arcosh z) for z>=1. Then

    K_c = 1_{z>1} T_nu(z)/(2pi sqrt(z^2-1))

has the same distributional equations, regular weak endpoint conditions,
and unit endpoint delta. The exact substitution identity is

    Q_c[T(z)/sqrt(z^2-1)]
      = ((z^2-1)T''+zT'-nu^2 T)/sqrt(z^2-1).

At the endpoint T_nu(1)=1, so the delta normalization is unchanged. The
preceding Sturm product proof applies. For integral nu=n, T_n is the
Chebyshev polynomial and the density is algebraic. The parameter dictionary
is c=4(n^2-1). Thus a notation c=4k(k-2) gives n=k-1; a notation
c=k(k+4) gives n=(k+2)/2. The proposed T_(k-1) needs the first dictionary;
it is not a universal meaning of k.

## 7. Exact Z-linear panel source — VERIFIED, explicit exports

Use the x16 conventions r1=x, r2=y, r12=c_dist, Q=x^2+y^2, and its normalized
moment rows F0out,F2out. The symbol c_dist here is not the operator parameter.
For panel 1, H=x^2/y^2 and k=(c_dist^2-Q)/y^2. Extracting Z from the actual
psi30 coefficients gives

    [Z]c0 = sqrt(2)(bpanel^2-1)/72,
    [Z]c2 = sqrt(2)(1-bpanel)(9bpanel+16Y2-8)/144,
    bpanel=(x^2-y^2)/Q, Y2=c_dist^2/(2Q).

Since N_i=R/(sqrt(2)x) Fiout, these simplify to

    R^3 P1 = -x*y^2/18 * F0out
             + y^2*(x^2-17y^2+8c_dist^2)/(72x) * F2out.   (3)

Panel 2 exchanges x,y. This is a Z-coefficient extraction, not setting Z=1
in the helium Z=2 output. In the angular normalization R=1 no radial factor
remains. Multiplication by wv0/8 is algebraic:

    wv0/8 = sqrt(2)/4 sin(alpha)cos(beta/2),
    alpha=(X+eta)/2, beta=(eta-X)/2.

`x31_source.py` applies x16's exact endpoint/odd-moment reductions to (3).
`x31_source_output.txt` exports both rational output rows (20 terms each).
`x31_source_summary.json` verifies that each closes on exactly the SAME
23 graph-reachable states, with derivative nnz (84,84,69). This is a
sufficient module, not a proof that the 23 functions are independent or
that their integrated direct image has rank 23. The script reuses the
previously certified connection; it does not re-prove all its entries.

Continued distances are x=cos(alpha/2), y=sin(alpha/2),
c_dist^2=1-cos(alpha)cos(beta). Choose their physical germs and continue
on the prescribed strip. In particular x can change sign. A naive
positive-distance or principal-atan evaluation everywhere on the outside
square is not the same framing. The rational connection identities
continue, even when an undeformed real-q presentation crosses a branch cut.

## 8. Smallest exact telescoping route and realistic hours assessment

**PROVED setup reduction.** The exact rational cap polynomial is

    Fcap=N^2-(1+p^2)(1+q^2)(1+h^2),
    degree_(p,q,h)=(2,2,2),
    disc_q(Fcap)=16hp(1+p^2)(1+h^2).

On the positive p,q quadrant the support is Fcap>0, and

    K32 = N(4N^2-3D^2)/(2pi D^2 sqrt(Fcap)),
    dX d eta = 4 dp dq/((1+p^2)(1+q^2)).

The signs from the two reversed endpoint maps cancel. The radical D itself
has disappeared from K. This is one quadratic cover for the kernel, with
rational numerator and denominator. Source distances and wv0 introduce
additional finite algebraic covers; they must be retained or rationalized.
The displayed cap discriminant alone proves neither an essential elliptic
factor nor its absence in B_P.

**Concrete minimum work.** Use the existing 23-state panel connection,
pull it to the chosen algebraic (p,q) chart, and twist by the scalar kernel
and Jacobian. Over that algebraic coefficient field the twist adds
(d log kernel/Jacobian) times the identity. It does not introduce q-moment
or atan-representation integration variables. Work separately on the two
23-state blocks, or on their direct sum with the single observable P1+P2;
no unproved symmetry identifies the two physical periods. The observable
cyclic submodule may be smaller than the full 46-state sum, but its size
has not been measured.

Compute the relative direct image in the two variables p,q, with h as the
output parameter. A staged reduction in one variable followed by the other
can be benchmarked before a simultaneous reduction. The required certificate
has the form

    L_h(g) = partial_p(C_p) + partial_q(C_q)

on the pulled module/cover, together with all boundary contributions and
physical germs. Here g is the source observable including the scalar measure.
The cap front must be handled in the twisted/relative complex, equivalently
by continuation of Fcap_+^lambda to lambda=-1/2. Do not just set a divergent
primitive to zero on Fcap=0. Certificates with normal flux O(sqrt(Fcap))
have a vanishing front face; otherwise retain and reduce the face term.
At p,q=0,infinity retain the regular/log boundary data until they are shown
to vanish. Also retain source-contour faces if an expanded q presentation is
used. Apparent poles of the reduced connection are charts, not permission
to introduce discontinuities in the actual analytic source.

A homogeneous scalar annihilator follows only after annihilating any
inhomogeneous face periods as well. Then use x29's exact factor/Hom test
in the matching output coordinate/twist, after any algebraic descent. An
operator for just an interior integrand is not an operator certificate
for the physical B_P.

**Feasibility verdict — CO-SIGN exact route, QUALIFY timing.** Exact source
extraction and closure checks finished in seconds with peak RSS about
110 MB. The kernel differential proof takes seconds. These facts remove
the previous missing-kernel and missing-aggregate-source obstacles.
A bounded prototype for a first relative reduction is consequently a
reasonable hours-scale experiment. They do NOT bound the cohomology rank,
scalar output order, coefficient degrees, or the cost of the second direct
image. The 23-state input rank is not an upper bound on those quantities.
No direct-image reduction or actual telescoper was computed in this task.
There is therefore no evidenced promise of a complete exact B_P operator
within hours or 1.5 GB. Conversely, previous failed numerical ODE searches
for B_A do not supply a complexity lower bound for B_P.

The smallest justified next computation is this two-variable module
pushforward, beginning with one 23-state block and explicit face tracking;
a generic expanded (X,eta,q,atan) telescoper adds avoidable dimensions.
The exact source rows here make that experiment reproducible. It remains
an OPEN computation, not a demonstrated obstruction or a completed proof
of an elliptic constituent.

## 9. Co-sign, dissent, and verification

**CO-SIGN:** Claude's supported algebraic K32, delta normalization, and
regular-mode product identity, now proved for all regular modes. The
observation-circle degeneration is consistent with this formula. The
characteristic support and normalization are essential.

**DISSENT / qualification:** the untruncated raw rotation average remains
excluded by x29; that exclusion does not apply to this supported cap kernel.
Neither algebraicity of K nor availability of x16 proves B_P elementary,
nonelliptic, or small-order. The previous blanket “no exact holonomic
source” objection applies to an individually unresolved phi3^A, not to
this specified aggregate panel source. An exact B_P result still would
need its relation to the full hierarchy stated before deciding A.

Fresh verification commands (all bytecode disabled):

    /Applications/SageMath-10-3.app/Contents/Frameworks/Sage.framework/Versions/10.3/venv/bin/python3 -B codex/phase13/x31_check.py
    /Applications/SageMath-10-3.app/Contents/Frameworks/Sage.framework/Versions/10.3/venv/bin/python3 -B codex/phase13/x31_numeric.py
    /Applications/SageMath-10-3.app/Contents/Frameworks/Sage.framework/Versions/10.3/venv/bin/python3 -B codex/phase13/x31_source.py
    python3 -B codex/phase13/x31_audit.py final

Results are in x31_check.out, x31_numeric.out, x31_source.out and
x31_audit.out. Scientific jobs were serial. Peak self-measured RSS is
recorded in each scientific check; final audit preserves all protected
source/input hashes and both shared-log prefixes. The runner automatically
updates codex_31.log outside the task output prefix; this harness transcript
is reported separately and no task command writes it. The first audit
flagged that live transcript, then the audit was refined to classify only
that specific runner log separately; its original result is retained in
x31_audit_first_attempt.out. The initial default-Python run
failed because sympy is absent; the completed checks use installed Sage
Python. macOS time -l could not read kern.clockrate, so its output is not
a peak-memory measurement. No unverified run supplies a positive claim.
