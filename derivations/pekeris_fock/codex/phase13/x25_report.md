# TASK_codex_25 — collinear nearby-cycle fingerprint audit

2026-10-24. Independent serial audit [X]. Claude's inputs [C] were read
only. PROVED refers to the arguments below, VERIFIED to finite checks,
CITED to external inputs, CONJ to proposals, OPEN to missing evidence.

## 1. Disposition

**CO-SIGN (1), over rational cohomology and with an isogeny correction.**
The generic collinear degeneration is type II; its rank-five transcendental
part has Jordan blocks `(2,2,1)` and weight-graded dimensions `(2,1,2)` in
weights `(1,2,3)`. The residue curve is a QUARTIC in q. The Legendre curve
with lambda=(u-1)/(u+1) is its degree-two isogenous quotient, not in general
the same elliptic curve. Its periods require an algebraic prefactor.

**CO-SIGN (2) for the full specified old category.** The five order-two
atoms, eighteen classical atoms, all four sets of twenty-three panel atoms,
and algebraic factors have only Tate/finite and old-panel elliptic pure
constituents. All old-panel elliptic factors have mixed-Tate nearby cycles
at the generic circle. Tensor products and extensions preserve this
exclusion. None supplies the varying residue elliptic H1 there.

**CO-SIGN the sufficient fingerprint in (3); DISSENT from its literal
scalar-limit converse.** A certified elliptic constituent of a physical
branch's normalized log coefficient proves nonmembership in C3. Occurrence
of T5 in any composition position forces such an elliptic constituent in
the IMAGE OF N on nearby cycles. It does not force a nonzero zeroth normal
coefficient of N applied to the scalar function. Positive normal vanishing
can hide it; derivatives/leading coefficients repair this defect. The test
does not by itself uniquely identify the parent factor as T5.

**OPEN (D), enlarged square and actual psi40 fingerprint.** Physical-edge
analyticity follows from analytic elliptic regularity and rotational
invariance under the physical hierarchy hypotheses. That local theorem
does not establish holomorphy of the continued branch throughout
|tx|,|ty|<3. At the audit's final inspection claude11 contained only
cheb_decay.out; NOTES.md and claim_*.md did not exist. No physical log
coefficient or numerical certificate was available to co-sign. Conjecture A
for psi40 remains OPEN.

## 2. The actual limit curve and its form

Use delta=1-u²-v², P=q⁴-1, D=q⁴+1, A=1-v² and

    B=-uv(q⁴-1)+2(1-v²)q²,   C=D²-[u(q⁴-1)+2vq²]².

**PROVED / VERIFIED:** B²-AC=-P² delta, so

    Y²=X[(X+B)²+P² delta].                              (2.1)

At delta=0 its normalization is rational ruled:

    X=a²,   Y=a(a²+B).

The two preimages of the nodal locus form the conductor curve

    E_q: a²=-B=uvq⁴-2(1-v²)q²-uv,   u²+v²=1.           (2.2)

For uv≠0 this quartic has four distinct roots generically. Moreover its
roots avoid P=0: there B=2Aq²≠0. Thus the cover (2.2) is a smooth genus-one
curve on every generic arc considered here. Exclude u=0, v=0 and their
complex/algebraic equivalents when selecting the arc; those are further
degenerations, not a test of a generic type-II limit.

The normalized two-form is

    omega=dq wedge dX/Y=2 dq wedge da/(a²+B).

Its residue at a chosen conductor preimage is dq/a, with the opposite sign
on the other preimage. Equivalently, a positively oriented vanishing fibre
circle, swept over a lifted q-cycle gamma, has period tending to

    2 pi i integral_gamma dq/sqrt(-B).                  (2.3)

Orientation changes its sign. The q-cycle is a cycle of the double cover,
or a q-contour with the corresponding sign local system, not an arbitrary
ordinary P1 cycle with a globally chosen square root.

**PROVED isogeny correction.** The free involution `(q,a)->(-q,-a)` has
quotient

    r=q², Z=qa,
    E_r: Z²=r[uvr²-2(1-v²)r-uv],
    pullback(dr/Z)=2dq/a.

The involution is fixed-point-free on the smooth compact curve generically
(also at q=0 and infinity), hence this is an unramified degree-two isogeny.
The cubic's branch points are 0, infinity, r+=(u+1)/v, r-=(u-1)/v.
With x=r/r+ it is Legendre, up to algebraic scaling, with

    lambda=r-/r+=(u-1)/(u+1)=-tan²(phi/2), u=cos(phi).

In particular dq/a pulls back from

    [2 sqrt(uvr+)]^(-1) dx/sqrt[x(x-1)(x-lambda)].       (2.4)

The square root sqrt(uvr+)=sqrt[u(u+1)] must be continued on the chosen arc.
Cross-ratio permutations give equivalent Legendre presentations elsewhere.
The two j-invariants, independently checked, are

    j(E_q)=64(3-4u²)³/(1-u²),
    j(E_r)=64(u²+3)³/(1-u²)².

They are different generically. Both give the same rational H1 variation
up to isogeny; an integral cycle basis and its factors of two must not be
identified without specifying the isogeny map. E_q has the same j as the
second elliptic factor in the ENDPOINT K3 of x20. That auxiliary endpoint
K3 factor was never a constituent in the old order-three alphabet; its
appearance here does not refute item (2).

The periods in the normalization (2.3) satisfy the explicit tangential PF
equation

    [4u²(u²-1)D_u²+4u(3u²-1)D_u+(4u²+1)]p=0.          (2.5)

This follows exactly from the Gauss equation for K(lambda) with the gauge
p=[u(u+1)]^(-1/2)F(lambda); x25_exact.py verifies the gauge coefficients.
It gives a concrete proposed coefficient test, including its algebraic
normalization. A physical coefficient may instead be a rational/algebraic
differential transform of this system, so failure of this particular
scalar equation alone is not an exclusion of the H1 constituent.

## 3. Type II, weights and an independent operator check

**PROVED using CITED degeneration inputs.** One can see the standard
type-II.b construction directly in (2.1). After delta=w², blow up the
nodal ideal `(Y,X+B,w)`. The exceptional chart with
`X+B=wZ, Y=wV` has central conic bundle

    V²=-B0(q)[Z²+P(q)²].

It has the rational section Z=iP,V=0. Its compactified component and the
strict normalization are rational ruled surfaces. Their double curve is
the cover a²=-B0 above; the branch points are simple and avoid P=0.
Resolving the isolated special fibres supplies rational strata; the
genus-one double curve persists. This is the non-split nodal-fibre
type-II.b model. The standard semistable construction, including permissible
base change/birational modifications, is described in
[Lee–Weigand, §2.2, especially the type-II.b description](https://arxiv.org/pdf/2112.07682).
The chart here identifies its bisection in this family; an unreduced nodal
Weierstrass central fibre is not itself being called a Kulikov model.

For a type-II K3 limit the weight-one part is H1 of the elliptic double
curve, N²=0, and rank N=2. The weight filtration on full H2 has dimensions
2,20,22 in weights ≤1,≤2,≤3. These standard inputs follow from the
weight spectral sequence/Clemens–Schmid calculation in
[Morrison, §4(d)](https://web.math.ucsb.edu/~drm/papers/clemens-schmid.pdf).

The rank-17 algebraic lattice of x20 has finite monodromy, hence zero N
after a finite cover. Removing it gives, on T5,

    W1=Im N, dim 2; W2=Ker N, dim 3; W3=T5, dim 5,
    Gr1 ≅ H1(E_q), Gr2 ≅ Q(-1), Gr3 ≅ H1(E_q)(-1).     (3.1)

Use the usual cohomological H1 notation in (3.1); homology/period duals
give the corresponding dual identifications. Tate or finite twists from
markings do not affect the underlying elliptic differential constituent.
N maps Gr3 to Gr1(-1). Equation (2.3) identifies the limiting scalar
periods of Im N. Ramifying delta=w² doubles log monodromy: N_w=2N_delta.
It changes neither its rational image nor its Jordan block sizes.

**VERIFIED exact operator statement, independent of that geometric input.**
x25_exact.py reads op_e2_A.txt and decodes the actual op_e2_Q.sobj with a
restricted inert constructor decoder. They agree coefficient for coefficient;
no Sage installation/cache writes are needed. This does not newly certify
the reconstructed operator's characteristic-zero telescoping identity.

Put z=u², t=z-8/9 and divide the transformed L5 by u. The indicial
polynomial is

    -(32768/9) rho²(rho-1)²(rho-2).

For an analytic solution h=a0+a1 t+a2 t²+..., the all-orders recurrence
pivot for a_(n+3) is

    -(32768/9)(n+1)(n+2)²(n+3)².

Thus exactly three analytic germs exist. Compute the nonlogarithmic part
of L5(h log t). Its ONLY negative-power obstruction is

    [2/(9t)](23571a0+195840a1-65536a2).                 (3.2)

When this vanishes the same nonzero pivot constructs an analytic g with
L5(h log t+g)=0. The two-dimensional seed plane (3.2), plus the three
analytic germs, is a full five-dimensional Frobenius basis. Consequently
there are exactly two log blocks of size two, one block of size one, and
no log-squared term. This proves `(2,2,1)`, N²=0, rank N=2 for the supplied
operator; repeated exponents alone would not have proved it.

## 4. Two-cycle numerical comparison, beyond 30 digits

**VERIFIED numerical, not a ball enclosure.** At v=1/3, u0=2sqrt(2)/3 let

    a=sqrt((u0+1)/v), b=sqrt((1-u0)/v), lambda=-b²/a².

The roots of (2.2) are ±a, ±ib. x25_periods.py independently integrates
two q-contours, continuously tracking the square root through factored
analytic square roots. A is |q|=1 around ±ib. B is a circle around a and
ib, centered at (a+ib)/2, with radius the geometric mean of the second and
third root distances. Neither contour passes through a branch point.

The limits, INCLUDING 2pi i, agree with

    p_A=2pi i * 4K(lambda)/sqrt[u0(u0+1)],
    p_B=2pi i * 2iK(1-lambda)/sqrt[u0(u0+1)]

on mpmath's specified principal K branch. Their values are

    p_A = 28.9585710340215905509452610833194424353516008490... i,
    p_B = -29.0942725543807361061620339207208623241825062774...
           +14.4792855170107952754726305416597212176758004245... i.

To compare with the ANALYTIC Im-N subspace of L5, not just an arbitrary
constant multiple at one point, the calculation integrates normal jets
through t³. They come from the actual vanishing fibre-circle expansion

    2pi i F^(-1/2) sum_(k>=0)
      binom(-1/2,k)binom(-1/2,2k)[P² delta/F²]^k,
    F=-B, u=sqrt(8/9+t), delta=-t.                    (4.1)

The integrated jets satisfy (3.2) and the first analytic ODE recurrence

    a3=(8505a0+373590a1+550912a2)/131072.

Their (a0,a1) determinant is -627.877102776071358553395380108803252845...
i, hence they span the two-dimensional log-coefficient seed plane. The
cycle-A checks have errors below 1e-85; the tighter cycle-B contour,
4096 versus 8192 trapezoid nodes at 95 digits, has recurrence/refinement
errors far below 1e-70. Final raw residuals and K comparisons are in
x25_periods.json. These are numerical consistency checks with ample
30-digit margin, not certified integration error bounds. The geometric
residue argument supplies the identification; numerical jets do not
independently certify an entire two-variable Gauss–Manin connection.

## 5. Full old-category inventory and degeneration

Here C3 means the differential category generated by the stated old atoms,
allowing algebraic coefficients/finite covers, subquotients, finite tensor
operations and differential extensions. It does not allow arbitrary new
parameter direct images: those can introduce new pure constituents and
would change the definition.

| Entire generator family | Pure constituents, up to multiplicity/twist | Generic circle limit |
|---|---|---|
| order-two five-state `(1,L,alpha,V,U)` | finite algebraic characters and Tate; angles/logs/dilogarithms | finite/Tate, including logarithmic extensions |
| eighteen classical atoms | finite characters and Tate; L/T Clausen combinations, logs, angles, their stated squares | finite/Tate |
| panel 1, all 23 reachable atoms | H1(E1) copies and Tate/finite boundary factors | nodal normalization P1: mixed Tate |
| panel 2, all 23 | H1(E2old), algebraically isomorphic to E1 by s inversion, and Tate/finite | mixed Tate |
| panels 3 and 4, 23 each | H1(E3old),H1(E4old), related by sign/inversion, and Tate/finite | mixed Tate |
| rational/algebraic/radical multipliers and endpoint values | finite characters; Tate/logarithmic endpoint extensions | finite/Tate after finite cover |
| tensor/dual/subquotient/extension closure | constituents of tensor words in the above pure factors | only finite/Tate graded constituents |

**PROVED from the actual definitions, not just the displayed two cubics.**
For the first two rows, all exponentials of the angle/log arguments are
algebraic on a finite distance/radical cover. U and the eighteen classical
atoms are finite-length classical polylogarithmic iterated integrals of
dlog forms. Differentiation lowers logarithmic length, leaving rank-one
algebraic/Tate constituents. The elementary V is a radical times an angle.
No elliptic factor is hidden in these two rows.

For the panels, use x19_atoms.jl and the exact x16 definitions. Every one
of the 23 states is a relative integral on

    y²=s⁴+k s²+H

of a rational differential times A(s)^n/n! or B(s)^n/n!, n≤2, or an
endpoint value; A'=1/(1+s²), B'=H/(H+s²). The list includes the T-pole,
U/V-pole, B/s and endpoint atoms, not merely complete elliptic periods.
Filter by the log length of A,B. Its graded coefficient sheaves are
rank-one algebraic/Tate. Relative direct-image exact sequences then give
only H1 of this same elliptic curve, puncture residues and endpoint Tate
factors. Higher-order poles of 1/T alter de Rham representatives, not
the compact curve; punctures and relative endpoints add no new compact
genus. This applies to every one of the 4x23 states, including the removed
dependent states if they are retained instead.

With unit radial scale the four unscaled quartics are

    Q1=(1-u)s⁴-2vs²+(1+u),
    Q2=(1+u)s⁴-2vs²+(1-u),
    Q3=(1+v)s⁴+2us²+(1-v),
    Q4=(1+v)s⁴-2us²+(1-v).

Each quadratic in s² has discriminant -4delta. More strongly, if
Qj=aj s⁴+bj s²+cj, the exact identity

    aj Qj=(aj s²+bj/2)²+delta                      (5.1)

holds for ALL FOUR panels. Thus generically the quartic becomes a square
of a quadratic with two simple roots. The central curve consists of two
rational components meeting at two points. The s->s² quotient instead
has a single nodal cubic; these are isogenous multiplicative models of
the same H1 type. Neither has an elliptic component in its semistable
limit. The smooth elliptic H1 has limiting graded pieces Q(0),Q(-1).
The normalizations of punctured components are P1 with algebraically
moving punctures, so they add only Tate/finite pieces and their extensions.
This also handles pole collisions and endpoint incidences on a generic
arc. Rational factors can add poles or shift normal orders, but cannot
create an elliptic simple factor.

For the compact/open/relative weight and filtration arguments used here,
the CITED input is the residue spectral sequence and strictness for
geometric mixed Hodge structures, together with relative functoriality;
see [Deligne II, §§2.3, 3.1–3.2](https://www.numdam.org/item/PMIHES_1971__40__5_0.pdf)
and [Deligne III, §8.3](https://www.numdam.org/item/PMIHES_1974__44__5_0.pdf).
The application to these explicit panel definitions is the argument above.
Multiplicities/minimal ranks of the sufficient 23-state presentations are
not needed and are not newly asserted.

The residue elliptic variation is not finite/Tate: lambda is nonconstant
on the complexified circle and is a finite algebraic pullback of Legendre.
Its connected monodromy is SL2 and remains so after finite covers. Hence
it cannot occur in any of the generic-circle nearby-cycle objects just
listed. This proves item (2) for the FULL category, not only its pure
elliptic generators or one choice of physical branch.

## 6. Fingerprint, converse, and which branches

**PROVED sufficient theorem.** Let f be the actual physical psi40 sector
germ, with its finite-rank geometric regular-singular continuation module.
Choose a generic divisor arc, a path from the physical basepoint to a
punctured tubular neighbourhood and a continued branch b. Extract normal
Frobenius coefficients with algebraic powers and log powers separated.
If a coefficient of log-monodromy N b has the varying H1(E_q) constituent
as a differential function of the tangential variable, then Df is outside
C3. Indeed C3 is closed under continuation, nearby cycles, N-images and
these normalized coefficient subquotients, whose graded factors were
excluded in §5. This proves a new constituent outside the old category.
It does NOT alone prove that its parent is specifically T5: a different
new higher-rank object could have the same elliptic degeneration.

Here N=log(T) on the unipotent part, without a 2pi i normalization; on T5,
N=T-1. For a general mixed physical module higher powers may occur.
The coefficient of a single scalar log is then not automatically N b:
use `exp[(log t)/(2pi i) N]` and keep all log powers. Algebraic prefactors,
normal vanishing orders and meromorphic poles must also be retained.

**PROVED module-level converse, including middle positions.** On a regular
generic arc, nearby cycles are exact on local systems. Suppose a filtration
piece F/G of the physical continuation space is T5 (or its period dual).
Then

    N(F) is a subobject of N(V),
    N(F) -> N(F/G) is surjective.

These are elementary linear-algebra identities for an invariant filtration;
they do not assume exactness of the image-N functor. By (3.1), N(F/G)
contains the residue H1. Hence Im N(V) has that elliptic constituent even
when T5 is in a quotient or in the middle of an extension.

If V is spanned by actual branches of f, their N-images span Im N(V).
It follows that at least one branch has an elliptic constituent among its
normalized nearby-cycle coefficients (allowing sufficient normal jets).
This statement uses the geometric regular-singular/cyclic realization;
it is not a claim about arbitrary irregular differential equations.

**Counterexample to the literal zeroth-limit converse.** Take a nonzero
period Pi of T5 and f=delta Pi. Over the rational coefficient field Df and
D Pi are the same simple module. For EVERY branch,

    Nf=delta NPi,   (Nf)|_(delta=0)=0,

because NPi is an analytic invariant period at the generic type-II limit.
Nevertheless T5 is present. Replacing delta by delta^k hides it to any
prescribed normal order. Dividing by that factor, or extracting suitable
normal derivatives/leading coefficients, recovers H1. This is a concrete
obstruction to any claimed fixed zeroth-coefficient necessary test. It
does not assert that the actual psi40 has this vanishing factor.

No particular physical-access branch, real quadrant, or far grid sheet is
forced by the converse. At any chosen generic complex circle arc the
nearby T5 inertia is nontrivial; some branch reached by global complex
continuation must detect its image in the corrected sense. The accessible
branches on tx=±3 or ty=±3 need not exhaust that monodromy orbit. A
physical branch can remain analytic at every physical edge while other
branches carry the factor, as already explained in x21/x23. A negative
test on the far real lines, especially only their constant log terms,
does not decide non-survival in quotient/middle positions.

## 7. Analyticity and the Chebyshev evidence

**CO-SIGN the local physical-edge argument, under its stated hierarchy
hypotheses.** Inductively, analytic elliptic regularity on the real S3
applies to each angular coefficient away from {e1,-e1,e2}, because the
source and elliptic operator are real analytic there; pure harmonic
subtractions are analytic too. The CITED analytic-regularity input is
[Morrey–Nirenberg (1957)](https://onlinelibrary.wiley.com/doi/abs/10.1002/cpa.3160100204).
Near a generic point of the fixed circle,
an SO(2)-invariant analytic function is analytic in a tangential coordinate
and the transverse radius squared. Substituting radius²=1-u²-v² gives a
holomorphic germ in u,v across the generic edge. This is the ordinary
analytic elliptic-regularity theorem, not a numerical consequence of
Chebyshev decay. The invariant-power-series argument supplies the descent
from S3 to the disk coordinates. It excludes the Coulomb vertices.

**OPEN / not co-signed:** this does not give the enlarged square. Parts of
that real alpha,beta square map to the complexified sphere, where the
real elliptic-regularity argument no longer applies. At cosine critical
lines the target map also ramifies. Analyticity at resolved D0 vertices
would itself need a proof involving the actual Coulomb singularities.

Sliding meridians is a valid way to compare conjugate local inertias in
the SAME smooth divisor stratum while transporting the SAME germ. It
does not prove that all branches have trivial inertia, and trivial
inertia alone does not exclude poles (1/t is a counterexample). A proof
of square-wide holomorphy must establish the complete relevant singular
locus, removable extension on the chosen sheets, treatment of crossings
and ramified Coulomb lines, and consistency of the continued germs.
Normal crossings by itself establishes none of those removability claims.
These are specific gaps to audit in a future NOTES/claim file; they do
not constitute a refutation of square-wide analyticity.

**VERIFIED elementary slope calibration.** With the centered definitions,

    delta=cos(pi tx/2)cos(pi ty/2),
    log10(3+sqrt(8))=0.765551370675726156638799048384... .

Thus 0.773 is not that exact exponential rate. A logarithm at t=3 has
Chebyshev coefficients proportional to `(3+sqrt(8))^(-n)/n`; its exact
secant slope between degrees 40 and 70 is

    0.765551370675726... + log10(70/40)/30
      =0.7736526389652693046482903306888... .

This quantitatively explains the reported finite-degree slope without
inferring a new period class. As an explicit OLD example, set
r1=cos(alpha/2), r2=sin(alpha/2), r12=sqrt(2)sin(beta/2) on their continued
signed radical sheets. The order-two atom L obeys

    exp L = 2sqrt(2) sin[pi(3+ty)/8] cos[pi(tx+1)/8].

It is positive and nonzero on the open square and has logarithmic
singularities at ty=-3 and tx=3. Its elementary location/decay is already
an old-alphabet phenomenon. The other singularities of this example
include ty=5 and tx=-5; this is not a full psi40 singular-locus claim.

The supplied cheb_decay.out has no coefficient table, tail enclosure,
continued physical log coefficient or elliptic differential identity.
Its reported slopes are consistent with the old logarithmic calibration.
Finite slopes neither prove exact complex singularity locations nor
certify holomorphy on the whole square. Even a proven far-grid logarithm
must have its tangential coefficient tested; a nonzero logarithm alone
proves nothing new, since psi20/psi30 already show the same decay.

## 8. Reproduction, scope and next decisive evidence

Run serially from the task root, with bytecode disabled:

    PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 venv/bin/python codex/phase13/x25_exact.py
    PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 venv/bin/python codex/phase13/x25_periods.py
    PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase13/x25_audit.py

The exact checker verifies 20 assertions; the numerical checker verifies
nine, including two independent normal seed vectors, recurrence checks,
refinement and Legendre comparisons. Peak scientific RSS is below 70 MB.
They take seconds serially. No reconstruction, energy calculation,
large elimination, installation, agent, repository edit or commit was
performed. The attempt is shorter than the allowed 1–3-hour ceiling
because it reaches an explicit formula and a concrete converse obstruction.

Read-only prior context: COLLAB phases 11–13, THEORY phases 10–11,
x20/x21/x22/x24 reports; x23 exists as the COLLAB audit entry, not a report.
Actual old-atom definitions were checked in x19_atoms.jl, the x16 report
and phase8/report.md; the operator was checked against its binary file.
The resource/ownership audit preserves shared-log prefixes and scientific
input hashes. It is an application-level audit, not a syscall trace.

**Next positive certificate:** give the actual continued psi40 sector,
based path, generic arc, normal/log convention, coefficient with a
nonzero enclosure, and an exact or validated identification of its
tangential H1(E_q) factor, allowing the gauge/differential transforms of
(2.5). That would prove outside-C3 survival without completing the old
four-parameter relative-chain pairing. A fit to two K functions at a few
points is not an identification theorem. A negative result on those
particular sheets or only at normal order zero is not a converse theorem.

**Summary log line:** PROVED/VERIFIED type-II residue fingerprint with
quartic-to-Legendre 2-isogeny, explicit tangential PF and exact `(2,2,1)`
operator monodromy; full stated C3 has only finite/Tate generic-circle
nearby constituents. CO-SIGN the sufficient outside-C3 test, DISSENT from
the raw-limit converse (delta Pi counterexample); use Im N and normalized
normal coefficients for quotient/middle factors. Physical-edge analyticity
is supported by its local theorem; enlarged-square analyticity and an
actual physical elliptic log coefficient remain OPEN, with no claude11
notes/claim certificate available. Conjecture A remains OPEN.
