# TASK_codex_21 — physical T5 survival and two exclusion tests

2026-10-24. Credit: new arguments/checks below [X]; the displayed scalar L5
and spectral data are Claude's read-only inputs [C]. Status labels apply to
individual statements. All calculations were small and serial.

## 1. Result

**OPEN: Conjecture A for the physical C1, C2 and C3 sectors.** Neither a
nonzero physical T5 projection nor its vanishing has been certified. The
complete weight-graded relative Gauss–Manin object with the actual density
and chain has not been computed. This report uses the fallback explicitly
allowed by task item 3, and does not call the survival question decided.

There are, however, substantive new conclusions:

* **PROVED / VERIFIED:** the closed K3 periods satisfy `(Delta_S3-1)Pi=0`.
  An explicit rational operator of order three takes this very T5 module
  to nonzero solutions of `(Lambda^2-32)m=0`. Consequently the proposed
  necessary P-harmonic-section test PASSES; it cannot exclude T5.
* **PROVED / VERIFIED:** Claude's L5 has a three-dimensional space of
  holomorphic germs at each collinear endpoint on x2=1/3. There is even a
  nonzero solution analytic at BOTH physical endpoints of that line. Thus
  analyticity of a selected physical branch is not an exclusion invariant.
* **PROVED:** the bare relative K3 surface retains T5 in weight two, despite
  algebraic boundary divisors. **PROVED / VERIFIED:** the two density
  parameters have explicit double residues -1/4 or +1/4; their actual panel
  numerators do not vanish identically or become exchange-odd. These are
  ambient classes/cuts, not a nonzero pairing with the physical target chain.
* **PROVED:** Claude's cyclic quotient formula is correct. Testing only a
  quotient or only a submodule of Q does not test every Jordan–Hölder
  constituent when Q is a mixed extension. Pm=0 alone also omits Ann(psi40).

## 2. Actual period and the weight calculation that is justified

Use sphere source y, target x=(u,v,w,0), w^2=1-u^2-v^2. These u,v are angular
coordinates, not Pekeris distances. Put

    D=1+q^4, Pq=(q^4-1,2q^2), p=Pq/D,
    a=e2, A=1-v^2, B=a.Pq-v(Pq.x), C=D^2-(Pq.x)^2,
    S_x: Y^2=X^3+2B X^2+ACX,  omega=dq wedge dX/Y.

The other T5-candidate panels use the reflection specified in COLLAB; the
same argument transports to them. The bare surface has rank T=5 and Hodge
numbers (h20,h11,h02)=(1,3,1), by x20's rho=17 certificate.

### 2.1 Bare surface with boundaries

**CITED inputs:** the logarithmic residue/weight spectral sequence and
strictness of morphisms of mixed Hodge structures are in
[Deligne II, §§2.3, 3.1–3.2](https://www.numdam.org/item/PMIHES_1971__40__5_0.pdf).
Relative cohomology and its exact sequence carry compatible mixed Hodge
structures by [Deligne III, 8.3.8–8.3.9](https://www.numdam.org/item/PMIHES_1974__44__5_0.pdf).

**PROVED application.** Resolve the surface and its algebraic pole/boundary
divisors. For a normal-crossing pole divisor E on S, let U=S\E, and let
E^[j] be the disjoint union of its j-fold intersections. The terms relevant
to H2(U), after taking the indicated Gysin kernels/quotients, are

| Weight | Origin | Consequence |
|---:|---|---|
| 2 | H2(S) / span of classes of components of E | T5 survives; the removed classes are algebraic |
| 3 | H1(E^[1])(-1) | curve factors, depending on the actual divisors |
| 4 | kernel of the incidence/Gysin map from H0(E^[2])(-2) | Tate factors |

For a boundary curve Bnd contained in U, the exact sequence

    H1(Bnd) -> H2(U,Bnd) -> H2(U) -> H2(Bnd)

shows that T5 lifts as a subquotient of Gr2^W H2(U,Bnd): H2(Bnd) has only
Tate weight-two classes and receives no nonzero Hodge morphism from the
simple transcendental T5. Removing E killed only divisor classes in H2(S),
so it could not kill T5 first. Boundary H1 contributes weights 0,1,2 when
the boundary is a possibly open, reducible curve. This proves occurrence
in the **ambient bare relative surface**, with no assertion about a scalar
physical period. It also explains why “an open chain has only boundary
weights, hence no K3” would be wrong.

### 2.2 The density prevents identifying that object with the full answer

The actual moments are, with b a SOURCE coordinate and h^2=(1+b)/(1-b),

    N0=(2/pi) int_0^1 [Aq^2-Bq^2] dq/W(q),
    N2=(2/pi) int_0^1 [h^2 Bq^2/q^2-q^2 Aq^2] dq/W(q),
    Aq=pi/2-atan(q), Bq=atan(q/h).

The fixed endpoint 1 is the exact phase-7 presentation, not an approximation
to min(1,h). Inside the angular Green inverse,

    Bq^2=(1-b^2) int_[0,q]^2 dr ds/[Lr Ls],
    Lz=1+b+(1-b)z^2.

Thus r,s poles move with the integrated source y. They are not two punctured
lines in a product over the already integrated S_x. In particular, taking
the tensor product of H2(S_x) with a target-only logarithm module would
repeat the density-pull-through error identified in x19.

**PROVED fibre calculation.** At fixed b, the logarithm is a relative H1
period on `(P1\{+ih,-ih},{0,q})`. Its two graded factors are Q(0), Q(-1).
The symmetric square realizing {1,Bq,Bq^2} has graded factors Q(0), Q(-1),
Q(-2), of weights 0,2,4. If these factors were external to the K3 direct
image, candidate factors would be T5, T5(-1), T5(-2), of weights 2,4,6.
Here that last tensor-product description is only a **candidate**: the
angular pushforward and boundary maps still have to be performed. Green
logarithms supply further extensions. Hodge weight counts twice the usual
logarithm length for these Tate factors.

For the actual full period one must resolve its algebraic cover with its
pole divisor E and all chain boundaries Bnd, and compute the hypercohomology
of the relative logarithmic mapping cone (using a simplicial resolution of
Bnd when needed). This supplies an exact definition of the required Gr^W,
not its computation. The r,s pole incidences with angular, q and Green
boundaries, and their residue/Gysin maps, are missing. **OPEN:** their ranks,
the full Gr^W list/multiplicities, and the T5 part of the physical framed
subquotient. No Hodge-number table for that larger object is invented here.

The distinction between ambient cohomology and a class paired with a
specified chain is the same general distinction made for relative periods
in [Brown, §0.2 and §9.1](https://arxiv.org/pdf/1512.06409). His graph-specific
theorems are not being applied to this non-graph integral.

### 2.3 A computable density residue, including the actual numerators

**PROVED / VERIFIED.** At r=epsilon i h, s=delta i h,

    Res_r Res_s [(1-b^2) dr ds/(Lr Ls)] = -epsilon delta/4.

The four entries are (-1/4,+1/4,+1/4,-1/4). The two minus signs belong to
equal pole signs. These are explicit nonzero functions (constants here),
not a numerical fit. A positively oriented small two-torus at an equal-sign
pair integrates this density to `(2 pi i)^2(-1/4)=pi^2`. The sum of all four
residues is zero, as it should be; the four torus classes must not be
identified without a homology relation for the actual chain.

For the actual output `c0 N0+c2 N2`, its same-sign density residue is

    K(q,b,y)= [c0-c2 h^2/q^2]/4,                 (*)

apart from the common 2/pi, W, Coulomb and Green factors. With y denoting the
panel parameter used in phase3/structure.md, panels 1/2 give

    [Z]K = sqrt(2)(b+1)[2b q^2-9b-2q^2-16y^2+8]/(576q^2),
    [Z^2]K = -sqrt(2)[b^2q^2-20b y^2+10b-2q^2-20y^2+10]/(256q^2).

At b=0, y=-1/2, q=1/2 these are `7sqrt(2)/288` and `-9sqrt(2)/128`.
The paired panels coincide at this electron-exchange-fixed source, so an
exchange-odd-integrand argument cannot cancel them. Panels 3/4 give

    [Z^2]K=-(b+1)[bq^2-q^2-8y^2+4]/(72q^2).

At equal electron radii, y^2=1/2 for both panels, and this becomes
`(1-b^2)/72`, again nonzero away from b=+-1. Multiplication by -2v0 or
+2v1 changes charge and normalization as in x18; it does not make these
polynomial coefficients identically zero. These point substitutions are
exact tests of an alleged integrand parity, not a period calculation.

**LIMIT:** (*) is not yet the requested transcendental-cycle pairing.
It is a residue in INTERNAL density variables. To infer physical survival
one must identify it as the image of a physical TARGET continuation/cut,
reduce the remaining form to the T5 graded factor, and pair it with the
corresponding transported chain. None of those operations is implied by
the nonzero -1/4. A residue class can be nonzero in an ambient object and
annihilated by the chosen Betti framing or the subsequent direct image.

### 2.4 What the symmetries actually cancel

**PROVED:** the involution q->h(y)/q gives
`g(h/q)h/q^2=-g(q)` for each physical moment. It cancels the segment [1,h],
which is why the fixed endpoint works. The physical interval
`[0,min(1,h)]` is not invariant under it. At h=1 the cancelled segment is
empty while N0 is positive; an independent 45-digit quadrature with
W^2=1-q^2+q^4 gives N0=0.72983851865300320819285971304734786139...
Furthermore h(y) depends on the source, so this involution cannot simply
be declared an automorphism of the target-only K3 after angular integration.

Electron exchange transports panels and the target base. A charge sector
is preserved; it is not the operation of summing unrelated powers of Z.
No symmetry demonstrated here forces a zero T5 pairing. Conversely, the
failure of these proposed symmetry cancellations does not prove that a
cohomological or full charge-sector cancellation is absent. A nonzero
symbolic charge coefficient would establish a sector claim; survival at
Z=2 would additionally need its full charge sum checked.

## 3. Audit of Q and the P-harmonic test

Work over an algebraic angular field K with derivations du,dv and constants
C. Let D=K<du,dv>, f one actual psi40 charge sector, I=Ann_D f, s=Pf,
M=Df=D/I. The full hierarchy source includes

    s=12 psi41+2 psi42-2V psi30+2E psi20,

with pure-inverse projection if required. Harmonic subtractions are old
algebraic factors. Thus Ds is in the old category, under the already
specified old-module/holonomy hypotheses.

**PROVED algebra, no finite-rank assumption needed:** set
`J={R in D: RP in I}=Ann_D s`. Right multiplication by P induces

    0 -> D/J -> D/I -> D/(I+DP) -> 0,
           [R] |-> [RP].

So Claude's Q formula and constituent equivalence are correct, provided
T5 is absent from Ds. The equivalence is about a **composition factor**.
It does not say that Q is a direct sum or that T5 is a quotient of Q.

**PROVED on a regular connection locus:** applying the local solution
functor gives

    Sol(Q)={h in Sol(M): Ph=0};   h |-> Ph is the map to Sol(Ds).

This is solution duality, not an assertion that the scalar differential
operator P annihilates every element of the differential module Q. For
example it annihilates its cyclic generator, but need not annihilate its
derivatives. Identifying Sol(M) with the span of branches of f also uses
regular-singular minimality/cyclicity. On the geometric regular-singular
period realization it has the intended meaning; that identification is
not a general rule for arbitrary irregular holonomic functions.

If g fixes the source branch, `(g-1)f` is indeed P-harmonic. The whole
kernel can also involve linear relations among several continued sources.
One must not replace the kernel by only single loops that fix s.

**PROVED refinement of the quotient test:** a surjection Q->T5 would send
the cyclic generator to a nonzero m in T5 with

    Pm=0 AND I m=0.

The first is necessary; it is not sufficient for passage through Q.
Testing a submodule requires Hom_D(T5,Q), equivalently a quotient test
after differential duality, including the adjoint/gauge conventions.
Both endpoint Hom tests can miss a T5 in the middle of a nonsplit mixed
extension. For example a module with nonsplit successive layers
Tate, T5, Tate can have T5 as a constituent without T5 at either endpoint.
The weight filtration's pure graded pieces, rather than semisimplicity of
the entire relative object, are the appropriate setting for that case.

### 3.1 New explicit P-harmonic rational section of T5

**PROVED / VERIFIED.** For functions depending only on u,v on S3,

    Delta=(1-u^2)du^2-2uv du dv+(1-v^2)dv^2-3u du-3v dv,
    P=Lambda^2-32=-4(Delta+8).

In the original q,t coordinates of x20 put
`r=t a+(1-t)p(q)`, `Q0=r.r-(r.x)^2`. The holomorphic measure, up to
target-independent factors, is `dq dt/sqrt(Q0)`. For any fixed r,

    (Delta-1) Q0^(-1/2)=0.                       (1)

One direct proof sets z=(r.x)/sqrt(r.r) and computes
`(1-z^2)(1-z^2)^(-1/2)''-3z(1-z^2)^(-1/2)'=(1-z^2)^(-1/2)`.
This is also an exact rational identity checked by x21_verify.py.
Integration over a transported CLOSED K3 cycle therefore gives
`(Delta-1)Pi=0`. Differentiating an open chain would add boundary terms;
no such assertion is made for the physical weighted integral.

For a coordinate i define the tangential gradient Ai=grad_S3(x_i).grad and

    Ti(k)=Ai+k x_i,       Lk=Delta+k(k+2).

Homogeneously extend f by `F(z)=|z|^k f(z/|z|)` on four-dimensional ambient
space. Commutation of the ordinary Laplacian with partial_i gives

    L_(k-1) Ti(k) = Ti(k-2) Lk.                 (2)

For i=v, `Tv(k)=-uv du+(1-v^2)dv+k v`, over K. Combining (1),(2),

    Rv=Tv(-3) Tv(-2) Tv(-1),
    (Delta+8) Rv=Tv(-5) Tv(-4) Tv(-3)(Delta-1),
    P(Rv Pi)=0.                                (3)

The checker expands BOTH sides as rational differential operators and
tests every coefficient; it does not merely apply them to one sample.

The map is nonzero on T5. At u=1, generic v!=0, there is a period branch
`Pi=(u-1)^(1/2)c(v)+higher terms`, c nonzero, by the L5 half exponent.
The leading normal part of Tv is `-v du`. Hence

    Rv Pi=-(3v^3/8)c(v)(u-1)^(-5/2)+higher terms,

which cannot vanish. The SO5 irreducibility established in x20 means a
nonzero m=Rv Pi generates the same rank-five module. Thus D/DP has a
T5 quotient, explicitly. All coefficients of Rv are saved in
`x21_verify.json`. This settles the proposed necessary-section question
positively, without reconstructing a two-variable 5x5 connection.

It does **not** supply `I m=0`, a map Q->T5, or a physical period. The
P-harmonic functions in (3) are multivalued meromorphic functions, not
globally smooth physical spherical harmonics. Compact-sphere uniqueness
cannot be applied to them or to arbitrary continued branches of psi40.

## 4. Analyticity versus local monodromy

### 4.1 An exact counterexample from this same T5 operator

**PROVED / VERIFIED from the supplied L5:** at
`ell=+-sqrt(8/9)` its indicial polynomial, up to a nonzero constant, is

    rho^2(rho-1)^2(rho-2).

Repeated exponents alone do not determine logarithms. Here the analytic
solution count can be computed without that inference. Substitute
`f(ell)=g(z), z=ell^2`, divide Lf by ell, and expand about z0=8/9.
The coefficient that recursively determines a_(n+3) in `g=sum a_n t^n` is

    -(32768/9)(n+1)(n+2)^2(n+3)^2,   n>=0.

There are no negative-power constraints. Thus a0,a1,a2 are arbitrary and
give exactly three analytic germs. The rational recurrence verifies three
independent series through twelve coefficients; the displayed nonzero
pivot proves the statement at all orders. Regular-singular Frobenius
theory gives convergence. For example the seed (1,0,0) begins

    g=1+(8505/131072)t^3+(42023205/134217728)t^4+... .

The same dimension is three at the other endpoint. Transport both germ
spaces to an interior basepoint along the interval x2=1/3. Inside the
five-dimensional solution space their intersection has dimension at least
`3+3-5=1`. A nonzero solution is therefore analytic at BOTH physical
collinear endpoints on this line. Nevertheless it generates the full
irreducible T5 under global continuation. This directly disproves the
proposed inference from physical-edge analyticity on the available line.
It does not assert that an arbitrary such line solution is psi40, or solve
a global two-variable matching problem.

### 4.2 Which stronger assertion would be sufficient

**PROVED:** analyticity of one branch implies that the relevant local
inertia fixes that branch. It does not imply that it fixes its entire
global monodromy orbit. A branch obtained by going first around a different
singularity may acquire logarithms at the original edge.

If EVERY branch in the cyclic continuation space extended single-valuedly
across a generic collinear divisor, the local inertia would act as the
identity on that whole space, and hence on all its subquotients. This
would exclude T5, whose collinear local monodromy is nontrivial. A stronger
normal-closure formulation requires all conjugate inertias to fix the
cyclic vector. The supplied physical-branch analyticity/Chebyshev evidence
does not establish either stronger condition. No step from this local
counterexample to a theorem about all branches is being assumed.

### 4.3 Coordinate and numerical qualifications

In the uncentered coordinates actually used by the spectral evaluator,

    X=alpha-beta, Y=alpha+beta,
    u=cos((X+Y)/2), v=cos((Y-X)/2),
    1-u^2-v^2=-cos X cos Y.

A generic collinear edge is a simple divisor in X,Y, so logs there remain
logs. In contrast u=+-1 or v=+-1 pulls back QUADRATICALLY under the cosine
map. A half exponent there can become integral. Thus “unramified” is
incorrect if meant as a property of the map to u,v; if meant as absence of
branching of the pulled-back half power, that absence is precisely caused
by this ramification. At corners more than one vanishing factor must be
retained. Using Y centered by subtracting pi reverses the sign in the
cos-product identity.

**NUMERICAL / OPEN:** a shared empirical Chebyshev slope, including a
candidate ellipse radius 3+sqrt(8), probes the chosen real analytic branch.
It supplies neither the full complex singular locus nor monodromy of all
continuations. Finite slopes cannot certify exact singularity locations or
eliminate a pure constituent in a larger mixed system. No spectral data
were refitted or extrapolated in this task.

## 5. Cheapest remaining decisive computation

The necessary P-section search can now be skipped: (3) is an exact positive
answer to it. Another fit to real spec40 values would also not decide A.

**PROPOSAL, not a completed calculation:** compute one certified framed
T5 cut of an actual charge-sector integral. Start with C3's two nuclear
times panel-3/4 terms, or C1's two e-e times panel-1/2 terms, retaining their
actual c0,c2 coefficients. C3 is an attractive first pilot because its
numerators have no Z^1/Z^2 mixing before multiplying by the nuclear charge.

1. Keep the r,s rational-density representation, all q endpoints and the
   Green logarithm. Specify target loops on x2=1/3 (with a regular basepoint
   u=1/4), and transport the physical chain. Determine which combination of
   the four density torus residues and Green homogeneous cuts is actually
   produced. Record boundary terms. An arbitrary internal -1/4 torus is
   NOT a replacement for this step.
2. Reduce just that resulting closed two-form to the rank-five
   transcendental de Rham quotient of S_e2. Carry the boundary/Tate terms
   separately. One can use five closed-period classes and a small
   Griffiths–Dwork reduction, avoiding a cyclic closure of the whole mixed
   integral. The lowering operator supplies an explicit compatible
   P-harmonic T5 generator for checking the cut.
3. Assemble BOTH electron-related panels in the chosen charge sector;
   include any other term contributing to the same cut. Exhibit a nonzero
   component paired with the transported cycle, symbolically or with a
   rigorous enclosure. A nonzero D-linear T5 image of an actual variation
   proves survival by irreducibility. Internal residues alone do not.

**Cost estimate:** the present exact checks take about one second and
roughly 70 MB. Once the chain/cut map and five-state quotient are supplied,
the final vector arithmetic is negligible. A bounded finite-field pilot
for one numerator at one regular target should take seconds to minutes;
budget 10 minutes, work serially, and enforce 512 MB for that pilot. This
is an engineering estimate, not a benchmark of a missing reduction.
Keep unreduced elimination at no more than about 1,500 columns / 5,000 rows,
use sparse modular arithmetic, and abort rather than exceed the budget;
raw 64-bit dense storage at that size is 60 MB, but intermediate/storage
overhead still needs monitoring. One prime/point is only a screening step.
A positive identity needs a reconstructed exact certificate and an actual
cycle map; a controlled nonzero pairing can instead use validated numerics.

The difficult, currently uncosted part is step 1: identifying the physical
target variation inside the higher relative direct image. It may require
new residue/nearby-cycle work; no one-second or one-hour decision of A is
promised. If every tested cut vanishes, that is not a negative decision:
exclude T5 in ALL relevant weight-graded framed pieces, including possible
middle composition factors, to prove non-survival. A uniform complexity
bound for that complete negative calculation is **OPEN**. The task's
small-job limit is a reason to prioritize one positive cut, not to declare
the full problem inexpensive or solved.

## 6. Reproduction and scope

From this directory, run serially:

    PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
      TMPDIR="$PWD/codex/phase11" venv/bin/python codex/phase11/x21_verify.py
    PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase11/x21_audit.py

The checker verifies 25 scientific claims: complete lowering-operator
identities, the K3 measure equation, indicial polynomials, analytic-series
recurrences, density residues/numerators, a physical moment counterexample,
and coordinate pullbacks. Input hashes/shared prefixes and the final audit
are saved alongside it. No Claude checker/evaluator was executed or imported;
no claim of recomputing Claude's PF telescoper or irreducibility is made.
The prior x20 monodromy result supplies irreducibility where used above.
The baseline audit observed concurrent growth of Claude's jets_p1.txt,
an output not used for these mathematical checks. It is recorded separately;
the protected scientific inputs used here stayed unchanged. An initial audit
that required every hashed Claude output to be static correctly flagged this
concurrent write; the final audit checks scientific inputs and preserves the
full change record, rather than attributing another process's output to x21.

Only x21_* outputs in codex/phase11 and authorized append-only COLLAB.md /
THEORY.md Log entries are written. No agents, installations, commits, energy
solves, or recurrence expansions were used. Python bytecode is disabled.
All actual scientific jobs stayed far below 2 GB. The first optional OS
timing wrapper returned a sandbox sysctl error after the Python checker
had passed; the clean direct invocation is the verification evidence.
