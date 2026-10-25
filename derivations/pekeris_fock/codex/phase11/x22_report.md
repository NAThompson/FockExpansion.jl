# TASK_codex_22 — the source-independent pilot has no T5

2026-10-24. [X] denotes the present derivation; [C] denotes supplied Claude
results. Status labels apply separately to the statements below. All outputs
are x22_* in codex/phase11, apart from authorized append-only shared entries.

## 1. Result and task disposition

**PROVED, using the CITED mixed-Hodge inputs in §6; algebra VERIFIED:**
the source-independent part of the actual C3 nuclear-times-panel-3/4
contribution does **not** contain T5, including as a middle composition
factor. Here source-independent means the A(q)^2 part of the exact
fixed-endpoint split, A(q)=pi/2-atan(q), with the actual c0,c2, BOTH panels,
BOTH nuclear centres, and the pure order-four inverse. This is a
non-survival result for the fallback explicitly requested in TASK_codex_22.
It is not a vanishing of the function and not an old-elliptic-alphabet
membership theorem.

**OPEN:** a certified framed T5 cut of the full physical charge-sector
integral. Step 1 of x21 §5 was attempted as detailed in §2 but not completed;
no internal density residue is relabelled as a physical target variation.
Consequently steps 2–3 for that full-density cut have no input certificate.
The authorized simpler model is resolved instead.

The decisive cancellation was absent from the bare-measure analysis in
x19/x20: the angular-square-root term in the regular resonant zonal inverse
has a factor sqrt(1-|r|^2). Multiplying by the Feynman beta measure cancels
that radical. The logarithm still uses both square roots, so cancellation
of the prefactor alone would not prove the claim. An explicit filtered
four-state connection below proves that its graded characters are the
separate conic characters, never their K3 product. Their relative
two-dimensional direct images have only Tate and curve graded factors;
the rank-five K3 Hodge structure has Hodge level two and cannot occur.

**OPEN:** T5 survival in the complementary B(q,b)^2 part or in full psi40.
Since the A part has no T5, it cannot cancel a T5 constituent that is
actually certified in the B part. Other terms of the full charge sector
still have to be included. No positive or negative answer to Conjecture A
for physical psi40 is asserted.

## 2. Attempt at the full physical target cut

Let the sphere SOURCE be y=(y1,y2,y3,y4), target x=(u,v,w,0),
w^2=1-u^2-v^2, and T=(Lambda^2-32)^+ the pure inverse. These are angular
coordinates, not the Pekeris distances. Extract the Z^2 coefficient of the
two order-three nuclear panels. Their exact data are

    D=1+q^4,
    P3=(-2q^2,1-q^4),       P4=(2q^2,1-q^4),
    Wj^2=D-Pj.y,
    b=-y2,   hb^2=(1-y2)/(1+y2),
    c0=(1-y2^2)/18,
    c2,3=-2(1+y2)y1/9,     c2,4=+2(1+y2)y1/9.

The panels have the fixed interval 0<=q<=1, as proved in x16/x21. Define

    A(q)=pi/2-atan(q),             B(q,b)=atan(q/hb),
    Hj=c0-c2,j q^2,               Kj=-c0+c2,j hb^2/q^2.

Thus the ACTUAL panel coefficient is

    Ej(y)=(2/pi) integral_0^1 [Hj(y) A(q)^2+Kj(y) B(q,b)^2] dq/Wj.

In particular,

    H3=(1+y2)(1-y2+4q^2 y1)/18,
    H4=(1+y2)(1-y2-4q^2 y1)/18.

The selected full contribution is 2T[v1(E3+E4)], where
v1=sqrt(2)[(1-y1)^(-1/2)+(1+y1)^(-1/2)]. No nuclear centre is discarded.

### 2.1 The initial framing and candidate loops

**SPECIFIED, not transported:** use the real S3 Green chain, 0<=q<=1 and
0<=r,s<=q, with the fixed physical roots and the Green kernel of x18. Replace
B^2 by its exact rational-density period

    (1-b^2) integral_[0,q]^2 dr ds/(Lr Ls),
    Lr=1+b+(1-b)r^2,       Ls=1+b+(1-b)s^2.

For panel 3, reflect coordinates by x' = (-v,-u,w,0). This sends P3/D to
the phase-10 standard p(q) and the nuclear centres to +-e2. Use
x'2=1/3, x'1=lambda, base lambda=1/4. Candidate target loops circle
lambda=+-sqrt(8/9), and subsequently lambda=+-1, along paths avoiding the
other punctures; a local circle can have radius less than 1/100. Panel 4
is transported by electron exchange, including its centre signs. Giving
these target paths does not specify their lifts through the higher
relative direct image.

The missing calculation is the deformation of the angular Green chain
with its q,r,s boundaries. At fixed SOURCE y, the r,s poles +-i hb(y) do
not move with the target. Thus a target loop does not by definition select
one of the four internal density tori. A nontrivial such torus can arise
only after the Green/angular chain is deformed and its boundary
intersections are computed. Those intersection numbers and the resulting
boundary terms have not been certified here.

### 2.2 A concrete cancellation that the lift must respect

**PROVED / VERIFIED:** on the moving source arc y=Pj(q)/D,

    hb^2=q^4,  hb=q^2,  atan(q/hb)=atan(1/q)=A(q),  0<q<1,
    Hj+Kj=0.

Both actual weighted moment numerators therefore vanish at this Wj=0
incidence. The A and B pieces cancel there before taking a local leading
singular cut. This does not prove that every higher local term, or a full
relative cut, cancels. It does show why a cut of only A or an arbitrary
internal residue cannot be substituted for the full physical lift.

**STOP / OPEN for step 1:** no chain map from a candidate target loop to
a specific density/Green cut was obtained. In particular, no five-state
vector or nonzero Betti pairing for that full cut is reported. Proceeding
to a modular reduction of an arbitrary internal numerator would answer a
different question. The rest of this report resolves the permitted A model.

## 3. Exact fallback, including the actual numerator

Define Ej,A by retaining only Hj A(q)^2 in Ej, on the SAME fixed interval
[0,1]. This is the precise fallback, rather than a split performed on
the obsolete source-dependent min(1,hb) interval. Put

    C3,A = 2T[v1(E3,A+E4,A)].

For each j=3,4 and a=+-e1 let p=Pj/D and r=t a+(1-t)p. The beta identity
gives the exact representation

    C3,A = (4 sqrt(2)/pi^2) sum_(j,a)
        integral_0^1 dq A(q)^2/sqrt(D)
        integral_0^1 dt/sqrt(t(1-t))
        T[Hj(y)/(1-r.y)].                         (1)

The four terms in this sum include both panels and nuclear centres.
Polynomial Hj has degree at most two in SOURCE y. If Uinv, Ulog, Uxlog
denote the regular zonal inverses of 1/ell, log ell, ell log ell,
ell=1-r.x, then its numerator is handled by

    T[y_i/ell]=-partial_(r_i) Ulog,
    T[y_i y_j/ell]=partial_(r_i) partial_(r_j) Uxlog.

The pure projection commutes with these parameter derivatives. Fixing the
regular zonal solution instead changes only a target degree-two harmonic;
the pure inverse subtracts such a polynomial. It cannot introduce T5.

Set

    n=r.r, z=r.x, Q=n-z^2, h=sqrt(1-n), rho=sqrt(n),
    gamma=acos(z/rho),
    theta=2 atan[(1+rho)tan(gamma/2)/h],
    H=gamma-theta,         J=H/sqrt(Q),
    R3=z(4z^2-3n).

The formulas use compatible real branches for 0<n<1, |z|<sqrt(n), then
analytic continuation. Singular-looking n=0 coordinates are interpreted by
the actual zonal solution, or by another chart.

### 3.1 The cancellation in the actual regular zonal inverse

**PROVED / VERIFIED:** modulo a polynomial of degree at most TWO in target
z, each inverse is

    Uk = Ak(n,z) log(1-z) + h Ck(n,z) J.          (2)

The coefficients are

    Ainv=(4z^2-n)(3n-4)/(12n^3),
    Alog=-(3n^2 z^2+3n^2-12n z^2-2n+8z^2)/(24n^3),
    Axlog=-(-6n^3 z+15n^2 z^2+5n^2-20n z^2-2n+8z^2)/(120n^3),

    Cinv=-(n-4)R3/(12n^3),
    Clog=(1-n)R3/(12n^3),
    Cxlog=(1-n)^2 R3/(60n^3).

These are identities for the actual regular-antipodal construction, not
merely solutions of a necessary PDE. The checker reconstructs the
variation-of-parameters primitives from independent rational integration
recurrences. For U=W/sin(gamma), regularity at the antipode sets
kappa=-2 I2(pi)/pi. Its gamma/2 term cancels the rational gamma coefficient
in I2, leaving exactly the displayed h times H. The checker verifies all
three complete identities modulo explicit target polynomials of degree
two; their coefficients are saved in x22_verify.json. It also checks the
nonpolynomial terms of the angular PDE independently.

The polynomial coefficients can contain parameter-only logs or radicals.
After the parameter integrations they are still target polynomials.
When individual harmonic normalizations need endpoint regularization,
keep the expressions together or use the pure projected normalization;
these changes remain in that polynomial space. Boundary terms from
relative integration are retained in the argument in §5.

### 3.2 Why the remaining logarithm does not restore the K3 character

Prefactor cancellation alone is insufficient: theta still involves h.
The following full connection resolves this issue. Direct differentiation
on compatible branches gives

    d_z J = z J/Q - 1/Q + h/[(1-z)Q],
    d_n J = -J/(2Q) + z/(2nQ) - (z+n)h/[2n(1-n)Q],
    d_z h=0,             d_n h=-h/[2(1-n)].      (3)

The checker verifies its full rational curvature and compares BOTH
derivatives with the actual inverse-trigonometric function at independent
interior points. Add L=h log(1-z):

    d_z L=-h/(1-z),       d_n L=-L/[2(1-n)].

Thus (J,L,1,h) is a four-state rational connection. Its composition
characters, apart from Tate shifts from logarithms, are

    chi_Q : 1/sqrt(Q),   chi_h : sqrt(1-n),   1,   chi_h.

In particular chi_Q*chi_h does NOT occur in its semisimplification.
This statement concerns all graded constituents, rather than only a
Hom test at one endpoint of an extension.

Let mu=1/sqrt(D t(1-t)). For every actual panel/centre pair in (1),

    1-n = 2N t(1-t)/D,       N=D-a.Pj,
    N=(1+q^2)^2 or (1-q^2)^2,
    mu*h=sqrt(2N)/D.                            (4)

On 0<q<1 this last function is rational in q up to sqrt(2). Equation (2)
times mu consequently belongs to the connection (3) plus L, with rational
coefficients: mu log(1-z)=(mu*h)L/(1-n). Parameter differentiation for the
actual quadratic Hj stays in the same filtered connection. Differentiating
mu*h before imposing r=t a+(1-t)p only adds rational factors such as
1/(1-n); it does not create chi_Q*chi_h.

Finally A(q)^2 is a length-two logarithmic period on
P1\{i,-i,infinity}; dA=-dq/(1+q^2). Tensoring with (1,A,A^2) adds only
Tate extensions. A sufficient connection for the target-nonpolynomial
INTEGRAND therefore has at most 4*3=12 states per panel/centre, with
explicit matrices in x22_verify.json. This is not a rank bound or minimal
rank computation for the integrated target function.

## 4. The separate covers are rational; their product is the auxiliary K3

**PROVED / VERIFIED:** the chi_h cover is rational over the q line. Write
k(q)=2N/D. Its equation is h^2=k t(1-t), and a slope parameter a gives

    t=k/(a^2+k),             h=ka/(a^2+k).        (5)

The chi_Q cover is also rational after the finite TARGET extension
adjoining i and w=sqrt(1-u^2-v^2). Since a and p are unit vectors in the
same coordinate plane,

    disc_t Q=-4w^2 (a1 p2-a2 p1)^2.

The roots tplus,tminus of Q are therefore rational in q over that target
field. For Q=a2(t-tplus)(t-tminus), a rational parametrization is

    t=(a^2 tplus-a2 tminus)/(a^2-a2),
    sqrt(Q)=a a2(tplus-tminus)/(a^2-a2).          (6)

Degenerate parameter values use other charts or resolution. Generic
rationality of the smooth compactifications is all that is required.
Both (5) and (6), and the discriminants for all four physical pairs, are
checked exactly.

The PRODUCT character chi_h*chi_Q, in contrast, has cover

    Y0^2=D t(1-t) [D^2 Q],

up to rational squares and the constant factor 2N. This is precisely the
bare K3 of x19/x20 after their normalization. The normalized *bare*
measure can indeed have its T5 periods. Equations (2)–(4) show why that
product character is absent from the ACTUAL weighted A-model connection.
There is no contradiction with x20's Picard-number or monodromy theorem.
The present argument does not dispute those bare-surface results.

## 5. A framed exclusion, not just a vanishing one-cut test

Here is the Hodge-theoretic step that upgrades the exact connection to a
non-survival theorem. Work on a regular target locus and, if necessary,
make the finite target extension used in (6). The connected SO5 system
T5 from x20 remains irreducible under a finite cover and has pure Hodge
numbers (1,3,1), hence Hodge level two. Integer Tate twists leave that
level unchanged.

**CITED:** relative algebraic cohomology has a compatible mixed Hodge
structure; morphisms are strict; the logarithmic residue/weight spectral
sequence computes an open surface using its compactification and boundary
strata. The fixed-part and semisimplicity theorems apply to the pure
geometric variations. Exact locators and primary links are in §6.

**PROVED application.** Resolve the pole and chain-boundary divisors on
each of the two rational surfaces in §4. For their relative H2 objects:

* Compactification H2 is algebraic/Tate, since a smooth rational surface
  has h20=0 and all its H2 is generated by divisor classes.
* Non-Tate boundary or pole terms come from H1 of normalizations of
  curves, with integer Tate shifts. They have Hodge level one.
* Points, incidence kernels and curve H2 contribute only Tate factors.

The same assertions cover the trivial character on the rational q,t
surface. Equation (3) supplies an explicit differential-module filtration
with only these rank-one characters. Tensoring by (1,A,A^2) repeats the
same characters in an extension. This computational filtration need NOT
be identified with the weight filtration.

Apply the relative direct image to this finite filtration. Its long exact
cohomology sequences imply that every simple constituent in degree two
comes from a degree-two relative direct image of a graded character.
The latter are direct summands of relative cohomology on the rational
double covers, so their pure graded factors are Tate or curve H1 twists.

One further point is needed to exclude the underlying DIFFERENTIAL
constituent, rather than just a specified Hodge realization. A curve H1
variation cannot have the underlying irreducible T5 local system with
some other Hodge grading. Indeed let V=T5 and suppose V occurred in such
a curve variation W. By the geometric fixed-part theorem,
H0(Hom(V,W)) is a constant Hodge structure, and its evaluation with V is
a Hodge morphism. Absolute irreducibility and End(V)=Q make evaluation
identify V tensor this nonzero multiplicity space with its isotypic
component in W. Tensoring a nonzero constant Hodge structure with V has
Hodge level at least TWO: its extreme Hodge types differ by at least the
two steps already present in V. This is incompatible with a subvariation
of curve H1, whose level is ONE. The same argument applies after Tate
twists and finite target covers. Pure semisimplicity and the weight
filtration then rule out V in every relative rational-surface object.

Consequently T5 is excluded in ALL composition positions of the A-model
direct image, including the middle of a mixed extension. This argument
does not assume that an arbitrary D-linear map preserves Hodge filtrations.

The actual rectangular chain and its boundary maps define a framed
subquotient of this same relative direct image. A framing can remove
constituents, not create the excluded T5. This proves the result for
each panel/centre individually; their finite sum with the exact Hj and
the pure harmonic subtraction also has no T5. A catalogue of boundary
curves, or proof that they are only the OLD elliptic curves, is unnecessary
for this particular exclusion and is NOT supplied. Some curves could
have larger genus. Their H1 still has level one.

This is stronger than finding a single zero cut and different from
inferring non-survival merely from a rational ambient variety. The
filtered connection must first eliminate the PRODUCT character; without
that step the K3 cover really is present in the ambient presentation.

## 6. Sources, attribution and status of inputs

**CITED / inspected online 2026-10-24:**

1. Deligne, *Théorie de Hodge II*, Thm. 2.3.5, §§3.1–3.2, especially
   Cor. 3.2.17, and Cor. 4.1.2, Thm. 4.2.6, Cor. 4.2.8: strictness,
   residue/weight spectral sequence, compactification image, geometric
   fixed part and semisimplicity.
   [Primary text](https://www.numdam.org/item/PMIHES_1971__40__5_0.pdf).
2. Deligne, *Théorie de Hodge III*, §6.3, Ex. 8.3.8 and Prop. 8.3.9:
   relative cohomology, mixed Hodge structures and their exact sequences.
   [Primary text](https://www.numdam.org/item/PMIHES_1974__44__5_0.pdf).

The rational-surface and filtered-connection application is the present
argument, not a statement quoted from these papers. The rank-five T5
Hodge numbers, irreducibility and SO5 monodromy are the proved/cited inputs
of x20, not a new recomputation of the PF operator here. Claude's newest
four-edge common-invariant-space computation is not used or co-signed.

Local read-only inputs include COLLAB.md, TASK_codex_22.md, x21 §5,
phase3/structure.md, x18/x19/x20 reports, and the repository's
src/feynman_k4.jl and derivations/k4/README.md. The latter are attributed
sources for the regular zonal construction. No Claude/repository executable
was imported or run. The checker reconstructs the relevant primitives
independently by elementary integration recurrences.

## 7. What remains for the full density, and cost

**OPEN:** transport the INITIAL physical chain in §2, compute the target
nearby-cycle/boundary maps, reduce the resulting CLOSED form to T5, and
pair it with the actual transported cycle. The double residues in x21
remain useful once that map is known; they do not determine it.

The B piece cannot use (1) unchanged: b=-y2 is an integrated source
coordinate. With r,s retained it has two extra linear source denominators.
Combining the two half-power Coulomb denominators and these two simple
denominators needs a three-dimensional Dirichlet simplex; together with
q,r,s this gives a sufficient SIX-parameter zonal representation. Equivalently
use the beta parameter t and a further two-dimensional simplex. The density
factor (1-b^2) and actual c0,c2 give polynomial numerators of degree at most
four after clearing the displayed q-only powers. Target/parameter derivatives
are still finite. This is a representation and parameter count, not a minimal
rank or a cost bound for its direct image.

The specific rational-surface exclusion in §5 no longer applies to that
six-dimensional geometry. Nor can the A and B local source-arc cancellation
be promoted to a full T5 cancellation. Full Conjecture A remains **OPEN**.

**Engineering estimates, not benchmarks of an unimplemented algorithm:**

* Reproduce this A-model algebra/connection and quadratures: seconds,
  about 70 MB measured. No modular elimination is needed for its exclusion.
* Once a physical cut map is supplied, a five-state one-point/one-prime
  reduction should remain a seconds-to-minutes pilot. Budget 10 minutes,
  512 MB, and no more than about 1,500 columns/5,000 rows; stop instead of
  escalating blindly. A reconstructed identity or validated nonzero cycle
  pairing is still required after a modular nonzero screen.
* A first local full-density chain/nearby-cycle attempt can be allotted
  1–3 hours of serial development and small experiments under the same
  memory cap. There is no basis to promise a certificate in that time:
  the missing incidence/transport map is mathematical work, not an already
  available matrix operation. Exact reconstruction/pairing cost is unknown
  until its size is known. A full negative calculation over all relevant
  graded pieces has no justified finite engineering bound here.

**CONJ:** the complementary B geometry may carry a surviving new factor.
Neither its necessity nor its absence is implied by this result. Even a
negative T5 decision for B would not automatically prove membership in the
old alphabet or decide every possible new constituent.

## 8. Verification and scope

Run from this directory, ONE command at a time:

    PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
      TMPDIR="$PWD/codex/phase11" venv/bin/python codex/phase11/x22_verify.py
    PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase11/x22_audit.py

The verifier checks exact full connections, the actual regular zonal
identities modulo explicit polynomials, all four panel/centre beta
cancellations and split conic covers, numerator parameter identities,
source-arc cancellation, independent actual-angle derivatives, the beta
integral and the fixed-endpoint A/B decomposition. The Hodge-theoretic
exclusion is a proof in §5, rather than an assertion that numeric tests
certify a motive.

Logs retain the initial 40-check run, a SymPy-integer JSON export failure
in the expanded checker, and subsequent successful runs. The separate final
streams and JSON/audit give the final evidence. An abandoned fully symbolic
rational-primitive exploration was interrupted before any result; it is not
used as mathematical evidence. Input hashing briefly overlapped that
abandoned exploration; there were no concurrent scientific computations.
No agents, installations, commits, energy solves or stencil builds were used.
All measured verifier runs stayed around 70 MB, far below 2 GB.

No command wrote a Claude-owned file or any path outside this directory.
Shared entries are appended only. The audit records observed concurrent
input/output changes separately and checks all used scientific inputs and
shared prefixes. It is a hash/command-scope audit, not an OS syscall trace.
