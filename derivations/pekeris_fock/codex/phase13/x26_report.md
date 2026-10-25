# TASK 26 — phase 13 independent audit

2026-10-24. Status labels distinguish exact derivations, floating-point experiments,
and missing geometric identifications. All computations are serial; no Claude-owned
file is edited. The isolated corrected solver is `x26_chain_corrected.jl`.

## Verdict

**CO-SIGN, NUMERICALLY, with corrections:** the nonzero rational double variation
and termination in R3 survive removal of a real step-0 bookkeeping error. The
first variations in the original chain output are not all those of the stated
physical hierarchy. In particular its Z² psi41 first variation is doubled.
There is also an odd-variation sign convention error in both programs.

**DISSENT:** termination is not a test that separates equal from unequal T5
vanishing planes. Distinct planes with a common line give a nonzero double
variation and an identically zero triple variation, already on pure T5.
The proposed two-way alternative needs this third case.

**OPEN:** identification of the actual based strip pair on T5, and attribution
of the two rational coefficients to explicit original order-four C3 atoms.
Neither the chain arrays nor the saved line-A monodromies identify the required
transport through the u=-1 and v=1 divisors. No T5 survival or exclusion in
physical psi40 is claimed.

## A. Separable hierarchy, logarithmic equations, and code audit

Write u=cos(alpha), v=cos(beta), X=alpha-beta, Y=alpha+beta,
w=2 sin(alpha) sin(beta), Lambda²=-4 Delta_(S³). Direct change of variables gives

    w*(Lambda²-c) f = -16*(M_X-M_Y)f,
    M_z f = cos(z) f'' - sin(z) f' + (c/16) cos(z) f.

This is an exact identity, including the factor -16 and the potential sign.
The physical hierarchy source is

    h_kp=2(k+2)(p+1)psi_(k,p+1)+(p+1)(p+2)psi_(k,p+2)
         -2V psi_(k-1,p)+2E psi_(k-2,p).

Thus the modal right side is -w*h/16 as in the code. The transpose convention
for the Y eigenbasis, `F*iVY'` and reconstruction `modal*VY'`, is correct.
The repeated characteristic exponent is zero at each simple zero of cos X.

### Exact logarithmic sign correction

For the stated l=log(-cos X), l'=-tan X. With L=M_X-lambda,

    L(sum A_m*l^m) = sum_m l^m [L A_m
       -(m+1)*(2 sin X A'_(m+1)+cos X A_(m+1))
       +(m+2)(m+1)*sin²X/cos X*A_(m+2)].

Both programs instead implement the plus sign on the first cross term and use
+tan X for matching derivatives. They consistently solve the expansion in
**L=-log(-cos X)**. At Xm=+-pi, L=l=0 and both first derivatives vanish;
therefore the sign mismatch does not spoil their midpoint matching.
It reverses every odd normalized variation and preserves every even one.
A positive geometric meridian has log monodromy 2pi*i*d/dl; the code generator
is its negative. Zero/nonzero tests and Jordan degrees are unchanged.

### Endpoint conditions and matching

Assuming the A_m and modal sources are analytic, the only pole in the coefficient
equations is the displayed A_(m+2)/cos X term. Its pole forces A_j(X0)=0 for
j>=2. Its removable endpoint value is -sin(X0)*A'_j(X0), exactly the implemented
limit. These are necessary compatibility conditions, not an assumption that
all higher-log analytic coefficient functions vanish identically.

`sysmat` replaces the midpoint equation by the far/start endpoint value;
the regular-singular endpoint differential equation is retained. The unique
analytic homogeneous germ at a characteristic endpoint has one free value.
The two scalar constants in A0,A1 then match the ordinary value and derivative
at Xm. The two-element construction is consistent when the matching matrix
is invertible. Agreement on successively resolved elements is numerical
evidence for analytic continuation, not an exact analyticity theorem at
interior Coulomb/reflection loci.

### Real error in chain.jl: harmonic corrections counted twice

The physical loop first finishes the entire psi dictionary, including `harmc`.
Then step s=0 starts each coefficient from `axa'*psi[(k,p)]`, already corrected.
Later, at p-1, it again adds `harmc[(k,p)]` to VS/VF. For example psi21 is
already nonzero when the homogeneous continuation at (2,1) starts; its same
Y2 correction is replayed on entering (2,0). The same problem affects the
order-four log coefficients. `extend.jl` interleaves physical solve and
continuation and does not have this error.

The isolated corrected copy changes only the include path and disables this
replay. It keeps all final physical start values, harmonic data, sources and
matching conventions. Corrected n64/N48 and n80/N60 runs were performed at
256 bits, Z=0..4, E=0, three subsequent variation steps. Corrected step 0 is
compared directly with the saved independent `extend.jl` arrays.

The correction changes the Z² first variation of psi40 and halves that of
psi41; it preserves the reported double variation to the numerical resolution.
This preservation is a verified result, not a justification for leaving the
bookkeeping error in place.

### Commutation and harmonic additions

On the alpha,beta cover the coefficients w,wV and the harmonic polynomials
are single valued and analytic. With a formal log coordinate, differentiation
is d/dX=partial_X+L'(X)partial_L. Hence partial_L commutes with the hierarchy
operator, target derivatives, and multiplication by these coefficients.
Equivalently the local logarithm of monodromy commutes with the differential
hierarchy. A harmonic correction has zero variation, so no fresh correction
is to be imposed at subsequent steps. Variations of its forced descendants
are still part of the propagated hierarchy.

Passing only A1(X0,Y) to the next step is legitimate if the whole local
logarithmic derivative is analytic there. Endpoint vanishing of A_m for m>=2
alone is insufficient: delta*log²(delta) still has a logarithmic derivative.
The code's maximum-degree warning checks the coefficient arrays on the whole
element, which is the appropriate numerical safeguard.

### Independent scalar-mode reproduction

Use c=12 and the exact Y eigenfunction h(Y)=sin(Y/2), lambda=-1/2.
We normalize the modal functional as (2/pi) integral h(Y)f(Y)dY.
The corresponding homogeneous X solution is h(X)=sin(X/2). For a continued
modal solution f and r its projected source,

    [cos X*(h*f'-h'*f)]' = h*r.

Analyticity at Xa=pi/2 and evaluation at Xb=3pi/2 give the **positive-log**
coefficient B=integral_(Xa)^Xb h(X)r(X)dX/h(Xb).
The explicit inputs, without a psi20 evaluator or spectral continuation, are

    psi10=sin(beta/2)/sqrt(2)-Z*(cos(alpha/2)+sin(alpha/2)),
    psi21=-Z*(pi-2)*cos(beta)/(3*pi).

`x26_mode.py` integrates the Z-linear source using independent 45-digit
quadratures in X and Y. It obtains

    B = -0.771509959510067039452295971738812887078082...

The saved extend coefficient has the opposite sign, as predicted by the
negative-log convention, and agrees in magnitude through the displayed digits.
The cos(Y/2) mode is zero at this precision. This is a genuine independent
number check, although the quadrature is not an interval enclosure.

### What R1–R3 mean after the audit

* **R1:** no side-minus saved solve was supplied among the inspected inputs.
  The corrected side-minus spot check at Z=2, n64/N48, 256 bits finds
  maximum far-log coefficients 2.603e-48 (psi30), 1.582e-48 (psi40),
  and 1.969e-49 (psi41), with psi20 zero after thresholding. This supports
  numerical absence at that resolution; it is not an exact zero or an
  independent all-charge-sector proof.
* **R2:** co-sign the observed degree-one phenomenon as numerical evidence,
  with the odd-sign correction. Use extend or the corrected chain for first
  variations. Do not use the original chain's Z² psi41 or psi40 first variation
  as physical data.
* **R3:** co-sign the rational fit and termination numerically after the
  bookkeeping correction. Exact zero is stronger than the thresholded
  arrays establish. `chain.jl` drops trailing arrays below 1e-50; its stored
  zero triple variations are therefore censored at that scale.

## B. T5 blindness: the necessary three-way distinction

At a generic type-II point T5 is a polarized orthogonal rank-five system.
N²=0, rank N=2, and W=Im N is a maximal isotropic plane. Ker N=W-perp.
Consequently

    N_a N_b=0  <=>  W_b subset W_a-perp  <=>  W_a=W_b.

The second equivalence uses the fact that the nondegenerate one-dimensional
space W_a-perp/W_a has no nonzero isotropic vector.
But unequal planes have TWO possibilities:

| Intersection | Product | Alternating chain |
|---|---|---|
| dim 2, equal planes | zero | stops immediately |
| dim 1 | rank 1 | nonzero double, zero triple |
| dim 0 | rank 2 | generically nonterminating; paired-plane map is invertible |

The exact counterexample in `x26_exact.py` has basis e1,e2,f1,f2,h with
<ei,fj>=delta_ij, <h,h>=1. Put

    N_a=e1 wedge e2,    N_b=e1 wedge f2,
    (x wedge y)(z)=x*<y,z>-y*<x,z>.

Both operators lie in so(2,3), have rank 2 and square zero. Their image planes
are distinct and meet in <e1>. Direct exact multiplication gives
N_a N_b !=0, N_b N_a N_b=0. Thus the observed pattern is not, by itself,
evidence for T5 blindness or T5 absence.

### Which based pair was actually computed?

At the two strip edges, with fixed physical Y,

    (u_b,v_b)=(-v_a,-u_a).

This is not an automorphism of the labelled S_e2 problem that one may silently
use to identify its limiting planes. The strip transport also crosses the
alpha=pi and beta=0 lifts of u=-1 and v=1. The cosine maps are quadratic:
an alpha/beta semicircle gives a full meridian in u/v. The u=-1 half-exponent
reflection survives as an order-two action; it is not killed merely by using
cosine coordinates. The v=1 operator has repeated exponents 0,0,1,1,2;
one must include its unipotent action, not treat it as a harmless reflection.
A choice of detours and their order is essential on the bare T5 system,
even when the full physical hierarchy branch is analytic at those loci.

The saved `TA.pkl` contains the two circle meridians on **v=1/3**, based at
u=1/5. They are a different pair. Comparing their planes is useful as a
control and a proposed additional loop, but does not compute the strip
transport. Likewise a formal product of independently based reflection and
v=1 matrices is not a certified strip transport without the connecting paths
and braid identification. `x26_monodromy.py` explicitly labels such products
as experiments, not a solution of this identification.

**OPEN:** the actual W_a versus transported W_b for the specified strip.
A definitive certificate would give the five-state transport P along that
lifted path and the ranks of [N_a, P^-1 N_b P], N_a P^-1 N_b P and its triple.
This report does not replace this missing map with the rank test for another
pair. The source-dependent physical pairing remains a separate problem.

### Position of T5 and what a scalar chain can exclude

**Submodule of D.psi40:** a nonzero rational differential expression R psi40
realizing T5 would be an actual T5 scalar solution. Under the prior physical
analyticity hypotheses and the three-edge common-invariant computation, the
phase-11 submodule exclusion applies. The new scalar chain does not strengthen
that result without a certified extraction R and loop identification.

**Quotient:** a quotient of the differential module need not be realized as
one scalar R psi40. Solution local systems reverse the exact sequence:
a D-module quotient T5 gives a T5 subspace of the solution sheaf, and the
selected physical solution need not lie in it. Its tested word can terminate
without the word vanishing on that whole constituent.

**Middle:** one must first reach the relevant subobject and then pass to its
quotient, retaining the framing/pairing. A numerical word applied to the
physical scalar need not do either. In a monodromy-equivariant quotient of
the actual branch representation, a branch with nonzero projected component
must obey the projected word, but that component can lie in a proper word
kernel. Generic nontermination is a generic-position statement only, and
requires the disjoint-plane case above.

An operator identity for a word on the **entire** continuation space descends
to every subquotient and is genuinely position-independent. Testing one
branch and its one alternating chain does not establish such an identity.
Positive normal orders can also hide a constituent from a zeroth normal
limit, as in the earlier delta-times-period counterexample.

## C. Rational double variation and the next useful loop

Let delta=1-u²-v². The observed boundary value is the restriction of

    H=(2/9)*(u²+v²)-1/9-(14/45)*uv
     =(4u²-1)/27 - (7/90)*(4uv)
       +(4/27)*((u²+3v²-1)/2).

These are exactly the three implemented k=4 spherical harmonics.
**PROVED:** (Lambda²-32)H=0 and H|delta=0=1/9-(14/45)uv.
The double-varied lower coefficients and psi41,psi42 vanish in the numerical
chain, so the double-varied psi40 equation is homogeneous. If their vanishing
and the fitted boundary data hold exactly, analytic modal uniqueness at the
start edge makes the double variation equal to H throughout the strip,
not merely on the edge. H is analytic at the next edge, explaining termination.

This identifies explicit order-four polynomial terms that represent the
result. It does **not** identify the original order-four source atoms whose
separate variations give coefficients 1/27,-7/90,4/27. The initial added
Y4 harmonics have zero variation and cannot themselves be its origin.
The rational coefficients have been fitted and independently checked,
not derived from a term-by-term C3 source identity. Such an attribution
remains **OPEN**. The result is consistent with psi40 in the old category:
a polynomial harmonic is in that category. It does not prove membership,
and no elliptic constituent is detected by this scalar boundary value.

With t=tan(alpha/4), the fitted function is rational:

    1/9-(56/45)*t*(1-t²)*(1-6t²+t⁴)/(1+t²)^4.

**VERIFIED EXACTLY:** it satisfies the guessed third-order ODE. That ODE
annihilates the entire span of 1,sin(8 atan t),cos(8 atan t); it is not a
minimal rational-coefficient annihilator of this particular rational
function (which has a first-order rational annihilator). Its order three
must not be interpreted as a constituent rank.

### A concrete nonblind control and a separable follow-up

The saved opposite-circle meridians on v=1/3 provide a ready-made T5 control:
compute Im N_plus and Im N_minus in their existing common basis, their
concatenated rank, and N_minus N_plus N_minus. This pair is accessible from
the two physical X-edge arcs in the upper half circle, with distinct Y
positions. It avoids claiming that a generic or unspecified extra loop
will work. The numerical control is exported in `x26_TA.json`.

For the separable solver, the cheapest concrete next extension is to carry
the already nonzero far-edge variation back across the physical square and
test the **opposite physical X edge**, allowing the intervening logarithmic
edge as an additional element with a prescribed detour. A second option is
the X/Y-swapped continuation to a far Y edge. Before interpreting either
scalar zero, calibrate the actual new loop pair on T5 in one based frame.
The known opposite-circle pair supplies a nonblind operator control, but
there is no guarantee a selected physical branch has a nonzero projection
onto its detected T5 directions. An infinity loop is more expensive and
is not the first recommendation.

### Final numerical evidence

* Corrected n80 first variations versus the saved extend arrays: maximum
  absolute difference 1.562e-72 over the inspected coefficients/sectors.
* Corrected double-fit maximum over the n80 grid: 5.947e-60; other double
  sectors: <=3.016e-60; triple arrays zero after the 1e-50 cutoff.
* Original versus corrected psi40 Z² first log: maximum difference
  0.00932502912493509725429. Original psi41 Z² minus twice corrected:
  maximum 1.109e-48.
* Original versus corrected double variation: 1.239e-51.
* Corrected n64/n80 first-log interpolation comparison at the five specified
  test points: worst 3.186e-36. These are sampled convergence diagnostics.
* Saved line-A opposite-circle image concatenation has singular values
  1447.93, 99.12, 5.90, 0.1714, 9.66e-24: the two planes are disjoint
  numerically. The triple product has two nonzero singular values 5788.69
  and 29.95. Fourth/sixth alternating-word norms are 3603.70/14414.80.

An independent 40-digit, 65-term Taylor integrator reproduces circle rank 2,
u=-1 reflection rank 1 and involution residual 1.25e-20, and v=1 rank 2
with (T-I)^2 residual 1.66e-19. It also supplies two useful actual-T5 controls:

* Conjugate the Q2 circle inertia by the computed u=-1 reflection: image
  concatenation rank 3, product rank 1, triple norm about 2.4e-19. This
  realizes the shared-line terminating pattern in this very local system.
* Conjugate it by the explicitly based v=1 meridian: image concatenation
  rank 4 (fourth singular value 0.3466), product rank 2, triple singular
  values 23446.7 and 537.36. This is a concrete **VERIFIED numerical
  nonblind control** for the extra loop, before any scalar physical pairing.

The v=1 meridian is based by the path from (u,v)=(1/5,1/3), first upward
in the complex v plane by i/4, then horizontally to v=1+i/4, down to
1+0.004i, a positive small circle, and back. The C-line jets are glued to
the A-line jets by the exact exported 2-D jet matrix. The u=-1 meridian uses
the corresponding upper-u detour and radius 0.01. These fully specified
controls do not depend on assigning an uncomputed strip path to a word.
The additional `x26_strip.py` sign-symmetry/connector experiment also finds
rank-four image spans for several candidate words; this is **CONJ evidence
for nonblindness of the strip**, not a based homotopy certificate.

## Reproduction and limits

All numeric figures in the final metrics are non-rigorous floating-point
experiments. No energy calculation, repository edit, commit, recurrence,
agent, package installation, or original-script execution that writes into
Claude's directories is used. Importing `extend.jl` defines functions only;
its main-program serialization guard is false. Shared logs are append-only.
The literal actual-strip T5 decision and a source-term derivation of the
rational coefficients are unresolved; this is not a proof of Conjecture A.

### Resource exception

The first recorded Julia high-water RSS, on the side-minus n64 spot check,
was **2,518,417,408 bytes**, exceeding the requested 2 GB ceiling. Earlier
n64/n80 runs had no successful RSS measurement, so their compliance cannot
be asserted. A pending unrestricted repeat was interrupted as soon as this
was discovered. Subsequent solver runs use `--heap-size-hint=500M`; their
n80 measured peak is **1,698,938,880 bytes**, below 2 GB, and its numerical
results reproduce the earlier run. The final saved-data verifier is also
run with that heap target. This is a resource
constraint violation in this audit, not hidden by the later bounded rerun.

### Concurrent workspace changes

The final baseline comparison found concurrent changes to `codex_26.log`,
`claude11/chain.jl`, `claude11/guess_ode.jl` and new `claude11/lad*` /
`blind_u*` outputs. I did not write those paths. The immutable audited
chain text is preserved as `x26_chain_original.jl`, reconstructed exactly
from the copied solver before our two code edits and RSS instrumentation;
it agrees with the current chain text at final inspection, but differs
from the earlier baseline hash. The final audit lists these concurrent
changes rather than asserting the entire workspace remained unchanged.
The scientific computations and source derivations refer to the preserved
version. Later concurrent results have not been silently co-signed here.
