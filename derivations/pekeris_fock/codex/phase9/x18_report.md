# TASK_codex_18 — ψ40 function class and induction audit

2026-10-23. All computations serial. Read-only repository inputs and their hashes
are in x18_verify.json. No recurrence, energy calculation, or cyclic-vector
reconstruction was attempted. Status words apply to individual claims, not this
whole research problem.

## 1. Explicit pure-inverse representation and charge split

**PROVED.** Put L=Λ²−32=−4(ΔS³+8), P=P_degree2, T=L⁺, with the inverse zero
on P. With ordinary sphere surface measure dσ of total volume 2π²,

    Tf(x)=−1/4 ∫S³ G(arccos(x·y))f(y)dσ(y),
    T²f(x)=1/16 ∫S³ G2(arccos(x·y))f(y)dσ(y),
    G(γ)=[(γ−π)cos3γ+sin3γ/6]/(4π²sinγ),
    G2(γ)=[(γ²−2πγ)/(48π²)+1/72+1/(864π²)]sin3γ/sinγ.

Here LT=1−P; no assertion LT=1 on arbitrary sources. The verification checks
the off-diagonal equation, delta flux, purity through the resonant moment,
and nonresonant spectral moments. Thus the signs and factor 12/16 below
are independently checked, rather than inferred from numerical repository values.

Write V=v0−Zv1, v0=1/ξ, v1=csc(α/2)+sec(α/2), and
ψ30=Σ_(j=0)^3 Z^j A_j, ψ20=Σ_(j=0)^2 Z^j B_j.
A_j include their E and a21 dependence; missing indices are zero. Write
ψ41^pure=Z q1+Z² T f2+Z³ q3, where q1,q3 and f2 are the explicit
functions in src/fourth_order.jl. Then the charge coefficient of pure ψ40 is

    C_j = T[12 δ(j,1)q1+12 δ(j,3)q3
              −2 v0 A_j+2 v1 A_(j−1)+2E B_j]
            +12 δ(j,2) T² f2,                 j=0,…,4.

This is an exact, nonrecursive-in-ψ41 representation. The arbitrary singlet
harmonics a40 Y40+a42 Y42 must be added afterwards. The entire harmonic part
of ψ41 and ψ42 is killed by T, but is still needed for solvability of the
original unprojected equation. Do not confuse projecting an incompatible
source with proving compatibility.

The explicit hard source is, in particular,

    C2 = T[−2v0 A2+2v1 A1+2E B2]+12T²f2,
    C3 = T[12q3−2v0 A3+2v1 A2].

Let A_j=P_j+D_j+Σ_m E_(m,j), using exactly the classical and four-panel
charge decomposition in codex/phase3/structure.md. This specifies every
integrand, including coefficients and fixed real branches. The hard panel /
Clausen contributions are T[−2v0(D_j+Σ E_(m,j))+
2v1(D_(j−1)+Σ E_(m,j−1))]. There is also a j=1 contribution from A1;
restricting attention to Z²/Z³ must not silently discard it. This formula
retains the full charge split rather than assuming all transcendental terms
start at Z².

A concrete angular quadrature for any one of these terms is

    Tf(α,θ)=−1/4 ∫_0^π da sin²a ∫_0^π dt sin t
                    K(α,θ;a,t) f(a,t),
    K=∫_0^(2π) G(arccos(A+D cosφ))dφ,
    A=cosα cos a+sinα sin a cosθ cos t,
    D=sinα sin a sinθ sin t.

Replace −K/4 by K2/16 for T². Substituting the explicit fixed-endpoint panel
integrals gives a fully specified iterated definite integral for every hard
piece; this is a Green representation, not a claim of a two-dimensional
*elementary-integrand* representation. The kernel azimuthal integral is a
third integration before it is evaluated. Source-panel integration adds another.

**PROVED (resolvent algebra).** If (Λ²−12)χ=hχ then

    Tχ=(Thχ−χ+Pχ)/20.

Likewise ∂a21 ψ40^pure=−ψ41^pure/(ZB), B=(π−2)/(3π), under the
repo's a21 normalization; the expression has a finite Z→0 continuation.
This follows from differentiating the hierarchy and Tψ42=0. The E derivative
reduces through the displayed resolvent identity and algebraic sources.
These reductions leave the pure-charge V times classical/panel core.
**CITED/attributed:** the repository README supplies the Z⁰/Z⁴ closed forms
and numerical checks; they were read, not independently rerun here.

For algebraic two-centre sources the identity

    1/(d± ξ)=1/(√2 π) ∫_0^1 dt/[√(t(1−t))(1−p±(t)·x)],
    p±=(±t,1−t,0,0)

reduces T to the elementary zonal ODE W''+9W=−f sinγ/4+κ sin3γ.
It also handles Vχ through its defining inverse/iterated integration, but
same source letters alone do not prove closure in the five order-two atoms.

## 2. Curves, candidate singularities, and the strict classification limit

**PROVED.** The old panel curve in the rational presentation is

    E_(H,k): w²=s⁴+k s²+H,
    disc=16H(k²−4H)²,
    j=16(k²+12H)³/[H(k²−4H)²].

It has genus one away from the discriminant. Endpoint/pole incidences add
T(1)=1+k+H=0, H=0, and resultants with s²+1 and s²+H:
T mod(s²+1)=1−k+H; T mod(s²+H)=H(H−k+1).
These are singularities of a presentation; some cancel for physical moments.
The old physical discriminant is the collinear triangle locus already
certified in phase 8. Infinity belongs in any compactified singularity audit.

**PROVED (kernel curve).** Applying the Stieltjes identity

    (π−γ)/sinγ=∫_0^∞ dr/(r²−2r cosγ+1)

and azimuthal integration gives the curve

    E_K: w²=(s²−4)[(s−2A)²−4D²],   s=r+1/r,
    roots {−2,2,2(A−D),2(A+D)},
    disc=65536 D²[((A−D)²−1)((A+D)²−1)]².

To see the curve without interchanging divergent terms, divide
T3(C)/(s−2C), C=A+D cosφ. Its remainder is T3(s/2);
∫dφ/(s−2C)=2π/√((s−2A)²−4D²), while dr/r transforms into
2 ds/√(s²−4) on the two halves of r>0. Keep the polynomial
subtractions together: their separate integrals at infinity diverge.
This yields the regularized incomplete elliptic kernel recorded in the README.

Set c−=A−D,c+=A+D and
λ=((1+c−)(1−c+))/((1+c+)(1−c−)). Then

    j_K=256(1−λ+λ²)³/[λ²(1−λ)²].

Cross-ratio permutations give the same j. The generic genus is one;
D=0, c±=±1 are degenerations, not evidence of genus greater than one.
The kernel family depends on both target and source, unlike a fixed
source-panel modulus. For fixed generic source, varying the target varies j_K.
Consequently there is no generic identification with that fixed source curve.
This does NOT exclude an integrated identity that removes this auxiliary curve.

**PROVED candidate integration loci, not a complete Landau variety.** In the
angular integral the kernel diagonal/antipodal loci are x·y=±1;
the Coulomb source points are y=±e1,e2; the old panel discriminants,
angle/log letters and integration endpoints remain. Collisions of the moving
kernel singularities with source points project to target x=±e1,±e2
(including antipodal candidates that may be removable). Complexification also
requires pinches satisfying denominators=0 and stationarity in the integrated
variables, plus boundary versions; simply taking their union is insufficient.
On a rational sphere chart these equations are algebraic after adjoining the
explicit square roots. Eliminating them for the hard panel convolution is OPEN.
No complete list of new target branch divisors has been obtained.

**OPEN:** ψ40 reduction to iterated integrals on only the old elliptic curves;
a weight grading increased by exactly one; necessity of new elliptic curves,
higher-genus periods or K3 geometry. Integrating an elliptic family against
another parameter-dependent kernel defines a higher-dimensional period problem.
Neither counting integration variables nor multiplying two square roots proves
a K3 or a higher-genus *irreducible* contribution. The kernel's genus-one curve
is established; any claimed additional essential curve needs an elimination /
Picard–Fuchs / monodromy calculation. No such claim is made here.

## 3. Structured holonomy; rank feasibility

**CITED:** Oaku, “Algorithms for integrals of holonomic functions over domains
defined by polynomial inequalities”, §§1–3, particularly the boundary-indicator construction and
integration discussion, https://www.lab.twcu.ac.jp/oaku/Oaku2013.pdf.
It treats holonomic distributions and boundary indicator factors, rather than
assuming that integrating annihilators discards boundary terms. His elementary
example ∫_0^1 exp(xy)dy explicitly demonstrates the latter error.

**PROVED conditional function-class statement.** Each fixed-order Green
representation is a compact direct image of a holonomic integrand provided its
physical branches, locally integrable distributional extensions, products and
boundary indicators are constructed in that category. Lower-order finite
connections, algebraic Coulomb factors and the inverse-trigonometric kernels
supply the local holonomic factors. This places ψ40 in a holonomic direct-image
class under those conditions; it does not classify its motives or weight.

A defensible structured construction retains separate blocks for each classical
atom and panel module, tensors each with the scalar Green kernel module, and
then performs integration-ideal reduction including endpoint distributions.
Only the resulting target connection permits derivative-row ranks at finite-
field points. The rank of the *integrand* connection is not the rank of its
integral. Values of arccos, π, periods or real boundary constants have no
canonical finite-field reduction, although rational annihilating matrices do.

**OPEN / not feasible in this small-job audit:** no post-integration rational
connection is available for the hardest core, so no honest holonomic-rank
measurement is reported. The two exact finite-field kernel fibres in the
verifier test smoothness only. Sampling the existing 23-state source matrices
would measure the wrong object. No cyclic vector or large elimination was run.

## 4. H3: what survives the audit

**DISSENT from H3 as written being a complete proof.**

1. Coordinate field: use algebraic sphere charts (or Cartesian sphere equations),
   not trigonometric angle coefficients over Q(α,θ). arccos(x·y), sinγ and integer
   multiples of γ are holonomic after algebraic substitution on branch charts.
2. The quoted integer-frequency kernel is a resonant spherical inverse formula.
   For physical index k, Λ² eigenvalues are 4n(n+2), so even k resonates at
   n=k/2; odd k has no integer resonance. The Green frequency is (k+2)/2,
   not k+1 with the conventions in this repository. At k=4 it must be 3, not 5.
   Nonresonant kernels use the appropriate regular-antipodal solution; one
   cannot reuse a resonant subtraction formula unchanged at every k.
3. The inverse is the pure inverse: finite harmonic projectors, nonzero spectral
   denominators, and compatibility of each source must be accounted for.
   Projected inverses always exist but solve the unprojected recurrence only
   when compatibility holds. Free harmonic data at each fixed order are finite
   constants times polynomial harmonics and do not obstruct holonomy.
4. Compact sphere integration is proper. A rational chart makes it unbounded:
   use finite bounded semialgebraic charts and their boundary indicators or
   compactify with the requisite distributional extensions. Branch/endpoint
   terms must be included, not set to zero by a telescoping convention.
5. Products of arbitrary singular distributions are not automatically defined.
   Coulomb factors have 1/distance singularities at finitely many S³ points.
   Bounded lower coefficients make each source L² there, since r² dr/r² is
   integrable. Elliptic inversion maps L² to H²⊂C⁰ on the compact 3-sphere;
   hence boundedness propagates for the finite descending log-sector induction,
   assuming compatibility. This supplies useful convergence control for the
   physical sources, but not for arbitrary raw derivative-panel primitives.
   The kernel ~1/distance is locally integrable; even at coincident source and
   target Coulomb points its product behaves at worst like 1/r², integrable
   in dimension three. For D-module closure one must additionally justify that
   the chosen extension of this integrable branch is holonomic; use the cited
   distribution framework, not an unrestricted multiplication assertion.
6. Descend p at fixed k, starting p=floor(k/2); then increase k. There is a finite
   number of operations at every k. This is not induction from a “top” of the
   entire infinite tower. Fixed log^p R is holonomic for fixed p, and algebraic
   radial substitutions preserve fixed finite-order holonomy. An infinite sum
   over k,p has no holonomic-closure guarantee or global convergence theorem.

**CO-SIGN the corrected conditional H3:** once these distributional and spectral
hypotheses are supplied, compact direct-image closure gives finite holonomic
rank at each fixed order. This is a substantially stronger justified prospect
than “the PDE has polynomial coefficients” (which alone proves nothing about
holonomic rank in several variables), but it is not a controlled-size module.
Finite products have tensor rank bounds from known input ranks; integration
requires cohomology, singular divisor and degree data. H3 supplies none of
these quantitative data. There is no rank-growth bound in k derived here,
no uniform finite atom family, and no all-orders recurrence theorem.

**H1:** same letters are useful evidence, not membership in the old module.
Co-sign testing it; no independent proof or disproof is claimed (Claude owns it).
**H2:** co-sign the Green convolution formulation; the old-curves/weight+1
classification remains CONJ. Dissent from presenting it as established.

## 5. Verification and handoff

x18_verify.py is a standalone serial SymPy/mpmath checker. It checks exact
kernel differential identities, flux, quartic discriminants and old-panel j,
ten high-precision spectral moments and two exact finite-field smoothness tests.
It records input hashes and process peak RSS. These do not validate all of
ψ40 numerically or certify the open classification/rank questions.

Task coverage: charge split and explicit hard Green integral DONE; kernel curve
and candidate loci DONE with incomplete projected singular locus explicitly
OPEN; precise old-versus-new essential function class OPEN; optional hardest-
piece rank NOT MEASURED for the stated structural reason; H3 audit DONE.
No order-four recurrence should be launched on the basis of this report.
