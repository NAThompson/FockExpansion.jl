# TASK_codex_19 — order-four membership and line-charge geometry

2026-10-23. Status words refer to individual statements. Computations are serial,
with single-threaded BLAS and no energy, recurrence or stencil jobs. The only
writes are `codex/phase9/x19_*` and authorized append-only shared logs.

## 1. ψ41^(2): independent connection check and purity

**PROVED / exact computer algebra.** Claude's `closed41.pkl` is a polynomial
combination of the five Codex atoms `(1,L,alpha,V,U)` from phase 2, with
`V_Claude = r12 V_Codex` and the same U. The coefficient of the sixth
log(s/R) atom is zero. Independently reconstructing the integer relations in
`coef5_ident.txt` reproduces every coefficient of the pickle exactly.

`x19_check41.py` applies the Codex column-oriented rational distance connection
and the full six-dimensional distance Laplacian, including both mixed distance
terms. All six coefficients of

    Delta_6 F + R² [ B(V1(5xi³/6-xi)+vs(xi-1/xi)) + 24 psi42^(2) ]

vanish identically. Here F=R⁴ψ41^(2), B=(pi-2)/(3pi), xi=r12/R,
V1=R/r1+R/r2 and vs=(r1+r2)/R. The two arbitrary a_h,b_h additions cancel
exactly, and every coefficient has homogeneous degree four in the Codex
normalization. This is an identity over the rational function field with pi,
not a residual at a few points.

**CO-SIGN** membership in the already established rank-five order-two module;
no new atom is required. This establishes an upper bound five for the single
function's generated module; no new claim that this particular function has
minimal cyclic rank exactly five is needed.

**VERIFIED independently: purity.** Let Y0=4cos²alpha-1 and
Y2=sin²alpha(3cos²theta-1)/2. With the ignorable azimuth removed,

    <Y0,Y0>=pi, <Y2,Y2>=pi/8, <Y0,Y2>=0.

The variables X=2(alpha-beta)/pi and Y=2(alpha+beta-pi)/pi,
beta=acos(sin alpha cos theta), map the physical domain to [-1,1]².
The measure becomes (pi²/8) sin alpha sin beta dX dY. Tensor Gauss quadrature
of our own actual order-two atoms, at 256 bits and n=16,24,32,48 per axis,
gives the pure inverse constants for the pickle's polynomial particular part:

    a_h = 0.00082106684957556469415831715509798526043967205959...
    b_h = 0.01992734793080898549051527717103411701281763036468...

n=32 and n=48 agree to more than 48 absolute decimal places. This is
quadrature convergence evidence, not an exact identification in Q(pi).

The constants extracted from `coef5.txt`, which fit `spec41.txt`, are instead

    a_h = 0.00388038725134403242522098865866299174984534...
    b_h = 0.01455939150640111472454125565875339387374006...

Those data contain the hierarchy-required harmonic component. Subtracting the
pure constants leaves +0.00305932040176846773106267... Y0
and -0.00536795642440787076597402... Y2. They must not be called a pure inverse.
With these spectral constants the independent formula matches 81 scattered
spectral points with maximum absolute discrepancy below 9e-43. With the pure
constants it has the zero projections above, and appropriately differs from
that full spectral dataset. This agrees with Claude's newest purity entry.
The independent Feynman `data2.txt` comparison is recorded in `x19_purity.out`;
only two harmonic constants are fitted, with separate held-out points.

**CITED / local inputs:** the source and the resonant hierarchy normalization
are the read-only repository's `src/fourth_order.jl` and `src/feynman_k4.jl`.
No Claude executable checker was imported or run.

## 2. ψ40 charge-sector membership tests

**VERIFIED numerical tests, not a nonmembership theorem.** `spec40.txt`
contains C0,...,C4 at E=a21=0. The fitted targets are exactly its C1,C2,C3.
All five order-two atoms, all eighteen classical generators and all 4*23 x16
states were evaluated independently in `x19_atoms.jl` from their defining
functions/integrals. The full family has 115 columns before polynomial
multipliers; duplicates are retained and handled by rank-revealing factorizations.
No symmetry restriction is imposed on the fitting span.

On R=1 use all monomials r1^i r2^j r12^k with j=0 or 1 and i+j+k<=d.
This is the complete degree-at-most-d polynomial restriction, since
r2²=1-r1², and contains (d+1)² monomials. These are slice membership tests;
by themselves they do not assert a homogeneous global distance representation.

The fixed geometric cut is 0.25<alpha,theta<pi-0.25, selected without looking
at target values. Of the 4,000 supplied points it retains 2,801: 2,240 training
and 561 held-out. Float64 panel rules with 256 and 512 nodes agree over every
retained point to max scaled discrepancy 1.034e-14. The earlier 128-node rule
failed the guard; its values were never fitted. Its worst state's convergence
sequence is saved in `x19_quadrature.out`.

The following table gives Float64 **held-out relative L2 residuals**, with an
SVD cutoff 1e-12 relative to the largest singular value. Numerical rank is a
conditioning diagnostic, not a rational-function rank certificate.

| Family | Degree | Columns | Numerical rank | C1 | C2 | C3 |
|---|---:|---:|---:|---:|---:|---:|
| order2 | 0 | 5 | 5 | 2.115e-01 | 1.725e-01 | 9.417e-02 |
| order2 | 2 | 45 | 45 | 2.161e-05 | 6.800e-05 | 1.942e-05 |
| order2 | 4 | 125 | 81 | 2.735e-11 | 3.741e-11 | 1.685e-11 |
| order2 | 6 | 245 | 100 | 8.012e-12 | 6.556e-12 | 5.034e-12 |
| order2_classical | 0 | 23 | 19 | 3.815e-03 | 1.457e-02 | 8.194e-03 |
| order2_classical | 1 | 92 | 67 | 1.308e-07 | 1.137e-06 | 3.692e-07 |
| order2_classical | 2 | 207 | 89 | 2.520e-09 | 6.083e-09 | 1.762e-09 |
| order2_classical | 3 | 368 | 103 | 1.077e-10 | 1.947e-10 | 3.938e-11 |
| order2_classical | 4 | 575 | 117 | 4.179e-12 | 4.800e-12 | 1.930e-12 |
| all115 | 0 | 115 | 103 | 1.037e-05 | 8.099e-05 | 3.075e-05 |
| all115 | 1 | 460 | 195 | 1.054e-08 | 4.070e-08 | 1.349e-08 |
| all115 | 2 | 1035 | 214 | 8.679e-10 | 1.368e-09 | 3.480e-10 |
| all115 | 3 | 1840 | 230 | 2.458e-11 | 4.110e-11 | 8.741e-12 |

The SVD cutoff sensitivity at d=2 is material: 1e-10 / 1e-12 / 1e-14
produce different effective ranks and residuals. Even the known degree-four
psi41 control has held-out residual 1.833e-11 in Float64 at the main cutoff.
Thus residuals near 1e-11 in that table cannot distinguish an identity.
Every training residual, held-out relative maximum, cutoff and rank is retained
in `x19_fit40_results.txt`.

**256-bit follow-up.** A deterministic 1,000-point subset has 800 training and
200 separate held-out points. Each panel uses a Gauss rule on panels graded to
the arctan poles and the complex quartic-root distance. 96 vs 128 nodes per
panel agree to 1.906e-58 at five validation points. A second implementation,
adaptive mpmath quadrature of the prior independent x16 actual-integral
checker, compares all 92 states at two fresh cache points; see
`x19_panelcheck.json` for its errors. This guards against a shared finite-rule
convergence error.

Column-scaled pivoted BigFloat QR uses a relative R-diagonal threshold 1e-42.
The table shows **training / held-out relative L2 residuals**. The full degree-1
family is still ill-conditioned and extrapolates substantially; its small
training residual must not be treated as an identity.

| Family | Degree | Columns | QR numerical rank | C1 train / held | C2 train / held | C3 train / held |
|---|---:|---:|---:|---:|---:|---:|
| all115_big | 0 | 115 | 107 | 7.224e-06 / 1.642e-05 | 4.948e-05 / 3.747e-04 | 1.863e-05 / 1.372e-04 |
| all115_big | 1 | 460 | 418 | 2.448e-15 / 1.955e-10 | 4.955e-15 / 2.395e-10 | 1.264e-15 / 5.497e-11 |
| order2_big | 2 | 45 | 45 | 2.273e-05 / 2.429e-05 | 7.109e-05 / 7.438e-05 | 2.030e-05 / 2.101e-05 |
| order2_big | 4 | 125 | 125 | 1.791e-14 / 3.923e-14 | 1.318e-14 / 1.784e-14 | 5.767e-15 / 1.241e-14 |
| order2_big | 6 | 245 | 245 | 1.546e-20 / 5.059e-20 | 1.036e-20 / 3.230e-20 | 3.339e-21 / 1.150e-20 |

The degree-four psi41 control has full QR rank 125 and residual
1.720e-75 training / 1.974e-75 held-out. This confirms that the high-precision
pipeline can recover a known identity at the working precision; the psi40
residuals in the tested spans are much larger. The initial unrestricted
high-precision sweep stopped at degree one (460 columns). A symmetry reduction
then enabled the full singlet degree-two search below. Higher degrees in the
full family have only the Float64 evidence, where rank truncation is severe.

**Full degree-two singlet follow-up at 256 bits.** The atom family is closed
under electron exchange with constant coefficients. Center alpha by pi/2; the
order-two atoms have exchange parities (+,+,-,+,-). Classical L/T pairs
exchange with a minus sign, the two remaining L/log pairs with a plus sign;
classical a is odd and p,a²,p²,1 are even. The panel families exchange
1<->2 and 3<->4. Combining paired sums/differences gives 60 even and 55 odd
generators, including presentation dependencies. Direct actual-function checks
at two points verify these parities to 9.0e-77.

Let s=r1+r2 and d=r1-r2. Since s²+d²=2 on R=1, every degree-at-most-two
polynomial reduces to s^i d^j r12^k with j=0 or 1. The complete singlet span is
thus even atoms times six even monomials plus odd atoms times three odd
monomials: 60*6+55*3=525 columns. This is the singlet image of the full
1,035-column degree-two family, not a chosen subset of its atoms. The same
800/200 point split and QR cutoff give numerical rank 463.

| Multiplier | C1 train / held | C2 train / held | C3 train / held |
|---|---:|---:|---:|
| 1 | 8.017e-27 / 1.004e-17 | 1.333e-26 / 1.977e-17 | 2.520e-27 / 3.318e-18 |
| delta | 2.226e-20 / 1.442e-11 | 5.064e-20 / 3.890e-11 | 1.407e-20 / 1.132e-11 |
| delta2 | 8.466e-17 / 1.960e-08 | 2.895e-16 / 2.484e-08 | 9.859e-17 / 1.412e-08 |
| Dclear | 5.104e-13 / 1.693e-04 | 5.636e-13 / 3.714e-04 | 2.656e-13 / 1.978e-04 |

A separately manufactured member of this full panel-polynomial span has
training residual 1.328e-47 and held-out residual 2.425e-38. The fit remains
ill-conditioned and amplifies errors on held-out points, but the Cj residuals
are far above this control. The complete output is in
`x19_fit40_singlet_results.txt`; peak RSS was 806,502,400 bytes. Full-family
high-precision singlet degree>=3, unrestricted exact rational-rank certification
and unbounded-degree membership remain **OPEN**.


**Denominator controls.** Put

    delta=4 r1² r2²-(1-r12²)²,
    Dclear=r1³ r2³ r12³ (2-r12²)² delta²  (R=1).

The homogeneous version includes an extra factor R². The exact geometry checker
verifies that it is divisible by every one of the 92 exported first-order x16
row clearings. We fit delta*Cj, delta²*Cj and Dclear*Cj as well as Cj.
The following are held-out relative L2 residuals in the same 115-state span:

| Multiplier | Precision / degree | C1 | C2 | C3 |
|---|---|---:|---:|---:|
| all115_delta | Float64 / 3 | 7.845e-08 | 3.451e-07 | 1.343e-07 |
| all115_delta2 | Float64 / 3 | 3.838e-06 | 9.291e-06 | 3.928e-06 |
| all115_commonclear | Float64 / 3 | 2.370e-04 | 2.319e-04 | 1.019e-04 |
| delta | 256 bit / 1 | 4.783e-06 | 2.059e-05 | 6.791e-06 |
| delta2 | 256 bit / 1 | 6.643e-04 | 2.724e-03 | 1.153e-03 |
| Dclear | 256 bit / 1 | 9.830e-02 | 4.372e-02 | 2.062e-02 |

These are residuals for the multiplied targets, hence weighted residuals for
the original functions. Multiplication raises the required numerator degree;
failure at these low degrees does not exclude a rational representation with
these denominators. The existing derivative clearings need not be the only
possible denominators of a membership relation. The tests specifically try the
requested discriminant/row-clear options and **do not** conflate their failure
with a proof that every denominator choice fails.

**Assessment:** no exact membership of C1,C2,C3 in the tested old spans was found.
Their approximation improves quickly with degree. Neither an exact old-module
representation nor its impossibility has been established. The high-precision
control shows that these residuals are not merely the arithmetic floor. Exact
nonmembership even at a tested finite degree is not certified: discarded tiny
directions could matter with large coefficients. Higher degrees and other
denominator choices remain possible. No coefficient
PSLQ reconstruction or claimed closed form is based on these approximate fits.


## 3. Line charges, a concrete elimination, and its strict interpretation

Write the sphere source as y, target as x=(u,v,z0,0), with
z0²=1-u²-v². For panel 1 set

    D=1+q⁴, P=(q⁴-1,2q²), p=P/D.

**PROVED / VERIFIED algebraically:** P·P=D² and W²=D-P·y. Panel 3 uses
P=(-2q²,1-q⁴); panels 2/4 are reflections. Thus the radicand observation
is exact, not a numerical resemblance.

For a Coulomb point a in {e2,e1,-e1}, the beta identity is

    1/[sqrt(1-a·y)sqrt(1-p·y)]
      =1/pi int_0^1 dt/[sqrt(t(1-t))(1-r·y)],
    r=t a+(1-t)p.

The constants translating these factors to v0=1/xi and
v1=csc(alpha/2)+sec(alpha/2) are elementary. For a polynomial numerator of
degree at most two, applying T4 under this integral uses exactly the repository's
elementary zonal inverses of `1/ell`, `log ell`, `ell log ell`:

    T4[y_i/ell] = -d_(r_i) T4[log ell],
    T4[y_i y_j/ell] = d_(r_i)d_(r_j) T4[ell log ell].

The inverse here is the pure inverse, including its harmonic subtraction.
All manipulations are first on regular interior branches; the physical
integrability and projection qualifications from x18 still apply.

### The simultaneous radicals, with the density normalization retained

Put a_x=a·x, p_x=P·x, and define

    A=1-a_x²,
    B=a·P-a_x p_x,
    C=D²-p_x²=(D-P·x)(D+P·x),
    Q(q,t)=|tDa+(1-t)P|²-[(tDa+(1-t)P)·x]².

Then Q=D²(r·r-(r·x)²). Its exact t-discriminant is

    disc_t Q = -4 D² (a1 P2-a2 P1)² z0²,
    B²-AC = -(a1 P2-a2 P1)² z0².

In particular, Q alone is a conic in t with a split discriminant after adjoining
i z0. It is incorrect to infer a K3 just by counting two integration variables.
The additional beta measure and the external factor 1/sqrt(D) matter.
The simultaneous angular radical and this correctly normalized measure exhibit

    Y0²=D t(1-t) Q(q,t).

Use s=t/(1-t), chi=D s, X=A chi,
Y=A(1+s)²Y0. Direct substitution gives the Weierstrass cubic

    Y² = X³ + 2B(q) X² + A C(q) X,
    discriminant = -64 A² z0² C(q)² (a1 P2-a2 P1)²,
    j = 64(4B²-3AC)³/[A² C²(B²-AC)].

All these identities, and the four W² maps, are checked exactly by
`x19_geometry.py`. Retaining D avoids a misleading quadratic twist.

**PROVED (generic auxiliary geometry).** In q, the Weierstrass coefficients
have degrees at most 4 and 8. At generic target x:

* C has eight simple roots, disjoint from a1 P2-a2 P1; each gives a fibre I2.
* For a=±e1, a1 P2-a2 P1=±2q². There are I4 fibres at q=0 and infinity.
  Thus the singular-fibre Euler numbers total 8*2+2*4=24.
* For a=e2, the cross product is 1-q⁴. Its four simple roots give four
  additional I2 fibres, hence 12*2=24. Infinity is regular generically.

At these fibres c4 is nonzero, so the equations are minimal and the stated
multiplicative fibre types follow from the discriminant orders. The associated
Weierstrass line bundle on P1_q is O(2). Resolving the fibre rational double
points gives trivial canonical bundle by adjunction; H1(O)=0 because the base
is P1 and the line bundle has positive degree. This is a K3 surface.
The exact genericity conditions include z0!=0, A!=0, and absence of the listed
root collisions; the checker verifies simple/disjoint fibres at a fresh rational
target, in addition to the symbolic factorizations.

It descends under q -> -q to an elliptic surface over z=q² with coefficient
degrees 2 and 4, six I2 fibres and line bundle O(1): a rational elliptic surface.
The natural q parametrization is the quadratic base change. For a=±e1 the
base change is ramified at two I2 fibres, producing the two I4 fibres. For e2
the two branch fibres are regular. This explains the K3 here concretely.

**VERIFIED:** the j-invariant varies with q at fixed rational target (two exact
values for each Coulomb point are in `x19_geometry.json`). These exhibited
elliptic fibres therefore are not the fixed old target curve

    E_old(x): w²=(1+u)+(1-u)q⁴-2v q²,

whose modulus is fixed once x is fixed. This is a statement about this
fibration, not an obstruction to all other reductions or alternative fibrations.

### Why the complete panel is not yet a (q,t) integral of these U-functions

The nonconstant panel density is essential. In physical N0/N2 it contains
beta_q(b)²=atan(q/h)², h²=(1+b)/(1-b), with b a sphere coordinate.
It is not a coefficient depending only on q. Thus T4[v_a panel] cannot be
obtained by merely multiplying the elementary two-Coulomb U-function by that
density evaluated at the target. In particular T4[D(y)f(y)] is not
D(x)T4[f](x) for this density.

**PROVED explicit correction.** For fixed b on the physical interior,

    beta_q(b)² = (1-b²) int_0^q dr int_0^q ds /
      [(1+b+(1-b)r²)(1+b+(1-b)s²)].

The integrand identity follows by squaring
atan(q/h)=int_0^q h dr/(h²+r²), and is checked rationally.
After this substitution the density contributes two additional linear sphere
denominators. Combining them with the two Coulomb denominators by a Feynman
simplex gives a valid higher-parameter integral of elementary zonal inverses
(and their parameter derivatives). For example

    1/(ell Lr Ls)=2 int_(lambda_i>=0,sum lambda_i=1)
                    d lambda1 d lambda2 /
                   [lambda1 ell+lambda2 Lr+lambda3 Ls]³.

The cubic pole is an intercept derivative of the inverse-linear source.
Polynomial numerators are handled by parameter derivatives. This is an
explicit sufficient representation, not a minimal parameter count. The
fixed-endpoint involution used by x16 removes moving upper endpoints; without
that step, the original min(1,h) also introduces a source-dependent indicator.
Neither issue can be omitted when interchanging integrals.

**CO-SIGN** the line-charge identity and the elementary (q,t) reduction for the
source-independent angular-density pieces. **DISSENT** from extending that
reduction unchanged to the full weight-two panel density. The extra density
parameters can add further discriminants; the displayed K3 is not a complete
classification of that larger integral.

**CONJ:** the earlier old-curves/one-higher-weight description remains an
unproved hypothesis; these calculations do not promote it to a theorem.

**OPEN:** whether a nonzero K3-period component survives the complete physical
ψ40 combination, whether cancellations/identities reduce it to iterated integrals
on only old E_(H,k) curves, and the minimal post-integration differential rank.
The exhibited K3 and varying j do not prove necessity of a new essential curve
or a K3 period. Conversely, a line-charge interpretation or unchanged real
Chebyshev convergence does not prove old-curve closure. No order-four recurrence
is warranted from either inference.

## 4. Verification, resources and handoff

The independent Feynman comparison uses 80 training points solely to fix its
two harmonic gauge constants; all remaining data2 points have maximum absolute
residual 3.068e-47. The separate adaptive-mpmath panel check has 184 comparisons
and max scaled error 5.441e-63. Exact geometry checks all pass, including the
normalization twist, three cubic/discriminant eliminations, varying j and the
row-clear divisibility control. These checks do not certify unmeasured
post-integration rank or an essential K3 component.

| Serial calculation | Peak RSS, bytes |
|---|---:|
| Exact Codex-connection PDE / coefficient check | 172,310,528 |
| Purity, spectral and independent Feynman check | 323,878,912 |
| Float64 membership sweep, max matrix 2,801 x 1,840 | 858,161,152 |
| 256-bit membership sweep, max matrix 1,000 x 460 | 1,286,881,280 |
| Full 256-bit degree-two singlet test, 1,000 x 525 | 806,502,400 |
| Exact geometry / divisor eliminations | 75,644,928 |
| Independent adaptive panel check | 96,354,304 |

Every job was serial. The maximum measured RSS is below 2 GB. No subagents,
energy calculations, full recurrence expansions, installations or outside-directory
writes were used. Python imports suppress bytecode; Julia disables startup,
history and compilation of new cache entries and names an in-directory first
depot. The repository and Claude artifacts are read-only inputs.

Reproduce from this directory with `venv/bin/python` and `PYTHONDONTWRITEBYTECODE=1`
for `x19_check41.py`, `x19_geometry.py`, `x19_panelcheck.py` (the last consumes
the BigFloat fit cache). Run Julia scripts `x19_purity.jl`, `x19_fit40.jl`, then
`x19_fit40_big.jl` and `x19_fit40_singlet.jl`, one at a time, with

    JULIA_DEPOT_PATH="$PWD/codex/phase9/x19_depot:/Users/nathompson7/.julia"
    JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
    julia --project=/Users/nathompson7/FockExpansion.jl --compiled-modules=existing \
      --startup-file=no --history-file=no --threads=1 SCRIPT

The report tables are transcribed automatically from the saved fit outputs;
`x19_tables.py` refreshes the fit section from those outputs. Scientific checks
are rerunnable. The final audit records input
hashes, verifies saved metrics, shared-log prefix preservation and observed
concurrent changes separately from our authorized writes.

**Final co-sign/dissent:** co-sign the psi41 rank-five membership, with the
pure/full harmonic distinction above; co-sign the four line-charge radicands;
co-sign the source-independent two-parameter zonal reduction. No exact psi40
old-module membership was found in the tested spans. Dissent from inferring
unrestricted nonmembership from finite-degree or truncated-rank fits, from
pulling source-dependent weight-two densities through T4, or from declaring
a K3 contribution essential solely because the auxiliary surface is K3.
The full function class and minimal rank of psi40 remain OPEN.
