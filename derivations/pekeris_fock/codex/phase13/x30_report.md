# TASK 30 — ordered path audit and coefficient-class limitations

2026-10-24. **VERIFIED** the numerical distinction between C, the sum, and
individual A/B in the specified finite ansatzes. **PROVED** the ordered-word
sum identity and parameter-linearity of the moments. **OPEN** membership of
A over unrestricted rational coefficients; no exact nonmembership
obstruction is supplied or justified by these fits. Consequently this task
does not produce the requested global decision. **CO-SIGN** “no known
exact presentation of A”; **DISSENT** from treating this as a proof that
only the sum belongs to the order-three alphabet.

## Independent evaluator and conventions

Let W0=wv0 and W1=wv1 as in chain.jl (both positive Coulomb pieces), S_k
be the pure projected inverse taking the weighted RHS, and P2 the weighted
orthogonal harmonic projector at k=2. For an ordered triple (U,V,W), our
own unrolled recurrence is

    f10 = S1(-2U),
    f21 = P2((V/w) f10)/4,
    f20 = S2(8w f21 - 2V f10),
    f31 = S3(-2W f21),
    f30 = S3(10w f31 - 2W f20).

A=(W1,W0,W0), B=(W0,W1,W0), C=(W0,W0,W1). We share only
claude9/fock_spectral_copy.jl's low-level grid, inverse, harmonics, and
inner product. No fit_phi3 executable is imported. Values are interpolated
by a direct tensor Chebyshev transform, independently of Claude's
barycentric interpolator. BigFloat precision is 256 bits; n=64 and n=80
are recomputed on all 1,000 x19 points. Maximum absolute path discrepancy
is 1.228e-51. This is a numerical convergence check, not an enclosure.
A separate n=24 diagnostic checks polynomial interpolation and agrees
with full-Z sfock coefficient extraction to approximately 4e-74 at one
point. The latter is an implementation cross-check, not an independent
continuum proof.

The recurrence is trilinear, including its projected upper-log correction.
Expanding each slot W0-ZW1 proves, without a numerical fit,

    [Z] psi30(Z) = -(A+B+C)

for the same E=0 and zero-free-harmonic chain convention. This is a
polarization identity; it alone is not a proof of equality to an independently
chosen closed-form solution without matching the normalization.

## Important correction to the supplied exported “cl1”

fit_phi3d.jl creates column 1 of z1parts.txt as

    (classical(alpha,p0,1)-classical(alpha,p0,-1))/2.

But the repository classical formula has degree three in Z. Thus this
column is **cl[Z]+cl[Z^3]**, not cl[Z]. Its cubic coefficient is exactly

    e3 = -(r1+r2)*(5*(r1+r2)^2-3)/36.

The order-three panel pieces in columns 2:5 are the genuine four Z terms.
Therefore the correct fixed-coefficient control, without any coefficient
fit, is

    A+B+C = e3 - z1parts[:,1] - sum(z1parts[:,2:5]).

The fixed-coefficient check has max absolute residual 9.388e-64 and
held-out relative L2 residual 2.766e-63. No fitting step selects these
coefficients. All four moment coefficients are -1 in this presentation. The coefficient
on the saved first column is also -1, with the explicit classical correction
+e3. The uncorrected residual has maximum 0.2749854495; it is exactly the
missing cubic term to numerical accuracy. This does not invalidate claim
(ii)'s *span* inclusion: e3 is already in the degree-three classical span.
It does invalidate reading the odd extraction as a pure Z coefficient.
Also, the CONTROL print in fit_phi3e/f adds the unrelated general-moment
columns as if they were parts of the physical decomposition; that printed
quantity is not the decomposition residual. Our direct control uses only
the four actual Z panel parts.

## Finite fits (VERIFIED numerical evidence)

The original x19 point order is retained; rows 1:800 train and 801:1000 are
held out. Classical atoms are the first 23 x19 atoms, multiplied by
r1^i r2^j r12^k with j=0,1 and i+j+k<=d. The j reduction uses
r1^2+r2^2=1. Pivoted QR uses training-column normalization. Thresholds
1e-40, 1e-50, 1e-60 give identical retained rank for the degree-three tests.
All residuals below are relative held-out L2; outputs retain training
residuals too.

| ansatz | columns | rank | A | B | C | A+B+C |
|---|---:|---:|---:|---:|---:|---:|
| classical, d=3 | 368 | 304 | 2.950e-21 | 1.635e-21 | 7.816e-76 | 2.936e-21 |
| classical d=3 + saved cl1 + 4 Z parts + 24 general parts | 397 | 329 | 3.294e-24 | 2.097e-24 | 3.109e-75 | 8.232e-63 |
| classical d=4 + same 29 parts, threshold 1e-60 | 604 | 500 | 2.487e-29 | 1.583e-29 | 2.615e-74 | 1.391e-61 |

The general parts are both panels 1/2 and panels 3/4 as supplied. These
results reproduce the claimed scale separation, though not the literal
1e-22 residual cited for the larger ansatz. Different retained spans and
conditioning can change those digits. The C inclusion is exceptionally
well resolved. The sum's residual is compatible with precision of the
closed-form parts and the classical redundancy. A/B have no discovered
identity here; **“does NOT fit exactly” is not proved**. Neither QR threshold
stability nor a small finite residual supplies a certified algebraic rank
or an exact nonmembership witness.

## Higher parameter degrees and the rational coefficient question

**PROVED:** in the repository definition, for every fixed panel geometry,

    M(c0,c2) = c0*M(1,0) + c2*M(0,1).

The integration limits, h, b, y, and the denominator are independent of
c0,c2. This follows directly from moment_integrand and integration
linearity; x30_exact.py checks the integrand identity symbolically. Thus
higher-degree *parameter* polynomials in b,y only change coefficients of
the same two moment functions per panel. On a fixed physical chamber b
and y are rational/algebraic functions of the distance variables, with the
fixed square-root constants allowed. Folding swaps the panel labels and
must be respected on continuation. The general parts supplied in this
task are exactly of this parameter type. This statement does **not** apply
to introducing new q-dependent insertions inside the integral.

Over a rational-function coefficient field, multiplying the two generators
by higher parameter polynomials does not enlarge their span. It can enlarge
any *bounded polynomial-degree* search. Allowing denominators likewise
changes the finite coefficient ansatz and cannot be excluded by the d=3/4
results. A constant coefficient fit of a few parameter multiples is not a
rational-function span membership test.

An exact negative answer requires, for example, a differential-module
certificate: specify the coefficient field and analytically continued
panel module; prove its relations/normalization; show that every rational
coefficient vector fails the equations and the physical initial data for
A (including upper logs). An alternative would be an invariant
singularity/monodromy pairing that vanishes on the entire rational span
and not on A. No such invariant or complete rational solve is available
here. There is therefore **no known exact obstruction**, and membership
with general rational coefficients remains **OPEN**, rather than refuted.
The physical source -2*v0*psi30^(1) has the known order-three
presentation multiplied by a rational Coulomb factor (v0=1/r12 on the
unit-radius chart). This establishes an available exact source
representation, assuming the repository closed-form identity, without
settling either ordered summand.

The statement that ordered paths are “not natural” is an interpretation,
not such a certificate: ordered inverse/multiplication paths are perfectly
well-defined; charge polarization selects only their physical sum.

## Step (2), even assuming A+B is exact

By linearity, an exact A+B source would make the *aggregate* boundary
functional a valid input for a kernel/telescoping calculation. It would
not recover A's individual boundary functional without a further
separation argument. The classical C correction and upper-log sources
must also be accounted for when interpreting the physical C1 sector.

The missing operator is the **regular start-endpoint to logarithmic
far-endpoint source map** for

    M_z = d_z(cos(z) d_z) + 2 cos(z),
    L4 = M_X - M_Y,   I=(pi/2,3pi/2),   b=3pi/2.

For regular self-adjoint eigenmodes M h_j=lambda_j h_j, norms
n_j=int_I h_j^2, and a source R with zero start flux, its spectral
specification is

    B(Y) = int_I int_I K32(Y;X,eta) R(X,eta) deta dX,
    K32 = sum_j h_j(Y) h_j(X) h_j(eta)/(n_j h_j(b)).

The corresponding mode identity is B_j=int_I h_j(X)R_j(X)dX/h_j(b).
For the full k=4,p=0 hierarchy at E=0,

    R = (wv0/8)*Phi30 - (3w/4)*Phi41 - (w/8)*Phi42

for the selected final-v0 aggregate, with the appropriately selected
upper-log paths. One may omit these terms only after proving their
vanishing or treating them separately. A nonzero start flux adds its
endpoint term. The sign above is for +log(-cos X); chain.jl uses the
opposite logarithm convention.

**PROVED in x29 for the specified literal candidate:** the raw circle
average of Delta G restricted to X=b is not K32. The endpoint of K32 is
K32(b;X,eta)=delta(X-eta), whereas at (X,eta)=(pi,5pi/6) the raw angular
integral equals 1+sqrt(3), nonzero on a smooth off-diagonal neighborhood.
The exact radical is rechecked here. An evaluated, correctly normalized
relative-cycle boundary jump/log-residue or Riemann kernel realizing the
spectral map, with physical sheet, endpoint/corner conditions, and
convergence/interchange justification, is still missing. A raw average
cannot simply be substituted. The alternative is full coupled varied
transport with its physical framing; its free homogeneous endpoint trace
must be fixed. Knowing the order-three boundary variation alone does not
supply it. Thus step (2) has a precise *formal* source functional but no
completed evaluated-kernel route, even with A+B exact. No B_A annihilator,
elliptic factor decision, or Conjecture A decision follows.

## Additional coefficient-class searches

x30_rational.jl tests all 115 saved atoms with coefficient spaces
span{1,r1,r2,r12} and span{1,1/r1,1/r2,1/r12}, each 460 columns.
An interlaced split trains on rows not divisible by 5 and holds out every
fifth row. At threshold 1e-60, retained rank is 423 in each search. The
held-out A/B residuals are respectively 2.532e-10/9.172e-11 and
2.906e-8/6.465e-9. The very small training errors and much larger held-out
errors show poor extrapolation conditioning. These are deliberately
bounded, low-degree coefficient spaces and do not even contain all known
cubic classical controls; their failure is no certificate against rational
membership.

x30_denominator.jl instead keeps the 29 saved cl1/panel/general columns
plain and replaces the classical coefficient polynomials by P_d/r1 or
P_d/r12, d=4. It solves the cleared system divisor*Y against
[classical*P_d, divisor*29parts]. Each space contains the entire classical
d=3 + 29-parts space, since multiplication of a degree-three numerator by
the divisor has degree at most four. Each has 604 columns. Reported errors
are relative L2 on the *cleared target* divisor*Y, so they have a different
point weighting from the unweighted table. These tests allow a larger
class of reciprocal coefficients than constant simple-pole multiples, but
still do not allow arbitrary denominators or arbitrary numerator degrees.

| divisor | rank at 1e-60 | A | B | C | sum |
|---|---:|---:|---:|---:|---:|
| r1 | 500 | 7.157e-27 | 4.773e-27 | 1.456e-73 | 1.516e-61 |
| r12 | 500 | 4.191e-29 | 3.313e-29 | 2.335e-74 | 1.405e-61 |

No identity is found. These positive C/sum controls survive the enlarged
spaces, while A/B remain unresolved. The denominator run's peak RSS is
1,332,002,816 bytes, below 1.5 GB. All numerical/exact jobs were serial;
no agents or long external runs were started. Work completed in well
under the 1.5-hour allowance.

## Reproduction, ownership, and final status

Run from this workspace root, serially, with JULIA_NUM_THREADS=1 and
OPENBLAS_NUM_THREADS=1. Julia invocations use
`--startup-file=no --compiled-modules=existing --heap-size-hint=600M`.
Run x30_check.jl with default mode to regenerate n64/n80 targets; then
with X30_MODE=fit for the base fits. Run x30_control.jl for the corrected
fixed-coefficient identity. Run x30_rational.jl with X30_MODE=fit and
X30_SPLIT=interlace. Run x30_denominator.jl with X30_MODE=fit. Each output
is saved under the matching x30_* name. Exact checks run as

    PYTHONDONTWRITEBYTECODE=1 /Applications/SageMath-10-3.app/Contents/Frameworks/Sage.framework/Versions/10.3/venv/bin/python3 codex/phase13/x30_exact.py

Finally `python3 codex/phase13/x30_audit.py` checks the protected hashes,
shared-log prefixes, output/control thresholds, and recorded RSS. Its JSON
is saved as x30_audit_result.json. x30_log.md records interrupted diagnostic
runs and the final summary line. All intentional writes are x30_* files in
codex/phase13 or append-only additions to COLLAB.md and THEORY.md. No
repository source, Claude input, metadata, or commit was changed.

**Final verdict:** co-sign (i) and the corrected (ii) as numerical
verification, and (iii) as finite-search evidence only. The algebraic
polarization identity and the parameter-linearity observation are PROVED.
The interpretation that no exact A presentation is currently known is
co-signed. Its exclusion from the unrestricted alphabet is OPEN, and no
exact obstruction has been established. Step (2) remains OPEN as an
evaluated-kernel calculation even with exact A+B; the formal spectral
source map is specified above. Conjecture A remains OPEN.
