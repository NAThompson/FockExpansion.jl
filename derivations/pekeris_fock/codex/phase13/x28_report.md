# TASK 28 — branch kernel, elliptic reduction, and telescoping feasibility

2026-10-24. Serial, no agents, no repository edits or commits. Verdict:
**co-sign the practical fallback; qualify the proposed identification of K32.**
B_A's annihilator and Conjecture A remain OPEN. In particular this report does
not claim to have evaluated K32. The explicit elliptic formula below evaluates
the rotation average of the pure-inverse branch variation; proving that it is
the specified regular endpoint spectral kernel requires another step.

## B. Exact unreduced variation, with all normalizations

Use x18's sphere measure vol(S3)=2*pi^2, gamma=acos(z), and

    G(z) = ((gamma-pi)*cos(3gamma)+sin(3gamma)/6)/(4*pi^2*sin(gamma)),
    T(z) = -G(z)/4,  T=(Lambda^2-32)^+.

Start with -1<z<1, sin(gamma)>0. One counterclockwise circuit around z=1
not enclosing -1 takes gamma to -gamma. Consequently

    Delta G(z) := G_after-G_before
                = (4*z^3-3*z)/(2*pi*sqrt(1-z^2)),
    Delta T(z) = -(4*z^3-3*z)/(8*pi*sqrt(1-z^2)).

Here Delta means the actual difference, NOT division by 2*pi*i. If
N=(M-1)/(2*pi*i), then N G=(4*z^3-3*z)/(4*pi^2*i*sqrt(1-z^2)).
The same single circuit with opposite orientation has the same algebraic
sheet exchange; the second circuit restores G. Thus this zonal square-root
variation is not itself a unipotent logarithmic residue operator. A log
coefficient of an integrated physical source can result from the moving
integration chain. One must not identify these two operations without a
relative-cycle calculation. The finite resonant projector cannot change this
variation: its zonal polynomial is single valued.

Derivation: sin(3gamma)/sin(gamma) is unchanged under gamma -> -gamma;
((gamma-pi)cos(3gamma))/sin(gamma) changes by 2*pi*cos(3gamma)/sin(gamma).
A 2048-step continuous complex-acos lift at z=.37 checks the difference to
2.6e-17 in Float64. This is also a check of the sign of Delta, independently
of the code's -log(-cos X) convention.

## A. Explicit rotation average: Appell F1 / elliptic, and its exact scope

Set C=A+D*cos(phi), c=A+D. On the branch with -1<A-D<A+D<1,

    J(A,D) = int_0^(2*pi) (4*C^3-3*C)/sqrt(1-C^2) dphi.

Put x=-2D/(1-c), y=2D/(1+c), and define

    p0=4c^3-3c,
    p1=-24c^2 D+6D,
    p2=48c D^2,
    p3=-32D^3.

Then the closed hypergeometric formula is

    J = 2*pi/sqrt(1-c^2) * sum_(r=0)^3
          p_r * (1/2)_r/r! * F1(r+1/2; 1/2,1/2; r+1; x,y).

Therefore the angular integral of Delta G is J/(2*pi), and that of Delta T
is -J/(8*pi). No extra sphere-volume factor is implicit. Continuation to
other charts requires continuation of the displayed roots and integration
cycle, rather than evaluation on unrelated principal branches.

Proof: t=(1-cos(phi))/2 transforms the circle integral to twice the
integral on [0,1]. Its denominator becomes

    sqrt(1-c^2)*sqrt(t*(1-t)*(1-x*t)*(1-y*t)).

Expand T3(c-2D*t)=sum p_r*t^r, then use the Euler beta integral for each
term. The beta prefactor is pi*(1/2)_r/r!. The corresponding curve is

    W^2=t*(1-t)*(1-x*t)*(1-y*t).

It has four distinct branch points generically, hence genus one. In the
s-coordinate used in x18, an equivalent branch configuration is

    w^2=(s^2-4)*((s-2A)^2-4D^2),
    branch points -2,2,2(A-D),2(A+D).

The t and s descriptions have the same cross-ratio up to permutation.
For example lambda=((1+A-D)*(1-A-D))/((1+A+D)*(1-A+D)),
and j=256*(1-lambda+lambda^2)^3/(lambda^2*(1-lambda)^2).
The polynomial numerator produces complete Abelian integrals on this curve;
this is an explicit elliptic reduction, without needing to choose a K/E/Pi
basis. At A=.2,D=.15 the direct angular integral is
-3.5299025882181247 and the Appell expression is -3.5299025882181261.

**Scope correction:** this proves a closed elliptic/hypergeometric expression
for the branch-variation average, not a closed evaluation of K32. Genus one
of a presentation does not alone prove non-elementarity of the resulting
physical B_A, or even irreducibility of a special numerator period. Nor is
the full x18 G average identical to its Delta G average: the former retains
acos and has the regularized representation described there.

### What was checked for K32 itself

With a=pi/2, b=3pi/2 and regular eigenfunctions of
M=d(cos(z)d)+2cos(z), define K as in x27. Termwise, in a regularized or weak
spectral interpretation,

    M_Y K=M_X K=M_eta K,   K(b;X,eta)=delta(X-eta).

Also K is symmetric in all three arguments. In particular its required
weak moments, valid wherever the defining spectral map is justified, are

    int K(Y;X,eta)*sin(eta) deta = -sin(Y)*sin(X),
    int K(Y;X,eta)*(1 +/- sqrt(2)cos(eta)) deta
       = (1 +/- sqrt(2)cos(Y))*(1 +/- sqrt(2)cos(X)).

These use h0(b)=-1 and h_pm(b)=1. They are useful normalization tests for
ANY proposed elliptic formula, and stronger than fitting one point.

The read-only Claude SOps implementation was used at n=24,40,56. For its
nodal eigenmatrix V and endpoint/interpolation rows e_b,e_X,

    P_X=V*diag((e_X V)/(e_b V))*V^{-1}

is the discrete source-to-output kernel map. It reproduces all three exact
moments with maximum errors 1.24e-14, 1.94e-13, 4.12e-13, respectively;
the endpoint identity errors are <=2.3e-15. At (Y,X,eta)=(3.0,3.4,2.8),
point-kernel interpolation gives .9640,.9344,1.0933. Those values are NOT
converged continuum values or evidence against existence of an elliptic
formula. They demonstrate why a raw truncated delta-kernel point comparison
would not certify the proposed identification.

Missing: a map from the continued sphere cycle to this particular spectral
product, including its characteristic support, corner distributions, measure,
and endpoint normalization. Equations and a boundary delta alone should not
be substituted for that construction without an appropriate uniqueness
argument. I therefore neither confirm nor refute the asserted equality of
K32 with a specific rotation-reduced Riemann function. Its closed-form class
remains OPEN here. The displayed J is a concrete candidate building block,
not a secretly claimed solution of that missing problem.

## C. Creative telescoping: concrete input and stopping condition

No executable small-hours route to an exact B_A operator is established by
the currently supplied data. This is a practical assessment, NOT a proof
that no such route exists. In particular:

1. The two-variable route first needs an exact finite connection / holonomic
   ideal for the UNVARIED interior phi3^A, then tensors it with the algebraic
   factor wv0/8 and an evaluated K32 module. Restrict to rational/algebraic
   coordinates, integrate both strip variables with endpoint terms, include
   any actual upper-log hierarchy contributions, then restrict Y to arc b.
   The known tangential A3 operator cannot supply normal derivatives or the
   analytic part of phi3^A. The one scalar elliptic curve of J does not fill
   that gap.
2. A route need not literally compute an ideal for phi3 first: insert its
   explicit lower-order panel integrals and telescope the resulting iterated
   integral directly. The x24 four-parameter log-density/cap construction
   supplies a concrete starting presentation for relevant panel terms, but
   was not already a complete B_A relative-chain certificate. It still
   requires the correct chosen-ladder numerator, continuation, all classical
   pieces and face/corner distributions. A telescoper for an interior density
   alone gives an inhomogeneous boundary identity, not a B_A annihilator.
3. More than 150 coefficients/unknowns in an unsuccessful bounded arc guess
   is NOT a proved lower bound on the minimal physical annihilator, and
   does not bound multivariable holonomic rank. A better gauge, coordinate,
   block connection or telescoping order could change complexity radically.
   The stronger wording 'minimal operator exceeds ~150 unknowns' in the
   prior THEORY summary is therefore only heuristic, not certified.

A concrete staged attempt would first close the exact source connection
from the known panel blocks, retain separate charge and ladder sectors,
then telescope one panel variable and compare all boundary terms against
that representation. Only after this finite input exists would an
ore_algebra computation of an integration ideal, or a rational connection
implementation with certificates, be reviewable. No such input/certificate
has been produced in this task; no ore_algebra run or guessed annihilator
is reported. Sage is not on the current PATH; this is secondary to the
missing mathematical input and is not asserted to mean Sage is uninstalled.

**Co-sign:** no justified expectation of finishing the actual B_A telescoper
within a few hours from these inputs; use the explicit-stall fallback.
**Dissent/qualification:** neither ellipticity of K32 nor impossibility of
all short telescoping routes follows from the feasibility argument. B_A's
operator, its elliptic constituent and Conjecture A remain OPEN.

## Verification

Reproduce with JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia
--startup-file=no --heap-size-hint=600M codex/phase13/x28_check.jl.
Only the two scalar quadratures, complex branch lift, and weak SOps kernel
checks are certified by this script. Input hashes, append-only prefixes,
checker stdout and resource output accompany the report. No pointwise
continuum K32 error enclosure or physical B_A value was calculated.
