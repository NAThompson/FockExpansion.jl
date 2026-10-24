# TASK_codex_15 — panel 3 rational descent

2026-10-23. All calculations are small and serial. No agents, energy solves,
stencil expansions, Claude-file edits or writes outside this directory.

**PROVED:** an explicit real 55-state rational presentation exists, without
restriction of scalars. Its coefficient field is **Q(r1,r2,r12)**; neither i,
sqrt(2), pi nor a variable radical is required by the connection. The output
rows use pi, as does the defining psi30 formula. The matrices, all per-entry
numerator/denominator total degrees and sharp polynomial row clearings are
exported, not inferred from line samples.

**PROVED:** these particular 55 atoms are not all unweighted form-admissible.
An explicit quartic-denominator atom is not even locally L2. A sufficient
homogeneous angular weight delta^4 makes every atom form-admissible and
leaves its radial degree unchanged. Optimality of that weight is **OPEN**.
Thus the request for a *minimal* admissibility weight is not completely closed.

## 1. Atom definitions and exact rational connection

Set x=r1, y=r2, c=r12, Q=x^2+y^2, e=2Q-c^2>0, and

    H=c^2/e, k=2(x^2-y^2)/e,
    T(s)=s^4+k*s^2+H, U(s)=s^2+1, V(s)=s^2+H,
    A(s)=atan(s), B(s)=integral_0^s H/(H+t^2) dt,
    f(s)=sqrt(H)/sqrt(T(s)),           0<=s<=1.

In particular B=sqrt(H)*atan(s/sqrt(H)). This definition by a rational
integral makes its invariance under the sign of sqrt(H) explicit. Write
f_A,n=f*A^n/n!, f_B,n=f*B^n/n!, identifying the two zero jets.
All integrals below have the fixed upper endpoint 1. The kernel list is

    1, s, 1/T, s/T, s^2/T, s^3/T, 1/U, s/U, 1/V, s/V.

For each kernel retain jets A0,A1,A2,B1,B2, except that the kernel s
retains only A0,A1,B1. Add integral f_B,1/s, the five endpoint values
f_A,0(1), f_A,1(1), f_A,2(1), f_B,1(1), f_B,2(1), and F0=1.
That is 5+3+40+1+5+1=55 atoms. Machine order and names are in
`rational_connection.json`: J,K,Q0..Q3,U0,U1,V0,V1,FP_B1,E,F0.
These are a sufficient presentation, not a proved minimal independent basis.
Some endpoint atoms are constant multiples of others.

The exact first-order integrand identities are

    partial_H f_F,n = (1/(2H)-1/(2T))*f_F,n
                     + [F=B]*( n*f_B,n/(2H) - s*f_B,n-1/(2V) ),
    partial_k f_F,n = -s^2*f_F,n/(2T),
    partial_s f_F,n = -T'/(2T)*f_F,n + a_F*f_F,n-1,
    a_A=1/U,  a_B=H/V.

Terms with n-1 are omitted when n=0. These identities, including derivatives
of the rational kernels, determine all connection entries by the following
finite rational reduction; the script implements it without interpolation.
For a requested rational kernel r, subtract

    q' - (T'/2T)q

so the remainder has only simple T,U,V poles and polynomial part of degree
at most one. To remove N/D^2 choose q=p/D, where, modulo D,

    p = -(2/3)*N*(T')^-1       if D=T,
    p = -N*(D')^-1            if D=U or V.

Partial fractions of the remainder give the ten listed kernels. The exact
integral reduction is

    integral r*f_F,n = q(1)*E_F,n - [n=0]*q(0)*F0
                      + integral remainder*f_F,n
                      - integral q*a_F*f_F,n-1.

The last integral is reduced recursively; at most two lower jets are needed.
This proves closure, not just its evaluation at sample points. For FP_B1,
first separate the coefficient of 1/s; B/s is regular at zero. Endpoint
rows use E0=sqrt(H/(1+k+H)) and

    partial_H E_B,n = (1/(2H)-1/(2(1+k+H))+n/(2H))*E_B,n
                      -E_B,n-1/(2(1+H)),
    partial_k E_F,n = -E_F,n/(2(1+k+H)),

with the extra n term absent for A jets and no lower term at n=0. F0'=0.

`rational_connection.json` stores **dX=M X**, with exact reduced entries.
Its H matrix has 225 nonzero entries, maximum degrees (3,4); its k matrix
has 174 entries, degrees (3,3). `distance_connection.json` exports the exact
pullback M_r=H_r*M_H+k_r*M_k for each physical distance.

| derivative | nonzero entries | max numerator/denominator total degree |
|---|---:|---:|
| x | 210 | (7,8) |
| y | 230 | (8,9) |
| c | 233 | (7,8) |

Every individual degree is attached to its entry. `row_clearing.json` gives
the joint least common polynomial denominator across all three derivative
matrices, separately for every row; its maximum degree is 11. These factors
are sharp for the **exported rows and columns**, up to nonzero rational
units. They are not invariant under additional rational basis changes or
quotienting actual endpoint dependencies. Axis-specific LCMs are also in
the distance JSON. Apparent factors c^2-Q arise from U/V collisions at H=1;
the integrals extend analytically there and cleared identities apply by limits.

Panel 4 is obtained by exchanging x and y, including exchanging the x and y
derivative labels. No further state or conjugate block is needed.

## 2. Galois descent and globality

Let z_i be the four roots of T. The quartet kernel is obtained from the old
pole integrals by

    s^j/T(s) = sum_i z_i^j / (T'(z_i)*(s-z_i)),   j=0,1,2,3.

This is a Galois-invariant Vandermonde rebase, generically invertible. It is
equivalent to using symmetric root moments, with a rational change between
the resulting four kernels. Pairing the +/-i poles produces 1/U,s/U;
pairing +/-i*sqrt(H) produces 1/V,s/V. Replace the old D-angle n jet by
H^(n/2) times that jet, and multiply every old density by sqrt(1+b)=c/R.
The normalized zero jet becomes 1. These last two changes are essential:
symmetrizing only the quartic poles leaves h and normalization radicals.

The final definitions involve only real rational parameter kernels, their
positive square-root density, and fixed real intervals. All root labels,
complex conjugations, h-sign choices and normalization signs disappear.
More precisely f=c/sqrt(c^2+2(x^2-y^2)s^2+e*s^4), with c>0 and the positive
square root. The polynomial under that root is strictly positive on [0,1]
throughout the strict triangle. B also has the rational integral definition
above. Hence every atom is a unique real-analytic, single-valued function
on the strict physical triangle. “Galois invariant” concerns the auxiliary
algebraic cover/root labels with the physical density fixed; it does not
assert absence of complex analytic monodromy of elliptic integrals.
There is no moving-endpoint fold. The phase7 involution identity proves
that these fixed-endpoint output moments equal the original physical moments.

## 3. Exact moment and psi30 output rows

Define I2_F,n=integral s^2*f_F,n. Integration of d(s^3*f_F,n) gives

    I2_F,n = E_F,n - k*J_F,n - (2H-k^2)*Q2_F,n + kH*Q0_F,n
             - [n>0]*C_F,n-1,
    C_A,m = K_A,m-U1_A,m,
    C_B,m = H*K_B,m-H^2*V1_B,m.

Only m=0,1 occurs. Integrating d(f_B,2/s) gives

    Iminus = integral f_B,2/s^2
           = -2*Q2_B2-k*Q0_B2+FP_B1-V1_B1-E_B2.

Set

    F0out = pi/2*J_A0-2*J_A1+4/pi*(J_A2-J_B2/H),
    F2out = 4/pi*Iminus-pi/2*I2_A0+2*I2_A1-4/pi*I2_A2.

Then **N0=(R/c)F0out, N2=(R/c)F2out**, R=sqrt(Q).
In the existing phase5/7 convention the panel-3 contribution to R^3 psi30 is

    Z^2*c*e/18 * F0out - 2*Z^2*e*(x^2-y^2)/(9*c) * F2out.

This output is rational in distances over Q(pi), homogeneous of degree
three. The R radical in the individual N moment conversion cancels here.
`outputs.json` gives both moment rows and all 22 nonzero coefficients of
this homogeneous output. Panel 4 is its electron-exchange image.

## 4. Admissibility: counterexample and a radial-order-preserving bound

**PROVED obstruction.** Take x=1,y=2,c=1+epsilon, epsilon>0. The atom
Q0_A0 is integral sqrt(H)/T^(3/2). At epsilon=0,
T=(s^2-1/3)^2 has an interior double root. Local scaling at s=1/sqrt(3) gives

    Q0_A0 ~ (27*sqrt(3)/16)/epsilon.

For Cartesian configurations with these radial distances, epsilon~theta^2
as the angle between the electron position vectors tends to zero.
The six-dimensional volume contains theta dtheta, so |Q0_A0|^2 is not
locally integrable. This persists on open neighborhoods of the two radial
distances, not merely the single chosen radial slice. Polynomial or
exponential factors that do not vanish there do not cure it. Thus a
weight-free assertion for **all of these new atoms is false**. Singular
Galois-invariant rebasing need not preserve boundary admissibility.

**PROVED sufficient weight.** Define

    Delta=4*x^2*y^2-(Q-c^2)^2>0, delta=Delta/Q^2, w=delta^4.

This is a rational homogeneous angular weight, rather than the polynomial
W^18 previously used. On the normalized compact cone Q=1, 0<delta<=1,
and c^2,e>=delta/2, x^2,y^2>=delta/4. The quadratic form

    Tdist(s)=c^2+2(x^2-y^2)s^2+e*s^4

has positive matrix with trace 2 and determinant delta. Consequently
Tdist >= (delta/2)*(1+s^4). All distance derivatives of its polynomial
coefficients are bounded on this cone. The new density is c/sqrt(Tdist);
its values/first derivatives are bounded by C*delta^-1/2 and
C*delta^-3/2. The quartet integrands have the form

    e*c*s^j/Tdist^(3/2) * angle^n/n!,

so their kernel values/derivatives are bounded by C*delta^-3/2 and
C*delta^-5/2 before differentiating angles. U kernels are bounded. V kernels
are e/(c^2+e*s^2) and s times that expression, bounded by C*delta^-1,
with first derivatives bounded by C*delta^-2. B=integral c^2/(c^2+e*t^2)dt
satisfies 0<=B<=s<=1, with the deliberately loose derivative estimate
C*delta^-2; the same derivative bound applies to B/s. A is parameter
independent and bounded. Endpoint E0=c/(2x), its derivatives, and endpoint
angle jets obey these same aggregate bounds. F0=1. Therefore uniformly
for **all 55 atoms**,

    |X|<=C*delta^-3/2,  |grad_angular X|<=C*delta^-7/2.

Since delta has bounded angular derivatives, wX and its angular first
derivatives are bounded. Returning to physical Cartesian coordinates gives
|wX|<=C and |grad(wX)|<=C/R. Every polynomial times wX times
exp(-eps*(x+y)) is in H1(R6): the worst corner radial gradient integral is
integral_0^1 R^3 dR. All boundary strata are controlled by the same estimates.
One can justify weak derivatives by slicing away the lower-dimensional
collision/collinear sets; bounded values and gradients exclude interface
jumps. Coulomb form admissibility follows from Hardy's inequality.

Multiplication by w does not increase radial degree: R^3*wX still has value
O(R^3), gradient O(R^2), and so can carry order-three angular information.
It does modify that angular function. Exact recovery of the original output
requires dividing its output coefficients by w. Thus this result does **not**
prove that a polynomial-coefficient span of weighted atoms alone contains
the exact unweighted psi30 angular jet. Retaining its exact output as a
separate physical atom or proving admissibility for better output-adapted
combinations is a distinct step. The minimal sufficient power/weight,
and an admissible minimal rational presentation with that exact angular jet,
remain OPEN. No corner-accuracy or global recurrence completeness claim
is inferred just from w's zero radial degree.

## 5. Verification, Claude comparison, and status

`verify.py` compares the complete rational gauge connection to the independent
phase7 h,z export at three fresh H,k points: **2604 entry identities**, worst
absolute error <4.1e-56 at 55 digits. It separately checks 20 exact Hermite
certificates, all 110 H/k actual-integral derivatives (central difference,
error <1.9e-36), and both output functionals against direct physical-integrand
quadratures (error <1e-48). The explicit counterexample has
 epsilon*Q0=2.9237333953 at 1e-4, 2.9228492024 at 1e-6, tending to
27*sqrt(3)/16. No convergence from those two samples is used as a proof.

`audit.py` verifies all **673** reduced distance entries, per-entry degrees,
axis row LCMs and **1398** exact fresh rational-point pullback identities.
It exports all 55 joint sharp row clearings. `build.out`, `distances.out`,
`verify.out`, `audit.out` and `outputs.out` contain the actual serial runs.

**CO-SIGN (proved independently):** phase8's proposal that a symmetric-pole
presentation can rationalize the connection without the 220/440-state scalar
restriction. This report supplies the additional normalizations and exact
55-state matrices for panels 3/4. **DISSENT** from weight-free admissibility
of every symmetrically rebased primitive; the explicit Q0 counterexample
rules that out. This does not refute admissibility of physical N0,N2 or their
special output combinations.

The final read-only inventory currently has no `claude8/panel1_*` files.
`claude_inventory.json` records that absence. The requested conditional
panel-1 comparison is therefore **PENDING**, not passed. The compiled
classical generator and baseline numerics are Claude's work; no co-sign of
an uninspected phase8 numeric/generator result is asserted.

Deliverable status: rational connection, exact outputs, row clearings and
interior globality CLOSED; admissibility attempt gives a counterexample and
sufficient radial-order-preserving weight, optimal weight/exact enriched
angular-jet claim OPEN; Claude panel-1 comparison PENDING its artifact.
