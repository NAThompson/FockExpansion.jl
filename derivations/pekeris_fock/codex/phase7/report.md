# TASK_codex_13 — structured order-three connections

2026-10-22. All calculations serial; no agents or large linear systems.
This directory contains explicit sufficient connections, exact degree data,
a fixed-endpoint globality proof, sufficient admissibility weights and an independently
derived classical Hamiltonian subblock. Minimal actual-function ranks and
one unweighted global order-three Hamiltonian recurrence remain open.

## 1. Matrices, atoms and coefficient fields

`connections.json` is the machine-readable sparse export. Each entry is
[row, column, coefficient, [numerator degree, denominator degree]], with
**dX=M X**. State order is stored explicitly. Fractions are reduced exactly
over the indicated rational function field before measuring total degree.
Entries not listed are zero. Jets are Taylor coefficients, divided by n!.
These degrees are not common-row clearing degrees or Laguerre bandwidths.

### Classical block

For a=alpha-pi/2, p=p0-pi, w±=(a±p)/4, the 18 actual atoms are

    L(w+ +m*pi/4), T(w+), L(w- +m*pi/4), T(w-)  (m=-1,0,1),
    L(p/2), L((a+pi/2)/2), L((pi/2-a)/2),
    lc(w+), lc(w-), a,p,a²,p²,1.

The ordering of the first eight is three L then T for each sign, as the
JSON specifies. L'=s cot(s), T'=s tan(s), lc'=-tan(s).
With t=exp(ia/2), q=exp(ip/2), both matrices are rational over Q(i,pi)(t,q).
The dt matrix has 28 nonzero entries; dq has 25. Each has maximum
numerator/denominator degrees **(2,3)**. This is an exact sufficient
18-state connection; C is the output row already specified in the phase5
source. The specialized cyclic rank of that row is 12; 12 is not promoted
to a proved actual-function minimal rank. The selected 12 derivative atoms
and six omitted derivative targets are frozen in `line_inputs.json`.
Their 128-point independent line interpolants all fail holdouts; no reduced
classical connection degree is claimed from those samples.

### Each elliptic panel: explicit 55-state connection

Start with the 13 generating states J,K,J_i(i=0..7),FP,E,F0 in the exact
phase3 integrands, but choose **fixed upper endpoint 1** for the global basis.
The minus-chart matrices are additionally exported to audit the old min form. Retain the five jets 00,10,20,01,02, except:
K retains only 00,10,01; FP retains only 01; F0 retains only 00.
Thus the actual closed list has

    5 J + 3 K + 40 pole integrals + 1 FP_D + 5 E + 1 F0 = 55.

Duplicated zero jets, four identically zero F0 jets and unused K/FP jets
are removed. This is a sufficient presentation, not an irredundant basis.
No periodic/actual-integral independence theorem is asserted.

The new useful rational chart is **(h,z)**, where z is one quartic root:

    b=(h²-1)/(h²+1),
    Y²=-(z²-1)(z²-h²)/(2z²(1+h²)),   Y<0 on the physical chart,
    poles=(z,-z,h/z,-h/z,i,-i,ih,-ih),
    u=1 (h>1), u=h (h<1).

In this chart the complete connection is rational in h,z: although Y is
algebraic, the integrand depends on it only through Y². In the generating
connection replace pole velocities by their elementary h,z derivatives,
f0 by sqrt(1+h²)/(sqrt(2)h), and its logarithmic h derivative by
-1/[h(1+h²)]. This gives all entries without cyclic reconstruction.

| chart / derivative | nonzero entries | maximum degrees |
|---|---:|---:|
| u=1, dh | 282 | (6,7) |
| u=1, dz | 293 | (4,5) |
| u=h, dh | 311 | (6,7) |
| u=h, dz | 293 | (4,5) |

For comparison the unreparameterized exports, with h,z,u treated as
explicit algebraic generators, give (b derivative) 358/379 entries,
degrees (12,13), and (Y derivative) 293 entries, degrees (7,9).
These last degree measurements are in b,Y,h,z,u before algebraic reduction.

### Exact physical pullback for all four panels

Write x=r1,y=r2,c=r12, Q=x²+y², R²=Q, P²=2Q-c²,
A²=(Q-c²)²-4x²y² with A=i sqrt(-D). Then

    h1=x/y, h2=y/x, h3=h4=c/P,
    z1²=(Q-c²+A)/(2y²), z2=(y/x)z1,
    z3²=(y²-x²+A)/P², z4=i*c/(P*z3).

Choose the global roots described in phase5/goal2_global.md. For any distance
coordinate r, the explicit connection is

    M_r = (partial_r h_j) M_h + (partial_r z_j) M_z.

Every derivative here is explicit by differentiating the displayed rational
h_j and squared-root definitions: R_r=Q_r/(2R), P_r=(2Q_r-2c c_r)/(2P),
A_r=D_r/(2A), and z_r=(z²)_r/(2z). This is an exact algebraic distance
matrix formula, not a numerical chain-rule approximation. It also specifies
pullback to t,q by their distance formulas at R=1. The sufficient common coefficient field for THESE distance connections is
**K(P,A,z1,z3), degree <=16** over K=Q(i,sqrt(2),pi)(x,y,c).
R is unnecessary: h and z have the displayed squared definitions involving
only P,A and rational distances; differentiation stays in that field.
Individually panels 1,2 need at most K(A,z1), degree <=4, and panels 3,4
need at most K(P,A,z3), degree <=8. No necessity or minimality is claimed.
The earlier <=32 field remains sufficient for prior normalized output
presentations; the new <=16 statement is for the connection coefficients.

Restriction of scalars is constructive: use the at most 16 binary tower
monomials in P,A,z1,z3 and their four quadratic defining relations. Replace
each coefficient multiplication by its multiplication matrix and add the
derivatives of those monomials to the connection. If there are extra field
relations, first quotient the dependent monomials. This gives a rational
sufficient connection. We do not export its expanded matrices or quote their
degrees; the small degrees above belong to the algebraic cover, not to an
uncomputed rational restriction or a common physical row denominator.

### Reduced per-panel line connections

Claude's 200 saved exact points, prime 998244353, pi=101, along
(t,q)=(3+lambda,4+2lambda), were frozen before analysis. Independent
Euclidean interpolation uses **128 training points and 72 holdouts**.
All 438 nonzero panel boundary coefficients pass all holdouts. Exact
minimum-total-degree regular interpolants have maxima

    panels 1,2: (47,47); panels 3,4: (49,49).

The atom list is the 13 selected actual derivatives of N0,N2 specified per
panel in `line_inputs.json` (`meta`). `line_audit.json` supplies every
nonzero boundary coefficient as numerator/denominator arrays in increasing
lambda powers. Derivatives that are themselves basis atoms have unit rows;
zero omitted coefficients are zero. `line_connections.json` assembles both complete
13x13 dt/dq connection matrices **restricted to that line**. The entry
functions are restrictions of partial derivatives, not the total derivative
along the line (which would be dt+2dq).

Two additional fresh split points are checked by our own Taylor-row
construction and independent modular elimination; see fresh_line_check.json.

These are finite-field line reconstructions with holdouts, not characteristic-
zero bivariate identities or a proof of minimal rank. In particular a degree
along one line is not a certified bivariate degree bound. The 55-state
matrices above are exact symbolic sufficient connections independent of
those reconstructions. The existing 11/12/13 per-panel cyclic rank witnesses
remain specialized represented ranks; generic symbolic closure and actual-
function independence still require proof.

## 2. Interfaces and domain

`global_and_domain.md` proves that the old min-endpoint primitive atoms
require matching, but then removes the switch for the physical N0,N2:
their integrands are antisymmetric under s -> h/s, hence integral_1^h=0.
Both physical moments equal their **fixed-endpoint-1** integrals throughout
the strict triangle. Choose the 55 primitive atoms with that fixed endpoint;
they are globally single-valued analytic in the interior, including exchange
and c²=Q interfaces. Use the `plus_h`, `plus_z` matrices on both sides.
Apparent pole-collision singularities of the coefficient presentation are
handled by the actual integrals or cleared identities. No extra complement
states are necessary. The separate old-J derivative jump is a check on the
old primitive presentation, not an obstruction to this global construction.

The same note proves unweighted admissibility of all eighteen classical
atoms and proves explicit sufficient **W^18** weights for every one of the
55 panel atoms, where W is the homogeneous degree-12 polynomial defined
there. This supplies a valid form-domain enrichment, at a substantial loss
of corner order; it does not prove unweighted order-three jet preservation.

## 3. Independent classical Hamiltonian recurrence

`angle_recurrence.py` derives the invariant classical subblock {1,a,a²}
from the distance Hamiltonian, independently of Claude's phase7 files.
Let M=-32 eps³ r1 r2 r12, Q=r1²+r2² and eps=1. With s=r1+r2,

    e^s(H-E)e^-s(a^k F) = a^k H_g F
                           + k a^(k-1) B F - 2k(k-1)a^(k-2)F/Q,

where, writing x=r1,y=r2,c=r12,

    B = (2y/Q) d_x - (2x/Q) d_y
        - (x²-y²)(c²-Q)/(Qxyc) d_c
        - 2(x²-y²)/(xyQ) + 2(x-y)/Q.

The saved specialized physical recurrence uses eps=1,E=-1,Z=2,lambda=1,
and row multiplier uvw Q on the COLLAB normalization of P. Its three
distinct polynomial operators are

    D0=uvw Q P, D1=uvw M Q B, D2=-2uvw M.

For coefficient arrays F0,F1,F2 the rows are

    D0 F0 + D1 F1 + 2 D2 F2 = 0,
    D0 F1 + 2 D1 F2 = 0,
    D0 F2 = 0.

Each D is replaced by its exported Laguerre stencil. Equation-index
convention: an offset d at (l,m,n) multiplies A(l+d1,m+d2,n+d3).
The extra uvw cancels derivative tails; this is a chosen sufficient clearing,
not an optimized one.

| operator | offsets | shell shift range |
|---|---:|---:|
| D0 | 563 | [-7,7] |
| D1 | 415 | [-7,7] |
| D2 | 251 | [-6,6] |

There are six nonzero blocks, 3*563+2*415+251=2770 block-offset entries.
This is a full Hamiltonian recurrence for an invariant three-atom classical
subblock, not the complete C or an elliptic Hamiltonian recurrence. The
JSON contains every coefficient and differential operator. Nine exact
Laguerre actions, polynomial derivative-tail certificates and two actual-
angle Hamiltonian product rules pass. Neither a test-eigenvalue solution
nor energy computation is being asserted at E=-1.

## 4. Co-sign / dissent and status

CO-SIGN Claude's per-panel degree reduction as **exact sampled evidence**:
the independently reconstructed 47/49 degrees pass all 72 frozen holdouts.
CO-SIGN and strengthen Claude's numerical fold claim for N0,N2: the
involution proof establishes their fixed-endpoint analyticity exactly. This
does not prove unweighted admissibility of every higher cyclic derivative;
normal-path fits do not supply a uniform Cartesian H1 estimate.
The structured symbolic matrices also demonstrate why avoiding the full
45-atom cyclic conversion can keep sufficient connections small.
DISSENT from promoting sampled rank 13 or branch invariance to a proof of
minimal actual-function rank or a rational characteristic-zero connection.
DISSENT from global analyticity of each original min-endpoint primitive;
the explicit nonzero derivative jump disproves it. The new fixed-endpoint
involution identity proves that this obstruction disappears for a sufficient
global presentation of the actual N0,N2 moments.

At the latest read-only claude7 inventory there are probes/scans but no
Hamiltonian recurrence artifact, so the task's conditional Claude recurrence
coefficient/count cross-check cannot yet be performed. The independently
supplied three-atom Hamiltonian recurrence is ready for that later comparison.

Strict open requirements: generic minimal-rank certificates; minimal radical
necessity; bivariate rational reduced matrices; sharper/unweighted panel
admissibility preserving order-three jets; a complete cleared elliptic Hamiltonian recurrence; and Claude's later full-block recurrence cross-check. No global
finite physical order-three recurrence theorem is declared complete.
