# Structured atoms: interfaces and the Coulomb form domain

2026-10-22. Distances are x=r1, y=r2, c=r12 below (the panel parameter
called y in the executable matrices is denoted Y here). Put Q=x²+y²,
R=sqrt(Q), P=sqrt(2Q-c²), D=(Q-c²)²-4x²y²<0. The strict physical triangle
has x,y,c,P>0. All roots have the branches specified in phase5/goal2_global.md.

## Exact interface matching for the old min-endpoint atoms

For a panel put T(s)=(1-b)s⁴+(4Y²-2)s²+1+b, h²=(1+b)/(1-b),
u=min(1,h), A(s)=atan(s), B(s)=atan(s/h). Each retained integral is
an integral from 0 to u of one of

    A^n/(n! sqrt(T)), B^n/(n! sqrt(T)),
    s A^n/(n! sqrt(T)), s B^n/(n! sqrt(T)),
    A^n/(n! sqrt(T)(s-z_i)), B^n/(n! sqrt(T)(s-z_i)),
    B/(s sqrt(T)).

The exact retained jet list is in connections.json. For a generic point of
b=0, Y<0 in the strict physical triangle, T is positive on [0,1] and no
pole lies on that interval. The integrands, all parameter derivatives and
the variable endpoints are analytic in a neighborhood. Thus each chart
extends analytically to the interface and the following are exact matching
conditions, not fitted observations:

* the two value traces coincide, obtained by setting h=u=1;
* the A and B jets coincide there;
* z_4=z_6=i and z_5=z_7=-i (zero-based pole labels), so the corresponding
  pole-integral traces coincide as well;
* if I=integral_0^u g(s,b,Y) ds, then

      (partial_b I)_plus - (partial_b I)_minus = -g(1,0,Y).

Here h_b(0)=1; the parameter derivatives under the integral cancel.
For endpoint jets e=g(u,b,Y), the derivative jump is -partial_s g(1,0,Y).
F0 has zero jump. Higher matching relations follow by differentiating the
two analytic integrals with their respective endpoints. This prescription
specifies every order, including the derivatives of the interface data.

In particular J's derivative jump is -1/(2|Y|), which is nonzero. Therefore
the original structured panel atoms are generally NOT analytic across the
interface. Changing algebraic root labels cannot remove this jump.
At h=1 the matrix partial fractions have apparent pole collisions; the
integral trace formulas, rather than substitution into those divided
expressions, define continuation.

The two charts per panel are sufficient and necessary for this min-endpoint
presentation (one analytic chart is ruled out by J). In physical distances
b1=(x²-y²)/Q, b2=-b1, b3=b4=(c²-Q)/Q. Consequently the joint presentation
uses at most four chambers, cut by x=y and c²=Q, not sixteen independent
choices. On exchange, panels 1 and 2 are interchanged, as are 3 and 4;
the integral labels and the algebraic pole branches are transported by
these explicit parameter maps. Equal actual integrals fix the root-label
permutation; no unrecorded sign choice is permitted.

## A fixed endpoint gives the SAME physical moments globally

The important additional identity is special to the two physical moment
integrands, rather than to every primitive atom. Set

    A(s)=pi/2-atan(s), B(s)=atan(s/h),
    g0(s)=(2/pi)(A(s)^2-B(s)^2)/sqrt(T(s)),
    g2(s)=(2/pi)((h*B(s)/s)^2-(s*A(s))^2)/sqrt(T(s)).

The symbol A here is the complementary angle used in the moment integrand,
not the uncomplemented angle in the generating jets. For s>0, the involution
s -> h/s has the exact identities

    A(h/s)=B(s), B(h/s)=A(s),
    T(h/s)=h^2*T(s)/s^4,
    gk(h/s)*h/s^2 = -gk(s)       (k=0,2).

Its fixed point is sqrt(h) and it interchanges endpoints 1 and h. Therefore
integral_1^h gk(s) ds is its own negative and is zero, for every h>0:

    integral_0^min(1,h) gk = integral_0^1 gk = integral_0^h gk.

Thus **the physical moments can both be represented with fixed endpoint 1**.
No complementary integral or extra state block is needed. This conclusion
is stronger than the generic primitive matching result and does not assert
that the old J min-endpoint atom itself is analytic.

Choose all 55 generating atoms with u=1 on the entire strict triangle. They
are globally single-valued real-analytic (complex-valued atoms are allowed):
the normalized quartic discriminant is strictly negative, so T(s)>0 on the
real path; the chosen quartic roots never meet it; the angle poles ±i,±ih
also avoid it; and the FP_D integrand is regular at s=0. Every compact
interior parameter neighborhood consequently permits analytic integration
under the fixed compact path, including all derivatives. The global root
branches were explicitly specified above. Pole coincidences in the matrix
partial fractions at h=1 are apparent: the actual integral atoms extend
analytically and their cleared differential identities extend by the
identity theorem. There is no primitive first-derivative jump in THIS
fixed-endpoint presentation, so no artificial interface flux term.

Use the **plus** connection matrices everywhere, including the physical
h<1 region, with limits/cleared identities at apparent coefficient poles.
The fixed-endpoint construction also respects electron exchange by the
explicit panel permutations; no min(r1,r2) fold is introduced. This proves
global interior analyticity and physical-moment recovery for a sufficient
structured basis, not minimality or boundary form admissibility without
weights.

Continuous matched piecewise C1 atoms have ordinary weak first derivatives
and are eligible for H1 estimates. Their first-derivative jumps create
surface terms in SECOND derivatives. Thus a global strong Hamiltonian
recurrence needs interface flux conditions or cancellation in the assembled
physical combination. Multiplying a cleared equation by an interface
polynomial does not establish these conditions.

## All eighteen classical atoms are unweighted form-admissible

Use a=alpha-pi/2, p=p0-pi, w_±=(a±p)/4. In the physical angular domain
|a|<pi/2 and -pi<p<0, hence |w_±|<3pi/8. Every L argument in the atom
list has absolute value <5pi/8; the other L arguments are in (-pi/2,pi/2).
Since L'(s)=s cot(s), its apparent singularity at zero is removable and
its values and first derivatives are bounded on these compact intervals.
Similarly T'(s)=s tan(s) and lc'(s)=-tan(s) are bounded for the arguments
used here, uniformly away from ±pi/2. Angle polynomials are bounded.

In Cartesian electron coordinates, a=2 atan2(y,x)-pi/2 and
p=2 atan2(c,P)-pi, with x²+y²=R² and c²+P²=2R². The distances x,y,c,P
are Lipschitz, with bounded Cartesian weak gradients, and differentiation
of atan2 gives |grad a|+|grad p|<=C/R. Thus every listed classical atom G
satisfies |G|<=C, |grad G|<=C/R. Polynomial times G times exp(-eps(x+y))
is in H1(R6): near zero the worst radial integral is int R³ dR; infinity
is controlled by x+y>=R. The Coulomb potential is form-bounded by Hardy,
so this is the physical Coulomb form domain. No normal-path extrapolation
is used.

## Explicit sufficient weights for every retained panel atom

This result is deliberately conservative; it does not claim the unweighted
55 atoms are admissible or the weight is minimal. Define the positive
homogeneous distance polynomial

    W = x² y² c² (2Q-c²) [4x²y²-(Q-c²)²]       (degree 12).

Then **W^18 G**, for every one of the 55 fixed-endpoint panel atoms G (and also the old matched atoms),
times any distance polynomial and exp(-eps(x+y)), belongs to H1(R6).
Here the jets are the Taylor jets used in connections.json, not higher
cyclic derivatives. Complex atoms can be replaced by real/imaginary parts.

Proof with explicit exponents: normalize R=1 and set beta=W/C0, with a
fixed sufficiently large constant C0 so 0<beta<=1. All distances and
polynomial derivatives on this compact normalized cone are bounded.
Each of x,y,c,P,sqrt(-D) is >=constant*sqrt(beta), because the other factors
in their product W are bounded above. The four physical maps give
1±b>=constant*beta, |Delta|>=constant*sqrt(beta), and
h,h^-1<=C beta^-1/2. A Cauchy root bound on T and its reciprocal gives
|z_i|,|z_i|^-1<=C beta^-1 for the quartic poles. Since
T'(z_i)=2 z_i (±Delta), implicit differentiation gives
|partial_b z_i|+|partial_Y z_i|<=C beta^-6. The angle-pole derivatives
satisfy this bound as well. Complete the square in s² to obtain
T(s)>=-Delta²/[4(1-b)]>=constant*beta. This estimate is valid on [0,1].
For a quartic root z_i, factor T(s); the other three factors have modulus
<=C beta^-1. Hence |s-z_i|>=constant*beta^4. The angle poles obey the
stronger bound |s±ih|>=h>=constant*sqrt(beta).

Consequently |G|<=C beta^-5 for J,K and all pole-integral jets. Their
parameter derivatives are bounded by C beta^-16: the worst term is
W(s)^-1 z_i'/(s-z_i)², bounded by C beta^(-1/2-6-8);
angle derivatives contribute at most beta^-2 and differentiated T^-1/2
at most beta^-3/2. Endpoint terms obey the same bound since
u<=1, |u_b|<=C beta^-3/2. For the retained FP_D jet, B(s)/s<=1/h;
differentiating this regularized integrand gives the same (much looser)
beta^-16 bound, with no logarithmic cutoff term. Endpoint and F0 jets
also satisfy these bounds. The physical normalized parameter Jacobians
are bounded, and R differentiation is supplied by homogeneity.

It follows that beta^18 G and its angular distance derivatives are bounded:
the two gradient terms are bounded by C beta^(18-16) and
C beta^(17-5). Therefore W^18 G=R^216 C0^18 beta^18 G has value O(R^216)
and Cartesian weak gradient O(R^215). The fixed-endpoint atoms are analytic in the interior; for the old version,
continuity at the interfaces excludes weak first-derivative surface masses. The resulting bounds
are integrable at the corner; polynomial growth at infinity is controlled
by the exponential. Hardy again supplies Coulomb form admissibility.

This sufficient weight sacrifices low corner orders. It proves an admissible
structured enrichment, not an order-three seed with the exact Fock jet.
Removing/reducing the weights needs sharper atom-specific estimates.
