# Global single-valuedness and coefficient fields (TASK_codex_10)

## Position on Claude's PHASE 5 statement

**Agree conditionally with the derivative-basis construction; distinguish it
from the saved panel basis.** Let Ω be the connected strict physical triangle
in positive distance coordinates. Assume that each of the actual functions
1, chi, and f=R³ psi30 is real analytic on all of Ω, including r1=r2.
Their derivatives are then global single-valued analytic functions by
definition. Simple connectivity is not needed for this implication: global
existence of these actual functions is the substantive hypothesis.

The saved finite algebraic connection on one chart supplies finite local
rational differential rank after restriction of scalars. Choose from the
actual derivatives a basis of their rational differential span as germs on
that chart. Such a finite basis exists because the span is finite-dimensional;
it can include derivatives of all three starting functions. For each basis
element g_j and distance derivation d_i, there is a rational relation

    d_i g_j = sum_k C_ijk g_k.

Multiply by a common polynomial denominator. Both sides are now global
analytic functions on Ω; the identity theorem extends their equality from
the original chart to all of Ω. This proves, under the analyticity premise,
existence of one global system of polynomial differential identities for
actual derivative atoms. It does not compute its rank, a basis, coefficients,
or the polynomial clearing degrees. It also does not prove that arbitrary
solutions of the cleared system satisfy the physical matching conditions.

This is stronger than requiring every *old panel generator* to cross its
moving-endpoint charts smoothly. Those old generators may fail to do so even
when the assembled physical function and its derivative atoms are analytic.
Thus the old interface issue is a limitation of that presentation, not a
proved obstruction to existence of a global derivative presentation.

If one instead selects a basis over an algebraic extension, one must specify
global branches for its coefficient field. A simply connected real domain
alone does not remove an algebraic branch locus crossing that domain (sqrt(x)
at x=0 is an elementary example). Local algebraicity of coefficients follows
from the original finite algebraic connection and linear algebra, not from
analyticity or simple connectivity alone. The rational construction above
avoids needing algebraic coefficient branches, at a possible rank cost.

Global analytic derivative atoms can still have inadmissible behavior at
collision boundaries. Analyticity in the interior gives no kinetic-form
estimate at those boundaries; the bare-V pole lesson remains relevant.
A polynomial differential identity consequently does not establish a Ritz
basis, form-core property, recurrence minimality theorem, or convergent
infinite Laguerre realization without further domain/convergence work.

## Reduced sufficient radical field

Write K=C(r1,r2,c), Q=r1²+r2², R²=Q, P²=2Q-c²,
A²=D=(Q-c²)²-4r1²r2²; constants include i and pi. The four pulled-back
quartic discriminants are 2A/Q, up to a branch sign. Choose one quartic root
per panel; their squares are

    z1² = (Q-c²+A)/(2r2²),
    z2² = (Q-c²+A)/(2r1²),
    z3² = (r2²-r1²+A)/P²,
    z4² = (r1²-r2²+A)/P².

These satisfy exactly

    z2 = (r2/r1) z1,
    z4 = i c/(P z3),

with compatible signs (the full root sets are unchanged by sign choices).
The remaining quartic roots in each panel are ±z_j, ±h_j/z_j. Therefore

    F = K(R,P,A,z1,z3)

suffices for the saved four-panel connection. Its degree over K is at most
32, improving the previous loose degree-128 construction. Endpoint roots,
panel h values and classical trigonometric connection coefficients introduce
no further extensions, by the formulas already checked in followup6.md.
The resulting old closure bound after restriction of scalars is at most
216*32=6912. Neither 32 nor 6912 is a minimal degree or differential rank.

No radical is proved necessary in a *rational coefficient connection*:
restriction of scalars and the actual-derivative construction can move
algebraic complexity into the atoms. In the smaller algebraic presentation,
R,P,A,z1,z3 are a sufficient list, not a list of necessarily surviving radicals.
The task's minimal algebraic rank and rational rank require further reduction.

There is also a direct interior branch choice for this sufficient field:
D=((r1+r2)²-c²)((r1-r2)²-c²)<0 in the strict triangle. Take R,P positive,
A=i*sqrt(-D), and the first-quadrant roots of z1²,z3², whose imaginary
parts are strictly positive. These are global real-analytic nonzero functions
on Ω, including the electron-exchange fold. The derived z2,z4 choices are
consistent as well. Thus these particular radicals have no interior real
branch obstruction. This does not regularize the old u=min(1,h) panel atoms,
and does not establish the analyticity premise for the physical combination.
