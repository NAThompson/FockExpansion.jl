# Explicit outer-variable connection for the psi30 panels

2026-10-22. New derivation; formulas below define every connection entry.
All statements are on a nonsingular open parameter chart, with b != 0,
|b|<1, y<0 and T positive on the integration path. Complex roots may be used;
the final moments are real. Root collision charts require continuation/limits.
This is an explicit finite closed system, **not a proved minimal system**.

## Generating integral and 13-state connection

Set h=sqrt((1+b)/(1-b)), u=min(1,h),
T(q)=(1-b)q^4+(4y²-2)q²+1+b and W=sqrt(T), positive on [0,u].
Let A(q)=atan(q), D(q)=atan(q/h),
f(q;t,s)=exp(t A(q)+s D(q))/W(q).
The normalization is f(0;t,s)=f0=1/sqrt(1+b).

The eight poles z_i are the four roots of T followed by i,-i,ih,-ih.
Writing Delta=sqrt((4y²-2)²-4(1-b²)), the four roots are
±sqrt((2-4y²+Delta)/(2(1-b))) and
±sqrt((2-4y²-Delta)/(2(1-b))). Their exponents a_i are
(-1/2,-1/2,-1/2,-1/2,t/(2i),-t/(2i),s/(2i),-s/(2i)).
Thus f_q/f=sum_i a_i/(q-z_i). Their sum is exactly -2.
No derivative of an unspecified function appears in the following system.

Use states
J=int_0^u f dq, K=int_0^u q f dq,
J_i=int_0^u f/(q-z_i) dq (i=1,...,8),
J_*=FP int_0^u f/q dq, e=f(u), f0=f(0).
Here FP means lim_delta->0 [int_delta^u f/q dq+f0 log(delta)],
or equivalently int_0^u (f-f0)/q dq+f0 log u.
All positive t/s derivatives of J_* are ordinary convergent integrals.

For x=b or y put v_i=partial_x z_i, u_x=partial_x u,
l0=partial_x log f0, kappa=l0-sum_i a_i v_i/z_i.
In particular h_b=h/(1-b²), h_y=0; u_b=h_b on b<0 and 0 on b>0;
u_y=0; l0_b=-1/(2(1+b)), l0_y=0.
For a quartic root, v_i=-T_x(z_i)/T_q(z_i), where
T_b=1-q^4 and T_y=8yq². For the other roots differentiate i,-i,ih,-ih.

Define the following linear forms in the states:
H_ij=(J_i-J_j)/(z_i-z_j) for i!=j,
H_ii=[e/(u-z_i)+f0/z_i-sum_{j!=i} a_j H_ij]/(a_i-1).
The denominator a_i-1 is invertible at t=s=0, including all angle poles.
The entire Pfaffian system is:

    partial_x J = kappa J - sum_i a_i v_i J_i + u_x e
    partial_x K = kappa K - sum_i a_i v_i (J+z_i J_i) + u_x u e
    partial_x J_i = kappa J_i + v_i H_ii
                    - sum_j a_j v_j H_ij + u_x e/(u-z_i)
    partial_x J_* = kappa J_* + sum_i (a_i v_i/z_i)(J_*-J_i) + u_x e/u
    partial_x e = [kappa+sum_i a_i(u_x-v_i)/(u-z_i)] e
    partial_x f0 = l0 f0.

These equations are rational in the explicit algebraic roots and t,s.
They follow by differentiation under the integral, partial fractions and
integration of d[f/(q-z_i)]. In particular the lower endpoint sign is PLUS
f0/z_i in H_ii. Omitting it gives a wrong inhomogeneous normalization.
The finite-part derivative uses the same fixed cutoff convention.

Take the five jets (0,0),(1,0),(2,0),(0,1),(0,2). Applying the product rule
(entrywise Taylor expansion of the rational matrices) gives a finite 65-state
padded connection with algebraic b,y coefficients. Mixed jets are unnecessary:
each pure derivative equation only involves lower derivatives in that direction.
The two copies of the zero jet agree. Numerical code uses Taylor coefficients,
so second jets are one half the second derivatives.

Initial data at ANY interior base point are the displayed convergent integrals,
the displayed finite part, and explicit endpoints, using real A,D,W. This fixes
a particular flat solution; arbitrary solutions of the connection are not psi30.
No E or a21 is fitted in this construction. Those are external constants.

## Recovery of N0 and N2

Here use derivative jets without factorials: J_A=partial_t J|0,
J_AA=partial_t² J|0 and similarly for D and other states.
Set

    Q(t,s)=u³ e - sum_i a_i[z_i K+z_i² J+z_i³ J_i].

The identity d(q³ f)/dq proves Q=int q² f dq, since 3+sum a_i=1.
With Q_A etc denoting the same jets and
H_D=(J_{ih,D}+J_{-ih,D})/2, the panel moments are

    N0 = (pi/2)J - 2J_A + (2/pi)(J_AA-J_DD)
    N2 = -(pi/2)Q + 2Q_A - (2/pi)(Q_AA-Q_DD)
         +(4h/pi)[K_D+J_{*,D}+4(y²-1)/(1-b) H_D]
         -2 T(u)e_DD/[pi(1-b)u].

Integration by parts using (W/q)'=(1-b)(q²-h²/q²)/W removes
the apparent q^-2 pole in N2. The lower endpoint term vanishes because
D(q)^2=O(q²). These agree with the seven-letter finite-word reduction in
the earlier workspace compact_elliptic_moments/RESULT.txt, but the outer
parameter connection above is a new explicit construction, not that word ODE.
A finite ODE in the dummy integration variable alone would not solve the task.

## Counts, minimality, and limitations

65 is an intentionally padded upper bound per panel. Four f0 jets vanish.
Only J_{*,D} is needed from the J_* family; its equation closes using
J_{*,D} and the zero jets J_i. K_AA and K_DD are not needed: sum a_i z_i
has zero constant term, so the recovery of Q_AA,Q_DD needs only K,K_A,K_D.
The ten families J,K,J_i therefore reduce from 50 to 48 states, plus
J_{*,D}, plus five e jets and f0: **55 states per panel suffice**.
The implementation retains the redundant 65-state system for easy auditing.

The endpoint jets are algebraic multiples of 1, endpoint A,D,A²,D².
For the four physical panels these endpoint angles are constants, alpha/2,
or half p0 and their complements, on the fixed folding/endpoint charts.
Consequently the 6 endpoint states introduce no functions outside the
classical module described in structure.md. Counting only integral states
is therefore 49 per panel, or 196 across four panels. Adding the classical
module and the existing five order-two atoms gives **at most 216 generators
over an algebraic coefficient field**, with substantial redundancy expected.
This count is a closure bound, not a minimality proof or a rational rank count.

The parameter connection denominators are explicit in the formulas:
z_i, z_i-z_j, u-z_i, u, T_q(z_i), 1-b², 1+b, and a_i-1 (jet constants).
Multiplying their products clears denominators over the indicated finite
algebraic extension. Restricting scalars to a rational field requires an
algebraic basis, including the distance radicals and quartic-root splitting
field; 216 is not the resulting rational rank. Row clearing after that
restriction would give polynomial coefficients of finite degree on a chart,
hence finite Laguerre bandwidth there, by the same algebra as phase two.

**Not established here:** a minimal combined module, its rational rank,
expanded Laguerre offsets/counts, numerical shell bandwidth, or a global
smooth atom realization across all folding/endpoint/discriminant charts.
Individual padded generators switch upper endpoints at b=0 and need not
be smooth there even when the physical combination is. One must resolve
that interface before asserting a global atom Ritz domain or one global
polynomial differential recurrence. A clearing prescription on a chart
is not a computation of those requested global recurrence counts.
