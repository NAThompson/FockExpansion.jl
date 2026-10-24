# R³ psi30 in physical distances, with general Z,E,a21

2026-10-22. Sources read only: src/FockExpansion.jl, particularly classical(),
psi30_parts(), moment_integrand(), elliptic(); earlier workspace
compact_elliptic_moments/ and round8_classical/compact_classical.py.
The last formula is attributed/copied in classical_reference.py. This is
an audit and new connection derivation, not a claim to rediscover its classical formula.

Write distances r1,r2,c=r12, R²=r1²+r2², S=r1+r2. Work first in r1>=r2,
and use electron exchange for the other chamber. Set
alpha=2 atan(r2/r1), v=(R²-c²)/R²,
p0=acos(v)=2 atan(c/sqrt(2R²-c²)),
s=S/R, d=(r1-r2)/R, x=c/R, eta=sqrt(2R²-c²)/R,
a=alpha-pi/2, p=p0-pi, w_±=(a±p)/4.
All normalized trigonometric quantities used below are algebraic in the
distances; the positive square roots are explicit. Quarter-angle sine and
cosine are obtained with the corresponding positive/continued half-angle roots.

The exact decomposition is

    R³ psi30 = P + R³ C + R³ sum_{j=1}^4 [c0_j N0(b_j,y_j)+c2_j N2(b_j,y_j)].

P is the following polynomial in physical distances (constants may contain pi,
Catalan G, log 2, E, a21). Put ell=log 2+4G/pi and Q=R²-c²:

    P = P0+Z P1+Z² P2+Z³ P3
        +E[Z S(2R²+r1 r2)/18-c(6R²-c²)/72]
        -a21[Z S Q/2-c(6R²-5c²)/12],
    P0=c(3R²-2c²)/72,
    P1=-c(5c²-6R²)ell/72-pi*c Q/36+c(91c²-106R²)/(144pi)
       +(4S³+4S²c-32Sc²+4SR²-27c³+26cR²)/288,
    P2=-S Q ell/12+53S Q/(72pi)
       +pi*(16S³-4Sc²-28SR²-45c³+45cR²)/576
       -(503S³-375S²c+384Sc²-1212SR²-656c³+960cR²)/864,
    P3=-S(5S²-3R²)/36.

Here C is the classical expression CL_compact(alpha,p0,Z)-E0p(s,x,Z,G,log2)
in classical_reference.py, fully specified in that file. More explicitly its
function alphabet consists of the following terms, with the polynomial/algebraic
coefficients K,N,e_lam,c9,qa2,qa,qp2,qp printed in that file:

    sum_{e=±1} [4 sum_{m=-1,0,1}(K_m+Z² N_m)L(w_e+m*pi/4)
                 -4(K_2+Z² N_2)T(w_e)+e_lam log|2cos(w_e)|]
    +c9 L(p/2)/2
    +2Z²/(9pi)[(r2/R)(8r2²/R²+x²-5)L(alpha/2)
                +(r1/R)(8r1²/R²+x²-5)L((pi-alpha)/2)]
    +qa2 a²+qa a+qp2 p²+qp p.

L(z)=z log|2sin z|+Cl2(2z)/2, T(z)=-z log|2cos z|+Cl2(pi-2z)/2.
Their derivatives are z cot z and z tan z. Thus all eleven L/T functions
in the formula, two log|2cos(w)| functions, a,p,a²,p²,1 form an
**18-generator closed algebraic differential module**. Each argument is an
affine combination of a,p,pi, and its differential and tangent/cotangent
are algebraic functions of distances on a branch chart. This is a closure
count, not a proof these eighteen functions are independent.

## Elliptic coefficients and arguments

Use the dimensionless x,s,v above. Four panels:

    (b1,y1)=((r1²-r2²)/R², -c/(sqrt(2)R)),
    (b2,y2)=(-b1, -c/(sqrt(2)R)),
    (b3,y3)=(-v,-r1/R), (b4,y4)=(-v,-r2/R).

For j=1,2,

    c0_j=-sqrt(2)Z[9Z b_j²-18Z-8b_j²+8]/576,
    c2_j=-sqrt(2)Z(1-b_j)[90Z y_j²-45Z-18b_j-32y_j²+16]/288.

For j=3,4,

    c0_j=Z²(1-b_j²)/18,
    c2_j=-2Z²(1-b_j)(2y_j²-1)/9.

For each pair, h=sqrt((1+b)/(1-b)), u=min(1,h),
W(q)=sqrt((1-q²)((1+b)-(1-b)q²)+4y²q²),

    N0=(2/pi) int_0^u [(pi/2-atan q)²-atan(q/h)²] dq/W,
    N2=(2/pi) int_0^u [(h atan(q/h)/q)²-q²(pi/2-atan q)²] dq/W.

The q=0 values are removable limits. W is positive in the physical interior.
The new explicit (b,y) connection and normalization data are in connection.md
and panel_connection.py. Composition with the four displayed rational/algebraic
maps and multiplication by R³ gives a finite distance-coordinate connection.
The chain rule is explicit: d g/dr_k=C_b b_{r_k}g+C_y y_{r_k}g;
multiplying an atom by R³ adds 3 R_{r_k}/R to its diagonal.

## Combined closure bound and meaning of “minimal”

The previous five order-two atoms contain alpha and V=h_old(beta+pi/2),
beta=asin(v). Since p=-pi/2-beta, a and p already lie in that module over
an algebraic coefficient field (h_old here is the phase-two radical, not the
panel h). Of the eighteen classical generators only fifteen are potentially
new. Their union with the old five is therefore bounded by twenty.
The 49 integral generators per panel described in connection.md give the
combined chartwise bound 20+4*49=216. This is an upper bound only.
No cyclic determinant or differential-Galois independence proof has been obtained
for this combined set. **The task's requested minimal closed set is unresolved.**
An order-two rational rank-five minimality theorem does not extend automatically
after enlarging the coefficient field to algebraic functions, nor after adding
these elliptic integrals.

## Verification and precision

verify_connection.py compares both parameter connections of all thirteen
padded states, all five jets, at b=.3,y=-.6 and b=-.27,y=-.7 (the shared zero
jet is checked twice). It also reconstructs the original N0,N2 integrals.
verify_repo.py reconstructs all three psi30 contributions and compares with
repo_values.jl, a startup-free, non-precompiling read-only repository invocation.
The classical formula is inherited, while the panel state/connection code is new.
The repository comparison includes both r1/r2 orderings and Z=0, general
E/a21 rather than only the helium values. Saved logs, not unrun formulas,
are the evidence. These are interior numerical checks, not endpoint/domain proofs.
