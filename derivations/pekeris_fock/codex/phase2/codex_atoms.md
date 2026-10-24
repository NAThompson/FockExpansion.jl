# Independent atom audit through radial order four

Codex, 2026-10-21. Read-only inputs: `src/FockExpansion.jl`,
`src/fourth_order.jl`, `src/feynman_k4.jl`, and `derivations/k4/README.md`
in `/Users/nathompson7/FockExpansion.jl`. No `claude2/` or `julia/` input.
PROVED means an identity derived here; VERIFIED means a finite computation;
SOURCE means transcribed from the indicated local implementation, not a new
proof of its Fock equation. No external literature theorem is inferred from
an implementation comment.

## Coordinates and what “polynomial part” means

Set a=r1, b=r2, d=r12, q=a²+b²=R², S=a+b, D=q-d²,
H=sqrt(2q-d²), and h=H/d. In the open physical triangle all are real and
H,d>0. Put

    s=2ab/q, c=(a²-b²)/q, v=D/q, g=dH/q,
    xi=d/R, sigma=S/R, alpha=2 atan(b/a), beta=asin(v).

Fold a>=b (0<alpha<=pi/2) when using chi and the classical psi30 formula;
exchange the electrons for the other chart. All assertions about
differentiation are on chart interiors, with boundary extensions requiring
limits. R is the physical hyperspherical radius, not sigma.

A polynomial/non-polynomial decomposition is not canonical: one may move
any polynomial between the pieces. The convention below exposes every
explicit polynomial in the inspected formulas and leaves the named special
functions/inverse operators in the remainder. Their harmonic subtractions
may themselves contain polynomial pieces. This is not a uniqueness claim.

## Orders one and two

**PROVED by substitution into the source.** With kappa=(pi-2)/(3pi),

    R psi10 = d/2-Z S,
    R² psi21 = -Z kappa D,
    R² psi20 = (1-2E)q/12 - Z S d/3
                +Z²(q/3+ab)+(a21-Z(pi+4)/(9pi))D + Z q chi.

The first two expressions are polynomial. The last expression has the
single non-polynomial named atom chi; even its algebraic radical part need
not be a polynomial. Thus “all terms except Z R² chi are polynomial” is
correct, but does not say that chi is an irreducible differential atom.

Here is an explicit distance-coordinate closed formula for chi, with no
implicit angular evaluator. Define

    L=log((S+d)/R),
    A=log(abs((g-c)/(g+c))),
    B0=log(abs((s+v)/(s-v))),
    J=log((c+S*d/q)² / ((1+s)*(g+c))),
    Delta=Cl2(alpha-beta)+Cl2(pi+beta-alpha)
             -Cl2(pi-alpha-beta)-Cl2(alpha+beta),
    T=beta*A+alpha*B0-Delta,
    U=T+pi*J,
    h0=-(3pi+10-16*Catalan)/(24pi).

Then

    chi = (g+s-2*S*d/q)/6 - v*L/3
               +(c*U+2*g*beta)/(6pi)+h0*v.                 (1)

**PROVED:** this is the repo expression with logarithms combined, using
s²-v²=(g-c)(g+c). Absolute values are necessary in the separate logs;
the combined J argument is positive in the folded physical chart.
At s=|v| separate terms may diverge while their combined limit exists.
Use the original limiting expression there. Cl2(t)=Im Li2(exp(it)).

Differential simplification, essential to part B:

    d Delta = B0 d alpha + A d beta,
    d T = beta d A + alpha d B0.                            (2)

These follow from Cl2'(t)=-log|2 sin(t/2)| and half-angle identities.
No infinite alphabet of successively differentiated dilogarithms is needed.

The quadratic-fit residual has a structural explanation: q chi has
non-polynomial angular dependence, and changing log R to log sigma adds
the additional angular term (-Z kappa D) log(R/sigma). A quadratic fit to
the log-free sector cannot absorb these functions. The earlier RMS
6e-3 value is Claude's evidence, not a value reproduced by this derivation.
It is not itself a proof of non-polynomiality; equations (1)-(2) display
the actual obstruction. Multiplying by exp(sigma/2) also adds polynomial
lower-order terms and does not remove the obstruction.

## Order three

**PROVED substitution:**

    R³ psi31 = Z(pi-2)d(5d²-6q)/(36pi)
                         +Z²(pi-2)S D/(6pi).

This is polynomial despite its angular square roots.

For psi30 the source gives `classical + state + elliptic`. An explicit
polynomial piece is P30=R³(e0+Z e1+Z² e2+Z³ e3+state), where after
homogenizing the source's e's it is

    P0 = d(3q-2d²)/72,
    P1 = -d(5d²-6q)(log(2)+4Catalan/pi)/72
         -pi*d*(q-d²)/36+d*(91d²-106q)/(144pi)
         +(4S³+4S²d-32Sd²+4Sq-27d³+26dq)/288,
    P2 = -S*(q-d²)(log(2)+4Catalan/pi)/12
         +53S*(q-d²)/(72pi)
         +pi*(16S³-4Sd²-28Sq-45d³+45dq)/576
         -(503S³-375S²d+384Sd²-1212Sq-656d³+960dq)/864,
    P3 = -S*(5S²-3q)/36,
    Pstate = E*[Z*S*(2q+ab)/18-(6q-d²)*d/72]
              -a21*[Z*S*D/2-(6q-5d²)*d/12],
    P30 = P0+Z*P1+Z²*P2+Z³*P3+Pstate.

**SOURCE, explicit remainder prescription:** define normalized
r=a/R, t=b/R, x=d/R, eta=H/R, delta=(a-b)/R,
alpha as above, p0=acos(v), A0=alpha-pi/2, P0angle=p0-pi.
For each sign e=+1,-1 set z=(A0+e*P0angle)/4 and eta_e=e*eta.
Use exactly the following finite sum (all K,N,elam are the elementary
polynomials displayed in `src/FockExpansion.jl`, functions K, Ncoef, elam):

    Cspecial = sum_e [4 kk(-1)Lfun(z-pi/4)+4 kk(0)Lfun(z)
                     +4 kk(1)Lfun(z+pi/4)-4 kk(2)Tfun(z)
                     +elam(delta,sigma,x,eta_e,Z)log|2cos z|]
      +c9*Lfun(P0angle/2)/2
      +2Z²/(9pi)*[t*(8t²+x²-5)Lfun(alpha/2)
                  +r*(8r²+x²-5)Lfun((pi-alpha)/2)]
      +qa2*A0²+qa*A0+qp2*P0angle²+qp*P0angle,

where kk(m)=K(m,delta,sigma,x,eta_e,Z)+Z²*Ncoef(m,...),

    c9=-eta*[Z*(2x²-5)/(9pi)-5Z²*(x²-1)/(2pi)],
    qa2=Z*x*(2x²+1)/(144pi)-5Z²*x*(x²-1)/(32pi),
    qa=7Z*x*delta*sigma/144
                         +Z²*(x*delta*sigma/64+delta*(x²-1)/(9pi)),
    qp2=Z²*sigma*(4sigma²+x²-9)/(36pi),
    qp=-Z*eta*(5x²-4)/(36pi)+Z²*sigma*x*eta/(18pi),
    Lfun(z)=z log|2sin z|+Cl2(2z)/2,
    Tfun(z)=-z log|2cos z|+Cl2(pi-2z)/2.

Then R³ psi30=P30+R³(Cspecial+Ell). This is a closed classical
remainder plus *defined definite integrals*, not a claimed classical closed
form for Ell. For reproducibility its two basic moments are

    Nj(b0,y) = (2/pi) integral_0^min(1,sqrt((1+b0)/(1-b0)))
         nj(z) / sqrt((1-z²)*((1+b0)-(1-b0)z²)+4y²z²) dz,
    n0=atan(1/z)²-atan(z/hp)²,
    n2=(hp*atan(z/hp)/z)²-(z*atan(1/z))²,
    hp=sqrt((1+b0)/(1-b0)).

Ell is the four c0*N0+c2*N2 combinations in source `elliptic_fixed`:
for b0=+c,-c and y=-xi/sqrt(2), use

    c0=-sqrt(2)*Z*(9Z*b0²-18Z-8b0²+8)/576,
    c2=sqrt(2)*Z*(b0-1)*(90Z*y²-45Z-18b0-32y²+16)/288.

For b0=-v and y=-r,-t, use

    c0=Z²*(1-b0²)/18,
    c2=2Z²*(b0-1)*(2y²-1)/9.
This specifies a finite collection of parameterized elliptic moments;
the implementation evaluates them numerically. It is inaccurate to call
psi30 wholly numerical: its classical and state pieces are closed.

## Order four

Define the polynomial harmonics

    H40=4(a²-b²)²-q²,
    H42=(3D²-4a²b²)/2,
    W=q²-8a²b²+2D²=(H40+4H42)/3.

**PROVED:**

    R⁴ psi42 = Z²(pi-2)(5pi-14) W/(180pi²).

**SOURCE / PROVED homogenization:**

    R⁴ psi41 = Z(pi-2)q*[3(32E-15)q-8(12E-5)d²]/(2880pi)
               -Z³(pi-2)(4q+10ab)D/(120pi)
               +c0 H40+c2 H42 + Z² q² Q41,

where c0,c2 are the repo's `psi41_harmonic` coefficients, including
Z[(3pi-8)E/(540pi)-(5pi-14)a21/(90pi)] times (1,4).
The pure-Z constants are stored numerical spectral results, not proved
closed constants. Their Z,Z²,Z³ triples are

    c0: (-0.0013554801797822578866, 0.0030593204017684670484,
           0.00070190494196250615738),
    c2: ( 0.0051053947955489353034,-0.0053679564244078672069,
          -0.0019003487766684624263).

Q41=psi41_z2(alpha,theta) is the numerical one-dimensional
Feynman integral in `src/feynman_k4.jl` (or the independent Green inverse).
Its mathematical definition is the harmonic-free inverse of Lambda²-32
on the source in `psi41_z2_source`. The source's harmonic component is
implicitly projected away by that inverse. No classical closed form for
Q41 is supplied by the inspected repository.

**SOURCE, psi40:** let G4=(Lambda²-32)^+ be the harmonic-free angular
inverse with the repo's regular boundary prescription. Then exactly at
the defining-equation level

    psi40 = G4[12psi41+2psi42-2Vpsi30+2Epsi20]
                  +a40 Y40+a42 Y42.

Here V=1/xi-Z*(1/r+1/t) is the dimensionless angular potential.
The inverse is a completely specified definite integral, not an opaque
name: for X=(cos alpha,sin alpha cos theta,sin alpha sin theta,0) in S³,

    (G4 f)(X)=-1/4 integral_(Y in S³) K(arccos(X dot Y)) f(Y) dOmega(Y),
    K(gamma)=[(gamma-pi)cos(3gamma)+sin(3gamma)/6]/(4pi² sin gamma),
    dOmega=sin²(alpha')sin(theta') d alpha' d theta' d phi'.

The same inverse gives

    Q41=G4[kappa*((1/r+1/t)*(5xi³/6-xi)
                                  +sigma*(xi-1/xi))].

It projects out the fourth-order harmonics. Thus it is valid to omit
the 24psi42 harmonic part of that source here. This supplies a closed
integral description for both hard fourth-order remainders while making
clear that no classical closed form has been established.

A specified exposed polynomial piece is

    P40 = (12E²-11E+6a21+1)q²/1152
          -(72E*a21+E-30a21-2)qD/720
          +Z⁴[ab*q/18+7q²/288+2H40/(135pi)]
          +a40 H40+a42 H42.

The non-polynomial remainder is q² times that displayed inverse minus
the first three polynomial terms of P40. This definite Green-integral
definition is explicit but not a finite classical closed form. The
repository implements the whole psi40 numerically (Float64). Its Z=0
and pure Z⁴ formulas above are SOURCE identities checked in the repository,
not independently proven here. The middle Z powers include chi, the
psi30 classical functions and elliptic moments under another inverse.

## Atom inventory and limitations of “minimal”

For *representing* order two one named transcendental atom chi suffices.
For rational differentiation use the five degree-zero functions

    1, L, alpha, V=h*(beta+pi/2), U.                       (3)

Indeed the h term in (1) combines with h*beta, giving

    chi=(s-2*S*d/q)/6-v*L/3+(c*U+2*(d²/q)*V)/(6pi)+h0*v.

Part B gives a finite rational connection. `minimal_rank.md` proves rank
five using actual-function monodromy and the exact generation witness. R itself is
not needed for its connection: d log R=dq/(2q), and L is scale invariant.
Slog=log(S/R) is an additional degree-zero log when converting
the *full* expansion to log sigma; log(sigma/R)=log(2eps)+Slog.
Part B delivers a six-atom audit system that retains chi and Slog
separately, and the smaller five-atom physical system obtained by replacing
L with Ltilde=L-(pi-2)/pi*Slog. The latter incorporates the conversion
without retaining two logarithms that occur only in one fixed combination.

Through order four a concise inventory of new *types* is: algebraic
radicals (R and H); logarithms and inverse angles; Clausen/Lfun/Tfun at
the finite arguments above; the two moment families N0,N2 at four panels;
Q41; and the G4 inverse applied to the psi40 source. The last two are
integral-defined angular atoms, not known finite elementary/dilog atoms.
Products of angles/logarithms already occur at order three.

There is **no justified globally minimal finite classical atom list
through order four** in these inputs. Treating an arbitrary numerical
evaluator as “one atom” makes minimality tautological and says nothing
about differential closure. In particular, neither the numerical Q41
nor the Green-defined psi40 proves or disproves D-finiteness of the
order-four content. That question remains open in this audit.

## Elementary coefficient polynomials in the psi30 remainder

To make the preceding prescription self-contained, K(m,D0,S0,X,H0,Z)
first rotates (D0,S0)->(-S0,D0) exactly m mod 4 times, and then equals

    Z*[D0*(-9S0²+6S0*X+8X²-8)-H0*(5-2X²)]/(144pi)
                  -5Z²*(X²-1)*(D0+H0)/(32pi).

Ncoef uses the *unrotated* arguments and m mod 4:

    Ncoef=D0*nu/(36pi),
    nu=X²-S0²                         (m=0,2),
    nu=9S0²+X²-2-8S0*H0               (m=1),
    nu=9S0²+X²-2+8S0*H0               (m=3).

Finally

    elam=Z*(7S0*D0*H0+9S0³+8S0*X²-26S0+10X³-12X)/72
         +Z²*S0*(137D0*H0+64S0²-154X²+26)/288.

These names distinguish the coefficient arguments from physical D=q-d²
and from the source's use of x for the normalized interelectron distance.
