# TASK_codex_24 — physical B-density cut pilot

2026-10-24. New calculations [X]; prior conventions and regular zonal
normalizations [X21/X22]; x23 is the append-only COLLAB entry, with no separate
x23 report. Claude's files were read only. All scientific runs were serial.

## 1. Result and stopping point

**OPEN: the physical T5 pairing and Conjecture A in C1, C2, C3.** No nonzero
five-state projection, physical T5 period, or rational operator R acting on
psi40 has been certified. Neither quotient nor middle positions are excluded.
The nonzero numerical cap variations below are not T5 pairings.

This attempt constructed more of the physical cut than x22:

- **PROVED / VERIFIED:** an exact density identity reduces the sufficient
  SIX-parameter representation to FOUR parameters, with the actual charge
  numerators, both electron panels, both nuclear panels and both centres.
- **PROVED / VERIFIED locally:** the complete rational zonal connection,
  all numerator cut rows, the moving-cap equation and the ordinary face
  inventory. On every compact interior parameter chain, the specified
  physical-germ continuations give the cap z>1, including its boundary.
- **VERIFIED numerically:** based physical-germ transport at eight rational
  fibres agrees with that cap prescription, including two inactive fibres.
  A four-dimensional quadrature screen is finite and nonzero in each sector,
  but has slow convergence and no rigorous error enclosure.
- **PROVED / VERIFIED:** the q=1 electron-panel boundary has an explicit
  rational chart to the SAME bare S_e2 K3 family, with its physical relative
  interval and a normalized measure identity. The density logarithm remains;
  this is an ambient boundary identification, not a physical T5 certificate.

**The remaining obstruction is precise:** complete the limiting physical
nearby-cycle/relative-boundary map, including logarithmic and exceptional
corner strata, then compute the q=1 log-extension transgression or another
closed T5 class and its Betti pairing. The cap is a relative four-chain;
its boundary is not zero. The ordinary q=1 surface does not by itself give
its coefficient in the full mapping cone. No closed two-form arising from
a CERTIFIED full physical boundary map was obtained, so no five-state
reduction of such a form is claimed.

## 2. Conventions and actual charge data

Sphere source y=(Y1,Y2,...), target x=(u,v,w,0), w^2=1-u^2-v^2.
These u,v are angular coordinates, not Pekeris distances. T is the pure
inverse of Lambda^2-32=-4(Delta_S3+8). Let D=1+q^4, 0<q<1, and

| Panel j | P_j | b_j | d_j, with b_j=-d_j.y |
|---:|---|---|---|
| 1 | (q^4-1,2q^2) | Y1 | -e1 |
| 2 | (1-q^4,2q^2) | -Y1 | e1 |
| 3 | (-2q^2,1-q^4) | -Y2 | e2 |
| 4 | (2q^2,1-q^4) | -Y2 | e2 |

W_j^2=D-P_j.y, p_j=P_j/D, h_b^2=(1+b)/(1-b),
B=atan(q/h_b), K_j=-c0_j+c2_j*h_b^2/q^2.
The B contribution to a panel is (2/pi) int K_j B^2 dq/W_j.

For panels 1/2 factor sqrt(2) out of the actual c coefficients. The remaining
rational coefficients in the two psi30 charges are

    [Z]   c0=(b^2-1)/72,  c2=(1-b)(9b-8Y2)/144,
    [Z^2] c0=(2-b^2)/64, c2=5(1-b)Y2/32.

Multiplication by -2v0 puts them in psi40 sectors C1/C2, with a=e2.
For panels 3/4 use

    [Z^2] c0=(1-Y2^2)/18,
    c2_3=-2(1+Y2)Y1/9, c2_4=+2(1+Y2)Y1/9.

Multiplication by +2v1 gives C3 and both a=e1,-e1. No charge sectors are
summed together. Electron exchange Y1->-Y1 exchanges panels 1/2 and 3/4,
and the nuclear centres. All eight term records, coefficient rows and
prefactors are exported in x24_exact.txt and x24_cut.txt.

## 3. Exact reduction of the initial density chain

The six-parameter chain of x22 starts from

    B^2=(1-b^2) int_[0,q]^2 dr ds/(L_r L_s),
    L_r=1+b+(1-b)r^2.

Combine its two half-power Coulomb denominators and two simple L denominators
by Dirichlet weights (1/2,1/2,1,1). Write

    lambda_a=tau*t, lambda_p=tau*(1-t),
    lambda_r=(1-tau)*eta, lambda_s=(1-tau)*(1-eta),
    R=tau*(t*a+(1-t)*p)+(1-tau)*(eta*kappa(r)+(1-eta)*kappa(s))*d,
    kappa(r)=(1-r^2)/(1+r^2).

The Jacobian is tau*(1-tau), giving the sufficient six-parameter zonal
kernel F6(y)/(1-R.y)^3 and measure

    (1-tau) dq dr ds dt dtau deta /
        [sqrt(D*t*(1-t))*(1+r^2)*(1+s^2)],
    F6=-c0*(1-b^2)+c2*(1+b)^2/q^2.

Its prefactor is -8sqrt(2)/pi^2 in C1/C2 and +8sqrt(2)/pi^2 for each C3
term, with the rational coefficient convention above. Setting r=q*xi_r,
s=q*xi_s adds q^2 and makes the domain a unit cube. F6 has source degree <=4.

**PROVED new identity:** on the physical branches,

    B^2=(1-b)*q^2 int_0^1
        xi*log((1+xi)/(1-xi)) /
        [1+b+(1-b)*q^2*xi^2] dxi.                    (1)

For a=q/h_b, the right side is
I(a)=a^2 int xi*log((1+xi)/(1-xi))/(1+a^2*xi^2) dxi.
An integration by parts, subtracting log(1+a^2) at xi=1, gives

    I(a)=-int_0^1 log[(1+a^2*xi^2)/(1+a^2)]/(1-xi^2) dxi.

Its a derivative is 2atan(a)/(1+a^2); I(0)=0. Thus I=atan(a)^2.
This proof retains the endpoints. Splitting unsubtracted diagonal principal
values in the double integral would lose a corner term and is not used.

Equation (1) preserves the actual physical scalar period and its analytic
germ. It is not an asserted isomorphism of the two ambient relative objects.
It leaves one simple L denominator. Dirichlet weights (1/2,1/2,1) now give

    c=kappa(q*xi), R=tau*(t*a+(1-t)*p)+(1-tau)*c*d,
    F4=-q^2*c0*(1-b)+c2*(1+b),
    mu4=xi*log((1+xi)/(1-xi)) /
        [sqrt(D*t*(1-t))*(1+q^2*xi^2)].

The corresponding contribution is

    sign*4sqrt(2)/pi^2 int_[0,1]^4 mu4 * T[F4(y)/(1-R.y)^2]
                                  dq dxi dt dtau,    (2)

where sign=-1 in C1/C2, +1 in C3. F4 has degree <=3 and no negative q
powers. Formula (2) is understood with the regular pure inverse; cutoff
identities and source normalization are verified independently. Endpoint
removal for the CONTINUED relative chain is the additional issue in §5.

## 4. Rational connection and complete local cuts

Put n=R.R, z=R.x, H=1-n, Q=n-z^2, ell=1-z, and

    J=(gamma-theta)/sqrt(Q),
    gamma=acos(z/sqrt(n)),
    cos(theta)=(z-n)/(sqrt(n)*(1-z)),
    K=sqrt(H)*J.

The column basis (K,log ell,1,sqrt(H)) has rational connection

    d_z K=z*K/Q + H/(ell*Q) - sqrt(H)/Q,
    d_n K=(-1/(2H)-1/(2Q))*K -(z+n)/(2nQ)
                                      +z*sqrt(H)/(2nQ),
    d_z log ell=-1/ell, d_n log ell=0,
    d_z sqrt(H)=0, d_n sqrt(H)=-sqrt(H)/(2H).

Full curvature vanishes exactly. All regular zonal inverse rows are exported
modulo target polynomials of degree <=2, which have zero target monodromy.
These are KERNEL connection rows, not a four-state connection for psi40.

If U1=T[ell^-1], Ulog=T[log ell], Ux=T[ell log ell], then the degree-d source
monomial in (2) is obtained by

| d | Apply d R-derivatives to | factor |
|---:|---|---:|
| 0 | U2=(1+E)U1, E=2n*d_n+z*d_z | 1 |
| 1 | U1 | 1 |
| 2 | Ulog | -1 |
| 3 | Ux | 1 |

D_Ri=partial_Ri+2Ri*partial_n+x_i*partial_z, including its action on the
basis. The export supplies ten monomial rows and all eight physical sums.

For 0<n<1, z>1 let s=sqrt(z^2-n)>0. Upper minus lower continuation in z is

    Delta Uk=-2*pi*i*Theta(z-1)*Gk,
    Gk=Ak(n,z)+Ck(n,z)*sqrt(H)/s.                    (3)

Ak,Ck are the exported inverse rows; G2=(1+E)G1. Exact checking gives
G1=O((z-1)), Glog=O((z-1)^2), Gx=O((z-1)^3). Therefore the d R-derivatives
above introduce NO delta terms in (3). This is specific to exponent two in
(2). In the original six-parameter exponent-three representation,

    Delta U3=-2*pi*i*[Theta(z-1)*(1+E/2)*(1+E)G1
                              -delta(z-1)/(8H)].

Discarding this face term before the density reduction would be wrong.

The four-parameter top character is proportional to

    sqrt(H)/sqrt(D*t*(1-t)*(n-z^2)).

Its cover has FOUR parameters. The beta cancellation used by x22 does not
hold throughout this cube. No rational-surface exclusion for full B follows
from x22's A-model theorem.

## 5. Target transport, cap and the boundary map still missing

Fix u=1/5, base v=1/3. Continue to v=6/5 through

    upper: 1/3 -> 1/3+i/5 -> 6/5+i/5 -> 6/5,
    lower: the conjugate path.

For fixed interior parameters the regular zonal germ is analytic across
z=+-sqrt(n); those are apparent for this germ. Its cut starts at z=1.
Let z0=(t*a+(1-t)*p).x, zd=c*d.x. Then

    z=zd+tau*(z0-zd),
    tau_c=(1-zd)/(z0-zd).

The continued compact parameter chain gives the cap z>1. The tau interval
is empty if both z0,zd<=1, all of (0,1) if both exceed 1, (tau_c,1) if
zd<=1<z0, and (0,tau_c) if z0<=1<zd. The rows in x24_cut.txt are the full
rational forms on that cap, including the rational cut character.

At this target, every cap point has R2>2/3 and n>25/37. Consequently its
upper/lower orientation agrees with that in z, and the n=0 coordinate
artifacts do not occur on the cap. A 75-digit Taylor transport at eight
rational fibres, including all six distinct panel/centre geometries and
two inactive samples, verifies the incoming physical branch and (3).
Orders 50/65 agree below 2e-34; this is numerical verification, not balls.

**OPEN for the full limiting chain:** simultaneous boundary/collision
strata must be resolved and their maps computed, rather than discarded.
The oriented list is in x24_boundaries.txt. In particular xi=1 has a
logarithmic boundary, t=0,1 has the beta cover, and q/tau endpoints intersect
both n=1 and the moving cap. The fibre transports do not certify the limit
at all these intersections or its class in a five-state quotient.

An exact sample at q=1/2, xi=2/3, t=2/5, x=(1/5,6/5) gives a smooth quartic
curve in tau:

    eta_tau^2=(1-n(tau))*(z(tau)^2-n(tau)).

Its coefficients and squarefreeness certificate are exported. Thus the
remaining tau reduction is not just substitution into a rational two-form.
An elliptic slice is not a classification of the four-dimensional direct
image and does not identify its T5 factors.

## 6. An explicit S_e2 boundary face

**PROVED / VERIFIED:** on q=1 in either electron panel, p=a=e2 and the old
t variable disappears from R. Its beta integral equals pi. Write

    R=(1-tau)*c*d+tau*e2, c=(1-xi^2)/(1+xi^2),
    Dxi=1+xi^4, Nxi=(1+xi^2)^2,
    tau_minus=-2xi^2/Dxi,
    tau=tau_minus+(1-tau_minus)*T.

The second intersection with the unit circle is

    panel 1: p_minus=((xi^4-1),-2xi^2)/Dxi,
    panel 2: p_minus=((1-xi^4),-2xi^2)/Dxi.

Therefore R=T*e2+(1-T)*p_minus and H=2Nxi*T*(1-T)/Dxi.
The physical interval is

    0<xi<1, 2xi^2/(1+xi^2)^2<T<1,

with z>1 retained on the cap. Panel 1 is the standard S_e2 family with
standard q=i*xi. Panel 2 has P=-P_standard(xi); its Weierstrass equation
is identified by X->-X,Y->iY. All changes in xi,T are rational; i is only a
constant quadratic extension for identifying the standard chart.

If

    omega_K3=dxi wedge dT/sqrt(Dxi*T*(1-T)*(n-z^2)),
    M=xi*Nxi*T*(1-T)/Dxi,

then the q=1 radical measure obeys

    xi/[sqrt(2)*(1+xi^2)] * sqrt(H/(z^2-n)) dxi wedge dtau
        = i*M*omega_K3,

choosing sqrt(n-z^2)=i*sqrt(z^2-n) on the cap. Include the old beta factor
pi and the term prefactor from (2). The actual radical coefficient is the
exported C_actual evaluated at q=1. This gives a plain-text, charge-correct
boundary form, not merely the bare surface equation.

**CITED input:** x20's rank-five/irreducibility theorem applies to this bare
family. It is not recomputed here. **OPEN:** its coefficient in the physical
relative boundary map, the reduction of its LOG-weighted class, and its
pairing. M is odd under xi->-xi whereas omega_K3 is odd, so the algebraic
product without the log descends to the rational quotient. But the retained
log((1+xi)/(1-xi)) is also odd, with nonzero derivative. One cannot infer
zero in a middle constituent from this parity of the algebraic part.
Nor can one infer nonzero from finding the K3 on an ordinary boundary.

## 7. Numerical screen, verification, and handoff

The following are candidate B cap variations divided by -2*pi*i at the
specified target. They are floating-point quadratures, not enclosures or
T5 projections. All panels/centres and correct signs are included.

| Main quadrature order | C1 | C2 | C3 |
|---:|---:|---:|---:|
| 18 | -0.01537126 | 0.04069145 | -0.00238719 |
| 26 | -0.01544313 | 0.04087221 | -0.00233042 |
| 36 | -0.01545974 | 0.04090339 | -0.00233485 |

The nonmonotone refinements do not justify extra digits. No spectral
psi40 evaluator, energy solve or Claude monodromy program was run.

Independent verification includes exact curvature and source derivatives,
actual angle/atan derivatives at 55 digits, actual zonal jumps and numerator
derivatives, Dirichlet quadratures, physical-germ Taylor transport, and
bivariate finite-field jets at primes 1000033 and 1000037. The latter match
all ten monomial rows and all eight actual rows with zero mismatches,
without using symbolic differentiation. None of these checks is called a
five-state Griffiths-Dwork reduction or a nonzero T5 Betti pairing.

Plain-text handoff:

- x24_exact.txt: six/four chain, variables, physical coefficients, full
  rational kernel/density matrices and inverse rows.
- x24_cut.txt: all rational cap/numerator rows, derivation conventions,
  target and cap intervals.
- x24_boundaries.txt: fixed/moving faces, orientations and unresolved maps.
- x24_q1.txt: rational chart to S_e2, Weierstrass data, physical interval,
  measure multiplier and logarithmic/parity qualifications.
- x24_transport.json: rational fibres and actual based target paths.

No new five-state connection or R psi40 witness is supplied: the necessary
physical map to T5 has not been certified. In particular Rv from x21 is
still only a P-harmonic section of the bare module, not an operator on
physical psi40.

## 8. Cost, co-sign/dissent, and reproducibility

Measured scientific runs take seconds each: algebra ~1 s, exact cut rows
~16 s, independent quadratures ~2 s, physical fibre transport ~6 s,
q=1 chart <1 s, two-prime checks ~6 s, cap orders 8..36 ~18 s. Largest RSS
is 394,018,816 bytes (about 376 MiB), well below 2 GB. A 1 GB RSS stop guard,
600 s CPU limit, 650 s per-pilot wall stop and exclusive scientific lock
are implemented in x24_run.py. Darwin rejects lowering AS/DATA/RSS limits
and the sandbox forbids ps, so the working guard checks ru_maxrss in-process
every 0.1 s. The early failed launcher produced no scientific result.

**Second-attempt cost:** reproduction and local row checks have the measured
bounds above. The q=1 chart makes one new finite pilot concrete: compute the
log-extension transgression and the q=1/cap relative boundary maps before
projecting its meromorphic two-form. But no finite engineering bound for
that mathematical construction, full five-state pairing or a decision of A
is justified. A proposed 1–3 h continuation would be another bounded
research attempt, not a promised certificate. A modular reduction of an
arbitrary q=1 numerator would still omit the physical framing. The original
small pilot limits (10 min, 512 MB, <=1500 columns/5000 rows) remain sensible
AFTER a certified closed form is available; no such elimination is launched
merely to fill the remaining budget.

**CO-SIGN:** x21's ambient/physical distinction and conditional quotient
framing, x22's A-model result under its cited inputs and restriction to that
model, and x23's conditional submodule audit. **DISSENT from any claim that** nonzero internal residues, nonzero candidate
cap values, a bare q=1 T5 surface, or single-branch edge regularity suffices
to decide quotient/middle survival.
The full resolved relative maps and physical pairing remain OPEN.

Run from this directory, strictly one command at a time:

    PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase12/x24_run.py x24_algebra.py
    PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase12/x24_run.py x24_cut.py
    PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase12/x24_run.py x24_numeric.py
    PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase12/x24_run.py x24_transport.py
    PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase12/x24_run.py x24_q1.py
    PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase12/x24_run.py x24_modular.py
    PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase12/x24_run.py x24_cap_pilot.py
    PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase12/x24_audit.py

Only codex/phase12/x24_* and authorized append-only logs are written.
No agents, installations, external writes, commits, stencil jobs or
Claude-owned edits. The hash audit records independent concurrent/harness
changes separately; it is not an OS syscall trace. Prior source formulas
come from codex/phase3/structure.md and x18/x19/x20/x21/x22, and COLLAB
phase 11–12. No Claude executable or external repository code was imported.
