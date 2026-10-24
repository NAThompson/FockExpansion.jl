# Global finite-part a21 estimator on Claude's exported states

Codex, 2026-10-21. **NUMERICAL**, not a certified value of the exact Fock constant. No exported state was changed.

## Convention and implementation

I used each file's full-precision E_variational header, eps=1.704, Z=2, lambda=1, Psi(0)=1, physical R=sqrt(r1^2+r2^2), and the Liverts total projection. The estimator is exactly the requested

    a_global = -D/4-b/8,
    D=FP integral_0^infinity P_v[C Psi/R^2-E Psi/R] dR,
    b=-2(pi-2)/(3pi).

Source: the supplied read-only a21_global/FINITE_PART.md. The ODE derivation is reproduced in codex_minimality.md. The independent implementation is codex_a21_functional.jl; it does not include Claude's numerical source, and reads only the exported states. All term arithmetic and radial factors use 128-bit BigFloat; angular Gauss nodes and weights originate from a Float64 symmetric tridiagonal eigensolve. Three quadrature orders are compared. No numerical radial cutoff is used for the reported estimates.

Put x=(r1-r2)/s>=0, y=r12/s, h_x=sqrt((1+x^2)/2), 0<=x<=y<=1. Then u=eps*s*(y-x), v=eps*s*(y+x), w=2eps*s*(1-y), R=s*h_x, and exp(-sigma/2)=exp(-eps*s). For an exported monomial define T=c[eps(y-x)]^a[eps(y+x)]^b[2eps(1-y)]^c, k=a+b+c. (Here c in the power and c as coefficient are distinguished by context in the file format.) Let V_s=-8/(1-x^2)+1/y. Projection is

    P_v[g]=(64/pi) integral_0^1 dy integral_0^y dx
          y(1-x^2)/(1+x^2)^3 * [1-2y^2/(1+x^2)] g.

The code maps x=y*t to a square and includes the Jacobian y. It independently checks P_v[v]=1.

For k>=2, the radial D term is

    T [V_s L_(k-1,p)-E L_(k,p)],
    L_(h,0)=Gamma(h) eps^(-h),
    L_(h,1)=L_(h,0)[digamma(h)-log eps+log(2eps)].

The degree-zero, p=0 radial finite part is

    T {V_s eps[log eps+gamma-1-log h_x]
       +E[gamma+log eps-log h_x]}.

The degree-one, p=0 part is

    T {V_s[-gamma-log eps+log h_x]-E/eps}.

These follow by taking the constant Laurent coefficient at z=0 after inserting R^z=s^z h_x^z in the radial Mellin integral. Keeping log h_x is essential: finite parts in s and R differ. All exported log terms have k>=2; the code asserts p=0 in the two exceptional cases. The same angular integration computes the direct local R^2 coefficient independently of Claude's alpha/theta quadrature.

## Manufactured check before evaluating the states

For f=R^2 log R exp(-R), h=exp(-R)[(R^2-9R)log R+8-2R]. Its Mellin transform is Gamma'(z+2)-9Gamma'(z+1)+8Gamma(z)-2Gamma(z+1). The constant term is -1, so -FP/8-b/8=0 for b=1, the exact local coefficient. The Julia script checks this constant analytically before its state loop. A separate 60-digit mpmath script, codex_manufactured.py, checks actual cutoff quadrature and the Mellin transform at small z. Both approach FP=-1 (see codex_manufactured.out); dropping -b/8 instead gives the erroneous 1/8. This is a normalization check, not a manufactured helium eigenfunction.

## Results

All 24-, 48-, and 96-point tensor quadratures agree at the 12 decimal places printed below. These are quadrature-stable estimator digits for the supplied finite states, not digits of exact a21. The independent direct local projection also reproduces Claude's values.

| N | global finite part | local projection |
|---:|---:|---:|
| 6 | 0.475417555401 | 0.476547974559 |
| 7 | 0.475287168694 | 0.472239026785 |
| 8 | 0.475576940052 | 0.471589805692 |
| 9 | 0.475838495821 | 0.471720962780 |
| 10 | 0.476031741949 | 0.471957403910 |
| 11 | 0.476172374479 | 0.472163193066 |
| 12 | 0.476276691876 | 0.472323373306 |
| 13 | 0.476355920012 | 0.472445423120 |

Relative to the supplied reference 0.476747879, N=13 global is low by 0.000391959 (about 0.0822%); local is low by 0.004302456 (about 0.902%). The global discrepancy is about eleven times smaller. N=6 happens to have a much better local value than later states; neither local monotonicity nor an error bound follows from energy convergence.

## Limits of the identity for trial functions

The fixed-log basis fixes the leading logarithmic coefficient but does not enforce the exact linear cusp coefficients. For example at N=13 the exported polynomial coefficients of w and u are -0.08682944165 and 0.05984444934; exact linear cusps would require (eps-2)/(2eps)=-0.08685446009 and (eps-1.5)/(2eps)=0.05985915493. Thus one must define the trial integral by its own Mellin finite part rather than subtracting exact-solution divergent coefficients and pretending an ordinary cutoff limit exists. I use the exact Fock b in the external -b/8 term as specified. Since the trial is not an exact ODE solution, that prescription is an estimator.

In particular the simplified ordinary residual integral in the reference requires exact lower orders as well as the log term. These states do not meet that hypothesis exactly; their comparison generally calls for the residual finite part. Energy accuracy alone supplies no bound on this singularly weighted residual functional. No extrapolation or new exact digits are claimed.

Reproduce from the working directory:

    venv/bin/python codex/codex_manufactured.py
    julia --startup-file=no codex/codex_a21_functional.jl
