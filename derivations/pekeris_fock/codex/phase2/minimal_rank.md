# Minimal differential ranks: a branch/monodromy argument

Codex, 2026-10-21. Rational coefficient field K=C(a,b,d); constants
include pi and the state parameters. Work first on a regular folded
physical chart, and analytically continue a putative rational relation
to the complex algebraic cover H²=2q-d². Loops below can be lifted or
repeated so that all algebraic functions return to their initial sheet.
Rational coefficients are single valued along every such loop.

## Independence, not just a formal matrix rank

Let X=c+i*s=(a+i*b)²/q and Y=g+i*v. On local branches,
alpha=log(X)/i and beta=log(Y)/i, up to constants. X and Y are nonzero
where q!=0; (g+i*v)(g-i*v)=1. Set

    A0=(g-c)/(g+c), B1=(s+v)/(s-v).

The already-derived identity is

    dU=beta*dlog(A0)+alpha*dlog(B1)+pi*dJ.             (1)

**PROVED:** U cannot be a K-linear combination of algebraic functions
and the single logarithms L, alpha, h*beta, Slog. To see this, take a
loop gamma about d=infinity with a,b fixed generic. On the algebraic
cover Y has a zero or pole of order two there, so gamma changes beta
by a nonzero multiple of 4pi, while alpha is unchanged. In (1) all
algebraic coefficients return to themselves. Consequently

    variation_gamma U = C*log(A0)+constant, C!=0.

Now take delta around the point d=a-b on the sheet H=a+b. There g=c,
so A0 has a simple zero (generic a,b). Therefore

    variation_delta variation_gamma U = C*2pi*i != 0.

J is a single logarithm of an algebraic argument, so its contribution
to this double variation vanishes. The same is true for all the other
listed single-log functions: their first variations are algebraic and
the second lifted variation fixes them. A rational coefficient times U
therefore cannot cancel against a combination of these lower functions.
This argument uses continuation of actual functions, not an assigned
formal transcendental weight.

For the remaining single logarithms, independence over K follows by
successive loops:

1. At d=-(a+b), L=log((a+b+d)/sqrt(q)) has a nonzero logarithmic
   variation. Y is nonzero and regular there generically, and alpha,
   Slog are independent of d. Thus the coefficient of L must vanish.
2. The preceding d=infinity loop then isolates h*beta: its change is
   C*h!=0, whereas alpha and Slog remain fixed. Its coefficient vanishes.
3. If present separately, Slog has a logarithmic divisor a+b=0 with
   q!=0; alpha is regular there generically. Its coefficient vanishes.
4. A lifted loop around a+i*b=0 gives a nonzero variation of alpha.
   Hence its coefficient vanishes, leaving only the constant function.

Thus {1,L,alpha,V,U} is K-linearly independent. The enlarged family
{1,L,alpha,V,U,Slog} is also independent. V=h*(beta+pi/2) differs
from h*beta only by an algebraic term, which does not affect the loops.
Apparent cancellations of U on the *physical* collinear face do not
invalidate the argument: the beta-changing gamma continuation is taken
before delta, on a different logarithmic branch.

## Generation and the two five-atom choices

`verify_atoms.py` gives an exact cyclic-rank witness at (a,b,d)=(3,2,2).
Modulo rational functions and a constant scale, chi has vector

    z=(0,-2pi*v,0,2t,c)

in the independent basis (1,L,alpha,V,U). The columns
1,z,partial_a z,partial_b z,partial_a² z, with covariant derivatives
using the connection, have determinant -4096*pi/21021. Thus 1 and chi
generate this entire differential module over K. Closure gives the upper
bound five, and independence plus this witness gives the lower bound
five. **The smallest rational differential module containing 1 and chi
has rank five.** A sixth separate algebraic radical is not needed.

For the agreed log(sigma) expansion, only the physical combination

    chi_sigma = chi+kappa*v*Slog, kappa=(pi-2)/(3pi)

is needed at order two, modulo polynomial terms. Define

    Ltilde=L-(pi-2)/pi*Slog.

Then chi_sigma has the same five-component vector as chi with L replaced
by Ltilde. DLtilde is rational. The only change in its connection lies
in the constant target row. The determinant just given includes the
constant basis column, so this change does not change its value.
Independence follows from the same double-variation argument and the
d=-(a+b) logarithmic divisor of Ltilde. Therefore **the smallest module
containing 1 and the physical order-two log(sigma) coefficient also has
rank five**, with basis {1,Ltilde,alpha,V,U}.

If one instead demands that chi *and* Slog be independently representable
in the same module, their joint closure has rank six. The full saved
six-atom stencil is useful for this audit and for checking log-convention
conversion. The five-atom physical variant removes this unnecessary
separation. These are different generation questions, not inconsistent
rank claims.

This proves a minimal rank, not a unique list of generators, an optimal
clearing factor, or a result about order-three/four integral atoms. It
also says nothing about independence with *arbitrary function-valued*
Laguerre coefficients, where the reconstruction gauges remain.
