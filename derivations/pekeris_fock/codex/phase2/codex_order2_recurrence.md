# Order-two atoms and their exact coupled Laguerre recurrence

Codex, 2026-10-21. Independent derivation before any `claude2/` comparison.
Inputs are the repository formula for chi and the previous independent
distance-Hamiltonian derivation `codex/codex_derive.py`. Neither `sym/` nor
`julia/` was needed. All operators use precisely COLLAB's P normalization,
with symbolic lambda in the Hamiltonian. The physical Fock atom formulas
use lambda=1; operator identities do not specialize lambda.

## Closed set and its scope

Use a,b,d,q,p,s,c,v,g,h,L,alpha,beta,U as defined in `codex_atoms.md`,
p=2q-d², h=sqrt(p)/d, t=d²/q, f=p/d². The order-two chi module closes on

    one=1, L=log((a+b+d)/sqrt(q)), x=alpha,
    V=h*(beta+pi/2), U=beta*A+alpha*B0-Delta+pi*J.       (1)

In particular

    chi=(s-2(a+b)d/q)/6-v*L/3+(c*U+2t*V)/(6pi)+h0*v.

The radical is included in V; a separate h atom is unnecessary. The
initial six-generator set {1,h,L,x,h*beta,U} reduces because the h and
h*beta contributions always occur in the combination V, including in
every derivative of U.

The *agreed log(sigma) ansatz* additionally needs

    Slog=log((a+b)/sqrt(q)).                             (2)

Indeed psi21=-Z*kappa*v and log R=log sigma-log(2eps)-Slog,
so its order-two log-free angular coefficient is
psi20+Z*kappa*v*(log(2eps)+Slog). The constant log(2eps)
multiplies a polynomial and introduces no atom. Thus the delivered
recurrence uses six atoms, in order (one,L,x,V,U,Slog). The five-atom
subsystem alone is appropriate if retaining log R or auditing chi alone.

**PROVED minimal ranks.** `minimal_rank.md` supplies the actual-function
independence proof by logarithmic monodromy, including a nonzero double
variation of U. The exact cyclic-rank witness in `cyclic_rank.json`
then proves rank five for the smallest rational differential module
containing 1 and chi. Rank six is necessary only if chi and Slog must
be independently representable. For the physical order-two log(sigma)
combination there is a smaller choice:

    Ltilde=L-(pi-2)/pi*Slog,
    {1,Ltilde,x,V,U}.                                  (2a)

Indeed chi+kappa*v*Slog has the same expression in (1) with L replaced
by Ltilde. Only its rational constant derivative changes. This five-atom
module is also minimal; see the generation and independence arguments
in `minimal_rank.md`. Thus no sixth atom is necessary for that fixed
physical combination. No failure of finite closure has been found.

## Rational connection, explicitly

Work on a folded open physical chart a>b>0, |a-b|<d<a+b.
For any coordinate differential D=partial_u,partial_v,partial_w, set

    kh=Dp/(2p)-Dd/d,
    kg=Dp/(2p)+Dd/d-Dq/q,
    z=c+(a+b)d/q, den=s²-v²,
    J0=2Dz/z-Ds/(1+s)-(t²*f*kg-c*Dc)/den,
    J1=-t*(Dc-c*kg)/den,
    Tx=2*(s*Dv-v*Ds)/den,
    TV=2*t*(c*kg-Dc)/den.

Then

    D one=0,
    D L=D(a+b+d)/(a+b+d)-Dq/(2q),
    D x=2*(a*Db-b*Da)/q,
    D V=kh*V+q*Dv/d²,
    D U=pi*J0+Tx*x+TV*V,
    D Slog=D(a+b)/(a+b)-Dq/(2q).                       (3)

All coefficients are rational in distances, hence rational in u,v,w.
The cancellation pi*J1=(pi/2)*TV eliminates the separate h term.
The apparent den and z singularities simplify considerably: the saved
`connection.json` contains the factored perimetric expressions.
For these homogeneous functions one may use a=2v+w,b=2u+w,d=2(u+v)
(common 4eps factors cancel).

**PROVED:** differentiate Delta with Cl2'=-log|2sin(angle/2)|,
obtain dDelta=B0*dalpha+A*dbeta, and hence
d(beta*A+alpha*B0-Delta)=beta*dA+alpha*dB0.
Use g=t*h, Dh=kh*h, Dbeta=Dv/g and
1/(g+c)=(g-c)/(s²-v²). These give (3), including the
radical cancellation. They also give an elementary finite proof of
D-finiteness on chart interiors. Neither differentiating a black-box
numerical chi nor an infinite derivative alphabet is involved.

In matrix notation D_j g_source=sum_target g_target C_j[target,source].
The checked integrability identity is

    partial_i C_j-partial_j C_i+C_i C_j-C_j C_i=0.

Branch continuation, exchange and collision limits require matching the
actual functions, not just solving this flat rational connection with
arbitrary integration constants.

## Atom and log-sector operator blocks

Write P=M_ij partial_i partial_j+b_j partial_j+c0 with symmetric M,
including both mixed entries in the summation. The exact atom blocks are

    P_ab=delta_ab P
          +2 sum_ij M_ij C_i[a,b] partial_j
          +sum_ij M_ij (partial_i C_j+C_i C_j)[a,b]
          +sum_i b_i C_i[a,b].                       (4)

Let T=sum_ij M_ij and sigma=u+v+w. Define

    A_ab=sigma² P_ab,
    B_ab=delta_ab C1 +2sigma sum_ij M_ij C_j[a,b],
    D_ab=delta_ab T,

where C1=2sigma sum_ij M_ij partial_j+sigma sum_i b_i-T
is the original scalar log operator. For every j,

    sum_b [ A_ab F_(j,b)+(j+1) B_ab F_(j+1,b)
                      +(j+1)(j+2) D_ab F_(j+2,b) ]=0. (5)

The saved names are S2P=A, C1=B, C2=D. Product differentiation
proves (4)-(5). They define a sufficient coefficientwise formal system;
unrestricted function-valued F's can have atom gauges, so reconstruction
injectivity must not be silently assumed.

There are 13/12/6 nonzero blocks respectively. All six diagonal blocks
occur. Off-diagonal source->target edges for A are

    L->one, x->one, V->one, U->one, U->x, U->V, Slog->one.

B has the same edges except U->x; D is diagonal. Set weights
one=0; L,x,V,Slog=1; U=2. Atom coupling never increases
this weight, while log powers couple j<-j,j+1,j+2. This finite atom
filtration supplies no highest radial log sector in an infinite tower.

## Clearing factors and finite stencils

Define

    Q2=2u²+2uw+2v²+2vw+w²,
    P2=(u-v)²+2(u+v)w+w².

Multiply equation (5) in each target atom row by

    d_one=uvw*(u+v)*Q2,
    d_V=uvw*(u+v)*P2,
    d_L=d_x=d_U=d_Slog=uvw.                           (6)

These are nonzero on the open perimetric cone. `cleared_operators.json`
contains the resulting polynomial differential operators. The additional
uvw factor is deliberate: it makes every monomial coefficient of
partial^d divisible by u^d1 v^d2 w^d3. This gives a simple finite-stencil
proof and larger counts than a specially optimized Laguerre-tail
cancellation. These counts must not be compared to the old 33-term
recurrence as if the clearing prescription were identical.

For completeness, every JSON entry can be reproduced without any
infinite-tail cutoff. In one variable,

    x^d partial^d L_n = n_under_d sum_(r=0)^d (-1)^r C(d,r)L_(n-r),
    x L_n=-(n+1)L_(n+1)+(2n+1)L_n-nL_(n-1).

For x^p partial^d with p>=d apply the second identity p-d times
to the first. Tensor the three resulting finite maps and sum over
the cleared polynomial monomials. For an output equation at (l,m,n),
the entry at offset delta multiplies A(l+delta1,m+delta2,n+delta3).
Negative source indices are omitted. For a monomial with powers p,
each coordinate offset lies in [-p_i+d_i,p_i]. The saved expanded
coefficients and their exact zero cancellations determine the final
offset lists. No guessed derivative-tail truncation enters.

For the **minimal five-atom physical log(sigma) recurrence**, delete the
Slog row/column, rename L to Ltilde, and replace its constant-target
blocks by

    A_(one,Ltilde)=A_(one,L)-(pi-2)/pi*A_(one,Slog),
    B_(one,Ltilde)=B_(one,L)-(pi-2)/pi*B_(one,Slog).

All other blocks are unchanged, and there is no C2 off-diagonal
replacement. To keep polynomial coefficients in pi, additionally scale
the one target row by pi. `minimal_sigma_corrections.json` contains the
two fully expanded replacement stencils **including this pi row scale**;
multiply the other one-target blocks in the six-atom JSON by pi.
`minimal_sigma_operators.json` and `minimal_sigma_counts.json` give the
corresponding cleared differential operators and counts. This prescription
fully specifies every block of the minimal recurrence without duplicating
the large unchanged stencil tables. It has 11/10/5 nonzero blocks.
`counts_minimal_sigma.md` lists all its counts and shell ranges. Its
replacement block counts are 1327 (S2P) and 841 (C1), unchanged by
E=-eps². Overall shell half-bandwidth is 11 at general E and 10 at
E=-eps² for this clearing, for either the audit or minimal system.

`codex_recurrence.json` uses the original offset-to-polynomial format
inside each `S2P/C1/C2` and `(target_atom,source_atom)` key. E is symbolic;
set E=-eps² for the fixed-energy recurrence. `counts.json` records both
counts, shell ranges, coordinate bounds and differential monomial counts
for every block. A companion Markdown table is generated in `counts.md`.

## Verification and remaining limits

- `order2.py`: exact rational connection flatness, radical cancellation,
  and rational block construction from the independent distance operator.
- `verify_atoms.py`: direct numerical differentiation of the actual
  Clausen expressions, source-chi agreement, and an exact formal rank
  witness. Numerical checks use interior points and are not boundary
  or branch-matching proofs.
- `verify_operator.py`: direct differentiation of random finite
  atom/Laguerre/log combinations against (5), including a fixed-energy
  test. This uses actual polylogarithms, not the connection to differentiate
  the left side.
- `verify_stencils.py`: independent exact monomial differentiation and
  Laguerre moment projection against every saved block stencil, with a
  random rational finite coefficient combination. Pi is specialized to a
  rational value for this exact algebra check; the direct action check
  uses true pi. The stdout files state the actual case totals/errors.

Final recorded checks: 108 exact flatness entries; 12 exact polynomial
homogenizations; 36 direct first derivatives and two source-chi comparisons
(maximum error 1.93e-45 at 45 decimal digits); three full random actions
(maximum relative error 1.93e-40 at 40 digits); all 31 six-atom stencil
blocks and both five-atom replacement blocks pass exact projection tests.
The generation determinant is -4096*pi/21021, separately from the
analytic independence argument in `minimal_rank.md`.

The finite recurrence is established for this order-two atom content.
It is not a closure result for the elliptic/inverse-operator atoms at
orders three/four or a global convergence
theorem. Clearing introduces no new interior solution constraint but
does not establish physical boundary/domain conditions.
