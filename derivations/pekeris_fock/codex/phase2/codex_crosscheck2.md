# Sync 4: independent recurrence cross-check

Codex, 2026-10-21. Inputs are the saved Claude eight-atom JSON and the
previous independent Codex distance Hamiltonian and six-atom connection.
No Claude verifier or recurrence generator is used in the direct check.
All writes for this task are confined to the requested allowlist.

## Explicit change of basis and connections

Here a=r1, b=r2, d=r12, q=a²+b², P=sqrt(2q-d²), s=a+b.
On the real open physical chart, with the branches in the saved formulas,

    one = one_Claude,
    L = la,
    alpha = al,
    V = (beP + (pi/2) P)/d,
    U = Om + pi lb,
    Slog = -lR.

For the minimal physical system replace L by
Ltilde=la+((pi-2)/pi) lR and omit independent Slog. The audit six-atom
system, rather than this physical five-atom combination, is the most
direct comparison with Claude's eight atoms.

To see U explicitly, write c=(a²-b²)/q, v=(q-d²)/q,
g=d P/q, z=c+s*d/q. The two coefficients of beta and alpha in Om are
A=log|(g-c)/(g+c)| and B=log|(2ab/q+v)/(2ab/q-v)|. Thus
Om=beta*A+alpha*B-Delta. The extra logarithm in U is pi*J,
where J=log(z²/((1+2ab/q)*(g+c)))
=log((a-b+d)²/(d P+s(a-b)))=lb. This fixes the additive constant as
well as the derivatives. These identities apply to the physical real
branches; arbitrary complex continuations require compatible branches.

`sync4_connection.py` constructs the 8-by-6 rational transformation T.
For every distance derivative it checks exactly

    dT + Gamma_Claude T - T C_Codex = 0.

All 144 symbolic entries vanish. Twelve actual-function identities at
two points, including both orderings of a,b, also pass (maximum absolute
error 5.35e-50 in that run). The connection check reads Claude's Gamma
formula but verifies its intertwining against the independently derived
Codex connection. The direct stencil check below does not use that Gamma.
See `sync4_connection.out`.

## Direct action of the saved recurrence

`sync4_direct.py` reads all 25,955 saved offsets, evaluates the Laguerre
coefficient polynomials using exact rational arithmetic, and reconstructs
the resulting functions after dividing each target row by its saved D.
It compares against mpmath differentiation of actual polylogarithmic
atoms under the independently derived P, multiplied by sigma².
For each of the eight source atoms it uses a rational combination of
L_0 L_0 L_0, L_1 L_0 L_1 and L_0 L_2 L_0, each multiplied by log(sigma)^p
for p=0,1,2. Both interior points pass: 48 complete actions, maximum
relative discrepancy 2.33572329193e-39 at 45 decimal digits. Parameters
are eps=1/4, E=-eps², Z=2, lambda=1, precisely the saved specialization.
The choice eps=1/4 identifies the distances with 2v+w,2u+w,2(u+v).
This is a strong independent finite-action check, not a symbolic identity
proof of every saved coefficient for arbitrary indices and eps.
See `sync4_direct.out` for every action and every block count.

## Counts, row clearing and the reporting correction

The JSON contains **16 sigma2P, 15 C1, and 8 C2 blocks**, totaling 25,955
offsets, and maximum shell half-bandwidth 10. This corrects COLLAB's
17/12/8 block claim. `sync4_counts.md` lists all 39 blocks individually,
alongside related audit-system counts. Counts in that related column are
the previous generic-Z/lambda fixed-energy counts, not claims of identical
operators after changing basis.

Let Q2=2u²+2uw+2v²+2vw+w² and
P2=(u-v)²+2(u+v)w+w². Claude uses

    D_Om=D_al=D_lR=D_la=D_lb=uvw,
    D_one=uvw Q2/2,
    D_beP=D_P=uvw P2.

Codex uses uvw on its four log/angle rows,
D_one=uvw*(u+v)*Q2 and D_V=uvw*(u+v)*P2.
Thus even the one row has an extra factor 2(u+v) in the Codex clearing.
Multiplication by a coordinate polynomial convolves Laguerre offsets;
counts cannot be compared before removing these factors. V additionally
has a coordinate-dependent 1/d and mixes beP with P, and U mixes Om
with pi lb, so its transformed differential operator also contains
product-rule derivatives of the transformation.

The directly comparable unmodified diagonal log/angle rows all have
sigma2P/C1/C2 counts 571/419/251: Claude Om, al, lR, la, lb versus
Codex U, alpha, Slog, L where appropriate. The alpha<-Om sigma2P
block and alpha<-U both have 461 offsets. The one<-beP and one<-V
counts happen to be 1371/931 for sigma2P/C1; equality of these counts
is not a statement that these differently scaled blocks are equal.
Claude's degree-one P and beP diagonals have 1195/947/683; Codex's
homogeneous-degree-zero V diagonal has 1617/1297/931. The scalar
one diagonal has 1209/967/683 versus 1617/1309/955 after the larger
Codex row clearing. These differences are consistent with the exact
connection transformation and direct Hamiltonian tests.

**Sync 4 position:** co-sign the saved eight-atom recurrence as supported
by these independent checks, with the block-count correction above.
Finite closure and compatible connections are established; neither the
count agreement nor these finite tests establish a global expansion or
an infinite-dimensional minimality theorem.
