# TASK_codex_14 — classical cross-check, panel budget, backward theorem

Codex, 2026-10-22. Read-only inputs: Claude's phase-7 files. New outputs use
`codex/phase7/x14_*`; shared logs are append-only. No agents, shell solves,
energy calculations or dense reconstruction were used.

## 1. Exact classical cross-check

`x14_check.py` reconstructs the invariant module `{1,a,a²,L(alpha/2)}`
from the distance Hamiltonian independently of Claude's Gamma/pipeline.
The last atom is g9 in the agreed order, and is a genuine L atom.
At eps=1, E=-1, Z=2, lambda=1, with pi exact, **all 21 blocks and 81
stored/expected differential coefficients agree exactly**, after multiplying
the independently derived operators by Claude's target-row D.
This includes sigma²P, C1 and C2, all outgoing rows from these four sources,
and checks for missing/unexpected terms. It does not check the other 17
source sectors or establish a general-eps identity.

The construction uses da=(-2r2 dr1+2r1 dr2)/Q and
`dL(alpha/2)=(a+pi/2) r1/(4r2) da`. For a coefficient F and atom G,
P(FG) has derivative cross term 2 M_ij (d_i G)d_j F and zeroth term
M_ij d_i d_j G+B_i d_i G, in addition to G PF. The C1 cross term is
2 sigma sum_i M_ij d_j G. This supplies the independent module calculation,
without importing Claude's symbolic code or executing its scripts.
Four further exact identities compare these operators to the previously
saved `angle_recurrence.json`, with its different uvw Q clearing.
The first checker run exposed a missing half-angle chain factor in this
new checker; correcting that derivative gave the exact comparisons above.

No `claude7/classical_stencils_*.json` files were present at the recorded
inventory. Thus no Claude stencil coefficient comparison is claimed.
The checker will perform small exact Laguerre-action comparisons for relevant
blocks if shards become available. `x14_check.json` records input hashes,
the actual shard inventory and check counts. Frozen prior independent
563/415/251 counts are not asserted to equal counts with different clearing.

CO-SIGN the saved classical cleared operators on this tested module and
specialization. The claimed complete 87/76/21 block counts are independently
counted in the verification audit; agreement of all their coefficients is
outside this small check. Real distance connections with an extra P sector
are consistent with the independent module; the full 21-sector reachability
is not newly rederived here.

## 2. Elliptic recurrence size: estimates, not an expanded recurrence

Use the fixed-endpoint-1 connection throughout, with global root branches.
For each panel the 55-state presentation is sufficient, not independent.
Restriction of scalars gives the following **safe allocated sector counts**:

| panels | coefficient field | degree bound | allocated states per panel |
|---|---|---:|---:|
| 1,2 | K(A,z1) (transport z2) | 4 | 220 |
| 3,4 | K(P,A,z3) (transport z4) | 8 | 440 |

Starting with all 55 primitives, these are upper bounds on reachable rational
sectors, not measured ranks. Differentiation of a quadratic tower stays in
its binary-monomial span, so repeated Hamiltonian/log products cannot exceed
them. Combined panel allocation is 1320 states; adding the 21 classical
sectors gives 1341 before duplicate removal. Shared constants/endpoints,
conjugacy and integral identities can reduce these numbers. Neither 13 nor
55 should be multiplied by a field degree and called a proved minimal rank.
Adding the sufficient W^18 weights changes the connection by a rational
logarithmic derivative; it does not require more radical sectors. The binary
scalar multipliers do not invalidate this domain weight: on the normalized
cone their values are O(beta^-1) and first derivatives O(beta^-6), using
the phase-7 root bounds (A,P and their products are bounded). Multiplying
the established G=O(beta^-5), dG=O(beta^-16) estimates by such a multiplier
leaves W^18 G times that multiplier bounded with bounded angular gradient:
the worst exponents are 18-16-1=1, 18-5-6=7 and 17-5-1=11.
Their radial degrees are nonnegative (P degree 1, A degree 2, z degree 0).
Thus the same sufficient weight covers these allocated physical sectors;
this is not an unweighted domain claim.

`x14_size.py/.json` measure common cover row denominators directly from
both exact h,z matrices. The sufficient first-order row degrees are at most **12**; including
neighbor rows, sufficient second-order cover row degrees are at most **30**.
The 13 distinct denominator factors are exported. Complex denominator
norms can overclear, so these are sufficient rather than minimal degrees.
These measurements are **cover degrees**. They
must not be quoted as physical-distance clearing degrees. A second-order
row needs denominators from its neighboring first-order rows as well as
its own differentiated denominator. Physical pullback further introduces
Jacobian/Hessian denominators, radical multiplication matrices and the
pole/discriminant factors. A sufficient rational distance clearing can be
constructed by taking norms of these denominators through the tower,
then multiplying target rows by uvw to cancel Laguerre derivative tails.
Clear after restriction, not by pretending an algebraic denominator is a
polynomial in u,v,w. Full norms and optimized distance row degrees have
not been expanded in this task.

Planning estimate: cover denominators of degree tens become physical
clearings of degree tens to hundreds after rational substitution, second
products and degree-4/8 norms; unoptimized norm products can be worse.
This is a budget range, not a certified maximum. If the final polynomial
operator coefficients have total degree d and differential order <=2,
a finite stencil has shell half-bandwidth at most d once tails cancel,
and an enclosing offset cube of at most (2d+1)^3 positions per block.
For d=30,60,120 these cubes have 226981,1771561,13997521 positions.
Actual sparse support can be much smaller. Even the connection has a sizeable
allocation: the two cover matrices together store 282+293=575 entries.
A k-dimensional scalar restriction can replace one scalar multiplication
by up to k² entries; including scalar-basis derivatives gives conservative
union budgets (575+55)k² = 10080 or 40320 first-order entries per panel.
Second-order path products can fill more blocks (absolute block ceilings
220²=48400 or 440²=193600). These are ceilings, not observed sparsity;
multiplying them by stencil cubes is intentionally not an actual term count.
Large entry degrees, common-row
degrees and number of blocks are distinct costs. The 216-degree weight
itself need not be expanded: conjugating by W^18 adds 18 dlog W and its
second-derivative products, hence fixed low powers of W in denominators;
it does not force expansion of a degree-216 polynomial into every row.

A full eager expansion of every rationalized 55-state panel is therefore
poorly suited to small serial jobs. A lazy sparse operator connection is
feasible; a selected physical output row/block with common factors retained
is the sensible next expansion. The exact 55-state connections already
avoid the >350-degree combined cyclic-vector bottleneck.

The 13-atom per-panel cyclic representation offers 52 panel states if its
rational descent is proved. Its sampled partial-connection entries have
line degrees about 48, but that is neither a bivariate bound nor a common-row
denominator degree. LCM growth and second derivatives may offset its state
saving. Prefer it only after obtaining exact bivariate certificates and
common-row degree measurements. Pairing pole integrals under +/-z and h/z
can remove redundant root characters and is a promising intermediate
presentation: derive exact paired identities and their parameter derivatives
before any rank/minimality claim. No such reduction is asserted here.

## 3. Step 6: extension and precise hypotheses

YES for the **physical weak-pencil theorem**, under the same isolated simple
positive ground-state assumptions, when the trial spaces retain the plain
singlet polynomial core (or have a separately proved dense union).
Take nested finite spaces U_N of exp(-eps s) times polynomials, classical
atoms times polynomials, and W^18 times fixed-endpoint panel atoms times
polynomials. Use symmetry projections of complex atoms as needed. The
phase-7 uniform domain arguments supply admissibility of every generator;
plain polynomials supply density. Pointwise normal-path fits alone do not
supply that domain premise for unweighted higher derivative atoms.

Choose a form-bounded ell with ell(psi0) nonzero, for example the positive
physical overlap with exp(-eps s). Choose one fixed admissible b0 in U_N
with ell(b0)=1. Center the other functions by subtracting ell(phi)b0 and
quotient exact physical dependencies. Solve all weak rows in
W_N=U_N intersect ker ell at the exact E0 with ell(u_N)=1 and zero exterior
coordinates. The spectral gap implies b=q-E0 is coercive on ker ell with
constant beta>0. Galerkin orthogonality gives

    ||u_N-psi0/ell(psi0)||_form
      <= beta^(-1/2) inf_{v in U_N, ell(v)=1}
                            ||v-psi0/ell(psi0)||_form -> 0.

The deleted head row equals b(u_N,u_N), decreases to zero, and nonreal Schur
convergence and Ritz convergence extend by the same argument. Small form-norm
exterior tails give stable corrections. Weights do not impair the argument
because the retained plain core supplies approximation independently of them.

**Admissibility plus a bounded positive normalization alone do not suffice.**
For example, repeating a fixed finite admissible space with a positive head
need not approximate the ground state. Nested density, correct symmetry,
the spectral gap/simple ground state, independent physical coordinates and
all nonhead physical weak equations are needed. Positivity of ell on the
positive ground state is enough; positivity on every generator is unnecessary.
If only weighted panel/classical functions are kept, density must be proved;
it is not implied by their form-domain membership. The existence of a dense
plain core cannot be inferred merely from a numerical Fock-seed replacement.

This proves no equality with the formal denominator-cleared Laguerre sector
recurrence. Clearing can hide boundary/interface conditions and raw atom
coordinates can be redundant. That identification, admissible infinite
coefficient reconstruction and a raw coefficient minimality theorem remain
separate. W^18 also removes the panel contribution to the low corner jet;
this extension does not establish an unweighted exact order-three Fock seed.

## 4. Co-sign / dissent

CO-SIGN the classical recurrence on the exact tested specialization/module,
and the panel split as a route to manageable structured connections.
CO-SIGN sampled energy/per-size observations only as Claude's numerical
observations; no new energy calculation was requested or run.
DISSENT from promoting normal-path regularity fits to uniform unweighted
admissibility of all panel/cyclic atoms, or sampled rank/branch invariance
to a proved minimal rational basis. The fixed-endpoint globality and
sufficient weighted domain proof remain the applicable exact statements.
The backward theorem extends with the hypotheses above; it does not by
itself certify the sparse formal recurrence as a physical selector.
