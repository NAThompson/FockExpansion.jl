# TASK_codex_10: cyclic rows and globality

2026-10-22. All jobs are serial. No Claude-owned source is edited. This is
an audit of the stated representations, with the distinction between
represented cyclic ranks and the smallest module of actual functions kept
explicit. Globality and a reduced coefficient field are in goal2_global.md.

**Result:** the exact two-point stabilized quotient candidate is **45**, with
an explicit derivative selection. Under the same quotient, all four panels
give **11/12/13** for N0/N2/joint; the classical piece gives **12**. The
sufficient radical degree bound is proved to be **32**. **Strict smallest
actual-function rank and generic closure remain unproved**, and the later
Claude recurrence cross-check awaits a recurrence file.

## Method and checks

For X satisfying d_i X=C_i X and f=lX, the rows are built with
r_(alpha+e_i)=d_i r_alpha+r_alpha C_i. Taylor coefficients of l **and** C
are retained; merely multiplying frozen point matrices is not used.
The Cauchy calculation uses branch-continuous quartic roots on each complex
polydisc. Principal square roots across the negative discriminant cut would
silently swap the root labels and invalidate a Fourier/Taylor calculation.

The Julia panel port matches 312 deterministic connection actions and four
moment outputs against the original mpmath implementation to <4.1e-16.
The classical port matches 72 independently differentiated atom components
and two physical classical outputs to <6.2e-16. The small matrices and one
active numerical job keep storage modest. BLAS uses one thread.

## Classical piece

For Z=2, the displayed 18-generator classical presentation gives rank 12
at (a,p)=(-.4,-1.3),(-.7,-.9), using orders through eight, radius .15 and
40x40 Cauchy samples. The rank history is 1,3,6,10,12,12,12,12,12.
The thresholds 1e-7,1e-9,1e-11 agree. This uses the nonalgebraic classical
piece C; its algebraic remainder is accounted for separately by 1.

An independent truncated-power-series calculation over F_998244353 gives
exactly the same history at (t,q)=(2,3),(5,7), with
 t=exp(i*a/2), q=exp(i*p/2), pi specialized to 101 and i,sqrt(2) represented
by verified roots in the finite field. The script is goal2_classical_exact.py;
selected derivative rows and parameters are saved in its JSON. An explicit
12-row candidate is

    {d_t^i d_q^j C: i+j<=3} union {d_q^4 C, d_t d_q^3 C}.

The nonzero finite-field minors are exact witnesses of rank at least 12 in
the formal presentation. The entries lift to rational functions over
Q(i,sqrt(2),pi); pi is transcendental, so its specialization is used as a
nonvanishing witness, not as an approximation to physical pi. Stabilization
at finite orders and points does not prove generic closure of these twelve
rows or independence of actual classical functions. Thus 12 is a strong
candidate, not yet a proved smallest actual-function rank.

## Elliptic panel pieces

The high-precision Cauchy run uses 160-bit arithmetic, 40x40 samples,
radius .06, and all mixed derivative orders through eight. At BOTH
(b,y)=(.3,-.6),(-.27,-.7), the ranks are:

| maximum order | N0 | N2 | joint |
|---|---:|---:|---:|
| 0 | 1 | 1 | 2 |
| 1 | 3 | 3 | 6 |
| 2 | 6 | 6 | 11 |
| 3 | 10 | 10 | 14 |
| 4 | 12 | 14 | 16 |
| 5 | 14 | 16 | 18 |
| 6 | 15 | 18 | 19 |
| 7 | 15 | 18 | 19 |
| 8 | 15 | 18 | 19 |

All three singular-value thresholds agree. These are ranks in the TWO
39-state padded families, before actual-function endpoint/zero-jet
relations are quotiented. They are not minimal actual-function ranks.
Each of the four physical panel maps is generically of full angular rank,
so this test describes the generic unweighted panel family on either
endpoint chamber. It does not add the four ranks to infer a combined rank.

Double-precision Cauchy runs, including a larger-radius/density check,
show threshold-dependent spurious high-order directions. Their logs are
retained as diagnostics and superseded by goal2_big_panel.out.

## Cross-check of Claude's current Goal 2 files

At the current snapshot, claude5/goal2_edges.jl/.out samples the original
psi30 atom, its gradient and Hessian along four approach paths. The saved
numbers are bounded over t=1e-2,...,1e-5 on those paths. This is finite
numerical evidence for those values and derivatives only. It is not a
uniform boundary estimate, a test of every derivative atom in a prospective
basis, or a form-domain proof. No goal2 recurrence/connection coefficients
are supplied by that file, so task item 3's independent recurrence rederivation
and coefficient spot checks remain pending its availability. No Claude
source is changed or executed by this review.

The 160-bit label refers to connection/Taylor/row arithmetic; the final
row-normalized matrices are converted to ComplexF64 for SVD. The three
threshold checks test that last rank decision. They are not rigorous
interval-arithmetic certificates for vanishing minors or infinite-order
closure. The exact finite-field calculation avoids that numerical rank
decision for the classical witness only.

## Combined physical output: exact stabilized candidate 45

The joint output contains 1, chi and the classical-plus-four-panel psi30
combination, with exactly the c0_j,c2_j coefficients of structure.md, Z=2.
Algebraic pieces of chi and psi30 may be omitted because 1 is included.
Consequently fixed E,a21 polynomial terms do not affect this generation
question. Both angular derivations are used; homogeneity supplies the radial
one. The four panel Jacobians are included by the chain rule.

The raw realization has 18 classical + 4*78 padded panel + 2 old (L,U)
states, totaling 332. A 332x240 actual-function substitution:

- identifies the A/D zero jets for every integral family;
- sets the positive F0 jets to zero and identifies F0 with its algebraic value;
- replaces endpoint E jets by E0 times 1, its endpoint angle, and one half
  its square, using the existing 1,a,p,a²,p² classical atoms;
- retains the old L,U atoms and all other integral states.

No new functional identities between different panel integrals are assumed.
The quotient is implemented explicitly in both Julia and Python.

The 160-bit Cauchy experiment at (a,p)=(-.4,-1.3),(-.6,-1.1) reaches raw rank
50 and quotient candidate rank 45 at order eight. The second point shows a
one-direction threshold ambiguity at the coarsest 1e-7 threshold, so these
floating-point data alone are not the final rank evidence.

The exact finite-field Taylor implementation (goal2_combined_exact.py) uses
split algebraic charts, 15 Taylor orders, and differential rows through
order 14. It retains only adjacent derivative orders. At BOTH (t,q)=(3,4)
and (4,7), its quotient rank history is exactly

    3,7,11,15,20,26,33,41,45,45,45,45,45,45,45.

The full 332-state connection has EXACTLY ZERO constant curvature at both
points, before the quotient, independently of the chosen output rows.
This is a useful port/chain-rule check. The finite field is F_998244353,
with pi specialized to 101; radicals are computed by exact algebraic Taylor
series, never by finite differences. Non-split points are skipped explicitly.
Both split points select the SAME 45 derivative rows, saved in
`goal2_combined_exact.json`.

An explicit pruning candidate is, with f=R³ psi30,

    {1, chi, d_t chi, d_q chi, d_q² chi}
    union {d_t^i d_q^j f : i+j<=7}
    union {d_q^8 f, d_t d_q^7 f, d_t² d_q^6 f, d_t³ d_q^5 f}.

There are 5+36+4=45 functions. Here

    a=alpha-pi/2, p=p0-pi,
    t=exp(i*a/2)=(1-i)(r1+i*r2)/(sqrt(2)*R),
    q=exp(i*p/2)=(c-i*P)/(sqrt(2)*R),

and angular differentiation holds R fixed. These are algebraic nonzero
angular coordinates; their derivations are algebraic combinations of
physical distance derivations. Replacing the computed output rows by actual
chi and f only adds algebraic multiples of 1 (and a nonzero constant scale
for chi), so the displayed candidate has the same selection meaning.

**Scope of the number 45:** it is an exact two-point, finite-jet stabilized
rank in the stated quotient representation, and an explicit candidate for
pruning that representation. It is NOT yet a proof of generic rank <=45,
closure of this basis as rational/algebraic functions on a whole chart, or
independence of these 45 actual functions. Remaining algebraic relations
between actual panel integrals could reduce the actual-function rank.
Neither a smallest actual-function module nor its rational rank is certified
by these experiments. A generic symbolic row-closure certificate and actual
function/period independence analysis are the next requirements; the
proved combined rational closure bound remains 6912, not 45*32.

## Per-panel quotient calculation (supersedes padded panel counts for pruning)

Using the SAME actual-function substitution as the combined calculation,
`goal2_pieces_exact.py` tests each pulled-back panel separately at both
(t,q)=(3,4),(4,7). All FOUR panels give the following exact histories:

| maximum order | N0 | N2 | joint |
|---|---:|---:|---:|
| 0 | 1 | 1 | 2 |
| 1 | 3 | 3 | 6 |
| 2 | 6 | 6 | 10 |
| 3 | 10 | 10 | 13 |
| 4 through 8 | 11 | 12 | 13 |

The same calculation gives classical rank 12 through orders 4–8. In these
(t,q) coordinates the per-panel N0 candidate is its derivatives of total
order <=3 plus d_q^4 N0 (11 rows). The N2 candidate is its derivatives of
total order <=3 plus d_q^4 N2,d_t d_q^3 N2 (12 rows). The explicit joint
13-row selections for each panel, including the small panel-dependent
choice of mixed rows, are in goal2_pieces_exact.json. These are represented
candidates subject to the same generic-closure/actual-independence caveat
as the combined 45; the padded 15/18/19 counts are not the final pruning
counts.

## Final verification and remaining work

The actual-function quotient is independently checked against two sets of
all 332 states evaluated by the original integral evaluator, with maximum
error <1.2e-16 after Float64 export. Both angular old-L/U pullbacks agree
with independently differentiated actual functions to <8.3e-17. Classical
and panel connection checks, exact radical checks, all two-point histories,
and the common 45-row selection pass goal2_verify.py. Its output is saved
as goal2_verify.out; goal2_input_hashes.json identifies the read-only
connection/formula inputs. Commands are listed in goal2_commands.txt.

**Task status:** the requested per-piece/combined stabilized experiments,
explicit pruning candidate, sufficient radical reduction and globality
assessment are supplied. The strict request for a *proved smallest*
actual-function module and an *established minimal* generator set remains
unresolved: numerical/exact-specialized stabilization does not substitute
for a generic symbolic closure certificate and independence proof. The
rational minimal rank is also uncomputed. The conditional later recurrence
cross-check is pending; only Claude's edge probe exists at the final file
inventory. No theorem of global recurrence convergence or atom form-domain
admissibility is asserted. All computations were serial; only the requested
Codex artifacts and append-only shared logs were written.
