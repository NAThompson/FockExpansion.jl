# TASK_codex_16 — 2026-10-23

All jobs were small and serial; no agents, energy solves, or stencil expansions.
Only x16_* deliverables and the authorized append-only shared logs were written.

## Four-panel exports

Use r1=x, r2=y, r12=c, Q=x²+y², e=2Q-c². Atom definitions and factorial
jet normalization are exactly those in `report.md`, section 1. Every panel
uses its own copy of the same integral family with its stated H,k pullback.
`x16_summary.json` gives all four maps, reachable names and removed names.
Panel 2 exchanges x,y in panel 1; panel 4 exchanges them in panel 3.

* `x16_raw55_connection.txt`: all four full 55-state matrices, strictly over
  Q(x,y,c), format `panel|axis|row_name|col_name|numerator|denominator`.
* `x16_connection.txt`: dependency-reduced, output-reachable matrices in the
  same format. Entries are rational over Q(pi)(x,y,c); pi appears only from
  eliminating constant-angle endpoint duplicates. Both files use dX=M X.
* `x16_output.txt`: `panel|output|atom_name|coefficient`, exact contribution
  to **R³ psi30**, helium Z=2. All coefficients are in Q(pi)(x,y,c).
  **No overall R, sqrt(2), or other algebraic factor remains.**
* `x16_clearing.txt`: `panel|clearing|row_name|polynomial`, joint sharp
  polynomial denominator across all three reduced derivative matrices.
  Rows with no derivative entries have clearing 1. Sharpness is up to a
  nonzero coefficient-field unit and refers to this presentation.

| panel | reachable states | d/dr1, d/dr2, d/dr12 nonzeros | max clearing degree |
|---|---:|---|---:|
| 1 | 23 | 84,84,69 | 11 |
| 2 | 23 | 84,84,69 | 11 |
| 3 | 23 | 76,84,84 | 10 |
| 4 | 23 | 84,76,84 | 10 |

These are graph-reachable sufficient presentations, not certified minimal
actual-function bases. There are no duplicated zero jets. Further functional
relations, beyond those explicitly proved below, have not been excluded.

### Exact reductions

The density endpoint E_A0 equals x/c in panel 1, y/c in panel 2, c/(2x) in
panel 3 and c/(2y) in panel 4. E_A1=(pi/4)E_A0 and E_A2=(pi²/32)E_A0.
Eliminate these three into F0=1. Their derivative rows are discarded;
substitution in retained rows needs no derivative-of-basis correction,
since those retained atoms themselves are unchanged.

For each reachable odd-quartic pair Q1_Fn,Q3_Fn, n=0 or 1, integration of
(f*angle^n/n!)' and (s²*f*angle^n/n!)' gives

    k Q1 + 2 Q3 = a,
    2H Q1 + k Q3 = b,
    Q1=(k*a-2*b)/(k²-4H), Q3=(k*b-2H*a)/(k²-4H).

At n=0, a=1-E_A0, b=E_A0. At A1,
a=-E_A1+U0_A0, b=E_A1-J_A0+U0_A0. At B1,
a=-E_B1+H V0_A0, b=E_B1-H J_A0+H² V0_A0.
Thus six further dependent states disappear. Retain only the closure of the
output support under the three derivative matrices: 23 states each.
Factors k²-4H or pole-coincidence factors can vanish in the physical
interior. At these loci use analytic limits/cleared identities; the
underlying real integral atoms remain regular. This is not a new singularity
of the physical functions.

### Output normalization

Let F0out,F2out be the two rows in `outputs.json`, after substituting H,k.
These are normalized moments, not the constant atom F0. Put b=(H-1)/(H+1)
and Y2=(1+k/(1+H))/2. For panels 1/2 set a=x/y respectively and

    N_j = R/(sqrt(2)*a) * Fjout,
    c0 = -sqrt(2)*Z*(9Z*b²-18Z-8b²+8)/576,
    c2 = -sqrt(2)*Z*(1-b)*(90Z*Y2-45Z-18b-32Y2+16)/288.

Multiplying R³(c0 N0+c2 N2) cancels sqrt(2) and R, giving the exported
coefficients. For panels 3/4 N_j=R/c*Fjout, and the rational output is
Z²*c*e/18*F0out - 2Z²*e*(a²-d²)/(9c)*F2out,
where (a,d)=(x,y) or (y,x). All exports specialize Z=2.

## Independent Claude comparison

`x16_compare.py` reads the saved `claude8/panel1_connection.json` directly,
including its real/imaginary parts and total-degree monomial ordering. It
constructs actual symmetric-pole integral states, applies those saved
matrices, and computes physical N0,N2 and both h,Y2 derivatives. Moment-row
parameter derivatives are included. The independent rational-family side
uses H=h², k=(h²+1)(2Y2-1), normalized density and B=h atan(s/h).

At fresh points (h,Y2)=(1.31,.37),(.79,.44), all 12 value/derivative
comparisons pass, worst absolute difference 1.079e-45 at 45 digits.
This is a fresh-point actual-function cross-check, not a proof that every
entry in Claude's reconstructed bivariate connection is an identity.

`x16_verify.py` additionally checks 5592 exact rational raw-entry pullbacks
at two distance triples per panel, all reduced matrix actions and physical
panel outputs against actual integrals, and every reduced joint row LCM.
Worst actual discrepancy is 2.067e-40 at 40 digits. Logs are
`x16_export.out`, `x16_compare.out`, `x16_verify.out`.

## Admissibility and closure: progress and precise limits

The existing delta^4 sufficient weight remains valid, with

    delta=[4x²y²-(Q-c²)²]/Q²,  w=delta^4.

It applies to all 55 atoms for all four maps. For panels 1/2 the normalized
density is x/sqrt(y²s⁴+(c²-Q)s²+x²) (exchange for panel 2).
On Q=1 the coefficient matrix has trace 1 and determinant Delta/4;
for panels 3/4 its trace is 2 and determinant Delta. The same uniform lower
bounds and B-angle estimates used in `report.md` therefore apply with only
constant changes. The reduced atoms are a subset of the original atoms,
so their admissibility follows without relying on singular reduction
coefficients.

**Necessary lower bound for power weights.** Along a generic collinear
stratum the earlier Q0_A0 example has X~C/delta. In a transverse Cartesian
radius rho, delta~C' rho² and the volume element is rho d rho. For
w=delta^a, 0<a<1, the weighted gradient has leading nonzero term proportional
to rho^(2a-3); its squared integral is integral rho^(4a-5)d rho, divergent.
The unweighted atom also fails. Thus a universal power in [0,1) cannot
work. At a=1 the leading derivative cancels, so this example does not rule
it out. The smallest sufficient integer power lies in {1,2,3,4}; its exact
value remains OPEN. No minimal-weight claim is made.

**A form-closure theorem (not a finite-span/jet theorem).** For every eps>0,

    closure_H1 { exp[-eps(x+y)] * w * P(x,y,c) }
      = the rotation-invariant H1(R6) subspace.

One may restrict P to invariant Cartesian polynomials, which are polynomials
in x²,y²,(Q-c²)/2, hence are allowed distance polynomials. Consequently the
polynomial-times-weighted-atom span already has this closure through the
constant atom alone. Exchange symmetry is imposed by averaging.

Proof: the zero set of delta is the collinear set, including collisions.
Away from the origin its smooth largest stratum has codimension two in R6;
its lower-dimensional strata have codimension at least three. It has local
H1 capacity zero. Explicitly, a transverse logarithmic cutoff from rho=eta²
to eta has Dirichlet integral O(1/|log eta|); lower strata admit the usual
higher-codimension cutoffs. Exhaustion and partitions of unity give cutoffs
for compact sets. Thus invariant smooth compactly supported functions
supported away from that zero set and the origin are H1-dense. For such a
function u, g=u/w is smooth and compactly supported in the regular region.
Multiplication by w is bounded H1->H1: w is bounded and |grad w|<=C/R,
so Hardy's inequality in dimension six bounds the extra term.

For completeness, exp[-eps(x+y)] times Cartesian polynomials is an H1
core. If an H^-1 functional annihilates them, its pairing with
exp[-eps(x+y)] exp(i z dot X) is holomorphic for sufficiently small imaginary
z. Every Taylor coefficient at z=0 vanishes. Analytic continuation gives
zero for all real z, hence the distribution obtained by multiplication by
the positive exponential is zero by Fourier uniqueness. The original
functional is zero (locally the exponential and its reciprocal are H1
multipliers). Hahn-Banach gives density. Rotation averaging preserves the
exponential and yields invariant polynomial approximants. Approximate g
by this core and multiply by w, proving the theorem. Coulomb form and H1
norms are equivalent after a fixed positive shift, by Hardy/form bounds.

**Meaning for exact psi30.** Every exponentially damped panel output that
is H1 belongs to this closure; an admissible full psi30 output belongs as
well. This removes any extra approximation obstruction caused solely by
using the weight. It does not prove H1 regularity of each individual
unweighted physical panel output: the special cancellations in those rows
still need uniform boundary estimates. That remaining panel-regularity
question is explicit, not inferred from finite-difference paths. The theorem
also does not assert a finite polynomial representation, a convergence
rate, or preservation of the exact R³ angular jet by approximants. A
homogeneous coefficient/jet topology is stronger here; dividing output rows
by w remains necessary for their immediate exact rational representation.

**CO-SIGN:** the same rational family covers all four panels; Claude's
independent saved connection passes the specified fresh physical-moment
checks. **DISSENT:** universal unweighted primitive admissibility and any
inference from form closure to finite exact angular-jet inclusion. Minimal
weight and uniform individual-output H1 estimates remain open research
items. The requested continuation has explicit new reductions, lower bound,
and a form-closure theorem, not a claim of complete minimal-basis theory.
