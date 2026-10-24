# TASK_codex_12: independent Goal-2 cross-check

The independent branch and constant-field tests pass. Actual minimal
numerator/denominator degrees for the psi30 boundary relations remain
unresolved: the saved samples reject the tested low-degree fits. No
bivariate reconstruction or generic symbolic closure theorem is claimed.

## Independent experiment

Prime **1004535809**, pi specialization **103**, points **(10,17)** and
**(13,22)**; these differ from Claude's prime/default specialization and
announced point. The phase-five connection is reused, with independently
written row elimination, target selection, sign/embedding controls and
rational interpolation. This is an independent cross-check of that
connection, not a second derivation of the connection itself.

The saved 45-row basis is reduced in the actual 240-state endpoint/zero-jet
quotient. All missing psi30 first-derivative boundary rows are tested:

* (source,a,b)=(3,a,8-a), a=4,...,8;
* (source,a,b)=(3,a,9-a), a=0,...,4.

Each of the ten relations has 45 coefficients, and each solve is checked
against all 240 represented coordinates. Bivariate Taylor order 9 includes
derivatives of coefficient rows and the connection. It is not a test using
powers of constant connection matrices.

At both points, **all eight** independent sign choices of sqrt(D),z1,z3
produce identical 450-entry arrays. This is strong specialization evidence
for descent through those radicals; it does not prove generic rationality
or independence of the actual functions.

## Constants

Vary pi through 103,107,109,113,127. Four specializations train each
rational fit, and the fifth is held out. At both points:

| Entries | Minimal training degree pair | Held-out check |
|---:|---:|---|
| 10 | (2,1) | pass |
| 40 | (1,1) | pass |
| 400 | (0,0) | pass |

Every nonconstant denominator is exactly pi after monic normalization.
The observed Laurent support is therefore {-1,0,1}. The 50 pi-dependent
entries are in basis columns **0,1,3,5,7**, the five old-module basis
rows. All i-conjugation tests are invariant. Conjugating sqrt2 negates
exactly those 50 entries and fixes the other 400. The conservative
reconstruction field Q(i,sqrt2)(pi)(t,q) is well posed; the experiments
suggest Q(sqrt2)(pi)(t,q). Divide the odd entries by the chosen sqrt2
embedding and reconstruct their pi Laurent components separately.
These are finite-field embedding tests, not automorphisms of F_p.

## Independent line-degree tests

Line **t=10+u, q=17+2u**, 261 distinct split regular points. No failed
point was replaced with a zero. All 450 coefficients are saved at every
point. Polynomial extended Euclid finds minimum-total-degree regular
training interpolants; a separate five-point set validates candidates.

| Training points | Entries tested | Balanced training-fit pair | Pass holdouts |
|---:|---:|---:|---:|
| 32 | all 450 | (16,15) | 0 |
| 128 | all 450 | (64,63) | 0 |
| 256 | 40: columns 0,1,5,44 in each boundary row | (128,127) | 0 |

These displayed pairs are **interpolation artifacts**, not true degrees.
The exact training searches give total-degree lower bounds 127 for all
450 entries and 255 for the selected 40. Actual minimal numerator and
denominator degrees have not been found, nor can a degree upper bound be
inferred from the data. This is the remaining incomplete numerical part
of step 1. The full data support continuing the probe without recomputing
existing points. A generic symbolic result would additionally need degree
bounds and reconstruction verification, not just held-out samples.

The independent results, code and this limitation were saved in
`goal2x_step1.md` and checkpoint hashes **before any claude6 inspection**.
The independent line is parallel to Claude's later line but has a different
intercept and characteristic; neither passes through the other's base point.

## Comparison with Claude, read only after the independent checkpoint

`guess_probe.py` chooses the first six rows outside the basis at total
orders 8/9. In its actual iteration order they are
(1,0,8),(2,0,8),(1,1,7),(2,1,7),(1,2,6),(2,2,6).
There is **no source-3/psi30 row**. Thus its original reported sign test
does not test the claimed psi30 boundary relations. Likewise the original
`guess_line.json` targets are (1,0,8),(2,0,8),(1,1,7): its zero/(0,0)
results concern trivial constant derivatives, while its chi coefficient-0
fit is unresolved. Those outputs cannot establish psi30 degrees.

The later **`guess_line2.py` corrects the target selection**: four chi rows
plus the ten psi30 boundary rows above, in a different order. No saved
bivariate reconstructed coefficient formulas were available on inspection;
`ratrec.py` is a reconstruction routine, not a reconstruction certificate.
The current corrected line sample is too short to settle the independent
high-degree obstruction. There is no contradictory validated psi30 degree
claim in the inspected outputs.

Cross-checked the corrected exported coefficients at two fresh points
**(177,352)** and **(187,372)** over Claude's prime 998244353, pi=101:
**630/630 match at each point**. Used this driver's unaccelerated Python
jets, Taylor order 9, independent elimination and aligned i/sqrt2 embeddings.
An initial comparison exposed the opposite sqrt2 embeddings: precisely
50 psi30 coefficients changed sign, as predicted by the independent tests;
after alignment the discrepancies vanish. Frozen input copies and SHA256
hashes are in `goal2x_claude_*` and `goal2x_comparison.json`.
These check saved modular solves, not unavailable reconstructed formulas.

## Verification and scope

`goal2x_verify.py` freshly checks saved branch/constant evidence, repeats
three 256-point degree searches, compares all 450 accelerated and pure-Python
entries, and tests full radical conjugation at a third independent point
**(25,34)**. The C accelerator performs only modular convolution/inversion,
reducing every eight products to avoid uint64 overflow. No agents, worker
pools or multithreaded linear algebra were used. Each sample is a small
serial computation, about 1.7 seconds after optimization.

All task writes are `codex/phase6/goal2x_*`, plus an append-only COLLAB entry.
Claude-owned source, outputs, task files and the repository were not modified.
Actual-function minimality, generic closure and actual minimal line degrees
remain open. The task's experiments and comparison are delivered with these
explicit limitations; it is not an exact reconstruction handoff.
