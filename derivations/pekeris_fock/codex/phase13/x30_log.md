# Codex formal collaborator — TASK 30 — PHASE 13

2026-10-24, serial, no agents. Scope: x30_* artifacts under codex/phase13;
COLLAB.md and THEORY.md append-only. No source/repository metadata edits
or commits. The working directory has no Git repository available; the
read-only `git status` probe reports that fact.

Read task, supplied fits/exports, x19 atom construction, low-level spectral
solver, repository classical/moment definitions, and x27/x29 endpoint
results. Applied using-superpowers, systematic-debugging for the direct
control discrepancy, and verification-before-completion guidance.

Implemented an independent unrolled three-step evaluator and a direct
Chebyshev-transform interpolation, sharing only existing low-level grid/
inverse/projector primitives. Computed all 1,000 points at n=64 and n=80,
256 bits. Max absolute n64/n80 path difference 1.228e-51; peak evaluator
RSS 838,598,656 bytes. Own n24 diagnostic verifies interpolation on a
polynomial and polarization against full-Z sfock at one point.

Detected an exact supplied-column issue: (classical(1)-classical(-1))/2
contains e3 as well as e1. A raw sum comparison therefore misses
-s*(5s^2-3)/36. After subtracting that explicit cubic contribution, the
fixed -1 coefficients give max absolute error 9.388e-64, held relative
2.766e-63. The claimed span identity survives. Supplied fit scripts' printed
CONTROL additionally sums unrelated general moments; it is not the
physical decomposition control.

Base QR checks (800/200 split): C held 7.816e-76 in classical d3; expanded
d3 sum held 8.232e-63, A/B 3.294e-24/2.097e-24. At expanded d4, 604 columns,
A/B held 2.487e-29/1.583e-29 at threshold 1e-60; rank 500. Rank thresholds
1e-40 to 1e-60 give rank 491 to 500 at d4, with strongly varying exact-control
residuals. Peak base fit RSS 1,203,175,424 bytes. These reproduce the finite
fit distinction, not a certified exact exclusion.

All115 degree-one and simple-pole coefficient searches use an independent
interlaced 800/200 split, 460 columns. They find no identity; narrow coefficient
classes and sizable extrapolation amplification limit the interpretation.
Peak RSS 979,501,056 bytes. Broader nested denominator tests are tabulated
in x30_report.md and saved in x30_denominator.out.

Exact parameter-linearity proves that merely increasing the b,y parameter
degrees of c0,c2 does not add generators over a rational-function coefficient
field. No exact rational nonmembership certificate is known. The requested
unrestricted A-membership decision and exact obstruction are OPEN. Physical
C1 has an available exact source through the known psi30^(1) formula;
that does not separate A and B.

Kernel conclusion: even with A+B exact, step (2) still needs the evaluated
regular endpoint-to-log source kernel K32 (and upper-log terms/start flux or
framed coupled transport). x29's literal raw average fails the delta endpoint
condition. An exact aggregate source suffices formally, not to evaluate the
missing kernel or recover A individually. B_A and Conjecture A stay OPEN.

Run notes: the first evaluator launch was interrupted to move the repeated
Chebyshev transform outside point evaluation; the corrected computation
completed. The first base fit was interrupted to investigate the nonzero
raw sum residual; the unrolled evaluator was checked against sfock, then
the cubic contamination was identified, and the complete base fit rerun.
A common-product-denominator trial was interrupted and replaced by nested
r1/r12 denominator spaces that retain the known classical d3 controls.
No interrupted run supplies a reported positive result. System python lacks
sympy; the exact checker uses the installed Sage environment's Python with
bytecode disabled. All scientific jobs ran serially. Final commands and
resource/audit evidence are in the report and saved outputs.

Final denominator runs completed: d4/r1 A/B held (cleared target)
7.157e-27/4.773e-27; d4/r12 4.191e-29/3.313e-29. Both retain C/sum controls
at <=1.52e-61. Peak RSS 1,332,002,816 bytes, under the 1.5 GB cap. Exact
checker rerun: all four identities PASS. Shared-log append and prefix/input
audit are performed as the final write/verification steps.

FINAL SUMMARY — PROVED polarization and coefficient-linearity; VERIFIED C and corrected sum (-1 parts), A/B finite residuals; CO-SIGN no known A formula, DISSENT exact exclusion; unrestricted A membership/obstruction and evaluated K32 OPEN; Conjecture A OPEN.
