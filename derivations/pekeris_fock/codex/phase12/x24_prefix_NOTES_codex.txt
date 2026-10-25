# Codex working notes — 2026-10-21

## Status and conventions

Task 1 workstream A: independent derivation completed before reading sym/. See codex/codex_literature.md for the primary-source audit, codex/codex_crosscheck.md for verification, and COLLAB.md for qualified ansatz agreement. Statements below are **derived sketches**, not a proved minimality theorem. Physical hyperradius is R=sqrt(r1^2+r2^2), sigma=2 eps(r1+r2), eps>0. Coupling lambda remains symbolic in code; the physical case is lambda=1.

## Shell grouping and dimensions

Let a_N=(A_lmn : l+m+n=N), using a fixed ordering. The unsymmetrized dimension is d_N=(N+1)(N+2)/2. Electron exchange r1<->r2 is u<->v, so singlets obey A_lmn=A_mln. The number of fixed triples is f_N=floor(N/2)+1. Thus

    d_N^+ = (d_N+f_N)/2 = floor((N+2)^2/4),
    d_N^- = (d_N-f_N)/2.

A normalized symmetric basis uses (L_l L_m+L_m L_l)/sqrt(2) for l<m and the unscaled basis for l=m. This distinction matters when comparing Euclidean block norms; merely eliminating duplicate indices changes the norm weights.

The exact recurrence is sum_(h=-b)^b B_(N,h) a_(N+h)=0, with rectangular B_(N,h) of shape d_N by d_(N+h), and unavailable negative-index coefficients zero. Counts and shell ranges, from codex/codex_counts.json:

| operator | general offsets | general h | at E=-eps^2 offsets | fixed h |
|---|---:|---|---:|---|
| P | 57 | -3..3 | 33 | -2..2 |
| sigma^2 P | 225 | -5..5 | 177 | -4..4 |
| C1 | 117 | -4..4 | 117 | -4..4 |
| C2 | 57 | -3..3 | 57 | -3..3 |

The infinite log hierarchy couples sector j only to j,j+1,j+2. At each N, including every log sector makes the shell space infinite dimensional: the local condition n_radial>=2j is a condition on corner jets, **not** a cutoff j<=N/2 in Laguerre index N (every Laguerre polynomial has a constant term). A finite K computational truncation is an approximation and needs an error model; it is not an exact finite tower.

## What minimality could mean

Scalar ratio minimality means f_N/g_N ->0 against an independent dominant solution. A vector-norm analogue is basis dependent and generally cannot uniquely select a vector when several decaying channels exist. With growing blocks, the outer B matrices are rectangular and need not be invertible; one cannot simply solve for a_(N+b) as a fixed-dimensional transfer recurrence. Grouping b consecutive shells makes a block-tridiagonal pattern, but the grouped dimensions still grow and off-diagonal blocks remain rectangular. Counting bandwidth is not a Pincherle theorem.

A defensible first definition is the subspace of coefficient sequences representing solutions in the physical operator domain with finite physical norm and the specified collision conditions. At a simple ground-state eigenvalue this physical solution space is one dimensional, after the PDE/recurrence equivalence and convergence are established. This is a spectral admissibility definition; identifying it with an asymptotically subordinate or backward-stable subspace is work still to do. For the full tower one must also fix its gauge and a topology ensuring summation, differentiation, and physical-norm convergence. Only the sum Psi need be physical; each sector need not separately satisfy every domain condition.

A possible rigorous route is to form the physical symmetric generalized eigenvalue problem in the Laguerre-function basis, establish a closed lower-bounded form, and use a tail resolvent/Dirichlet-to-Neumann map for nonreal spectral parameter. Finite-dimensional Schur complements give a concrete matrix continued-fraction construction even with rectangular inter-shell coupling. Prove convergence by form/core or resolvent arguments, then take the boundary limit toward an isolated real eigenvalue. At an eigenvalue the resolvent has a pole, so the physical eigenvector is extracted by residue/kernel, not by assuming every real-energy inverse exists. An energy-dependent choice eps=sqrt(-E) also changes the basis with E: fix eps while constructing the operator pencil, or explicitly track this dependence.

For the raw polynomial P recurrence, multiplication by the positive interior factor Q has changed the representation and is singular at collision boundaries. In fact P=-Q exp(sigma/2)(H-E)exp(-sigma/2) is formally symmetric for exp(-sigma) du dv dw when boundary terms vanish, because the physical measure is proportional to Q du dv dw. This formal symmetry is useful but does not settle the closure or domain. Formal bandedness alone does not imply that P is a self-adjoint Jacobi matrix in ordinary coefficient l^2. Domain equivalence, boundary forms, and the proper metric must be established before importing deficiency-index conclusions.

## Relevant theory and precise applicability limits

- **CITED primary source:** W. Gautschi, *Minimal Solutions of Three-Term Recurrence Relations and Orthogonal Polynomials*, Math. Comp. 36 (1981), 547–554, Sec. 3, Eqs. (3.1), (3.5)–(3.8), [author PDF](https://www.cs.purdue.edu/homes/wxg/selected_works/section_04/075.pdf). Gives the scalar ratio definition and identifies a minimal solution using a determinate moment problem, outside the supporting interval. These are scalar three-term hypotheses, not a proof for this growing-block problem.
- **CITED abstract only:** P. Levrie and A. Bultheel, *Matrix continued fractions related to first-order linear recurrence systems*, ETNA 4 (1996), 46–63, [journal abstract](https://etna.ricam.oeaw.ac.at/volumes/1993-2000/vol4/abstract.php?pages=46-63). Establishes a matrix Pincherle-type result for finite-dimensional systems; exact theorem hypotheses/number have not been inspected. No extension to the present rectangular growing blocks is claimed.
- **CITED abstract only:** V. S. Budyka, M. M. Malamud, K. A. Mirzoev, *Deficiency Indices of Block Jacobi Matrices: Survey*, Contemporary Mathematics. Fundamental Directions 67 (2021), 237–254, DOI [10.22363/2413-3639-2021-67-2-237-254](https://journals.rudn.ru/CMFD/article/view/28864/en_US). Its stated setting is m by m blocks with fixed m. Relevant operator questions are self-adjointness, deficiency spaces and boundary conditions at infinity, but no numbered theorem from it is asserted here.
- **RECALLED-UNVERIFIED theorem details:** J. Wimp, *Computation with Recurrence Relations*, Pitman, 1984, is the standard recurrence/minimal-solution reference; exact chapters/theorems were not inspected. Do not infer a growing-block theorem from that citation. Classical Weyl limit-point/limit-circle and block-Jacobi theory are relevant analogies only until their hypotheses are matched.

## Laguerre l^2 versus physical normalizability

Plain L_n are orthonormal for exp(-x) dx. Hence Parseval gives

    sum |A_lmn|^2 = integral_(octant) exp(-sigma) |F|^2 du dv dw,

whenever F denotes that L^2 sum. The physical S-state volume is 8 pi^2 r1 r2 r12 dr1 dr2 dr12. The absolute Jacobian is 1/(8 eps^3), while r1 r2 r12=Q/(32 eps^3), Q=(u+v)(2u+w)(2v+w). Therefore

    ||Psi||_physical^2 = pi^2/(32 eps^6)
        integral exp(-sigma) Q |F|^2 du dv dw.

These norms are not equivalent: Q is unbounded at infinity and vanishes at collision boundaries. For an angular cutoff supported in the cone interior and sigma=rho, F~rho^(-2) near zero has finite physical norm (integrand rho^(5-4)) but infinite unweighted norm (rho^(2-4)). After a smooth cutoff away from zero, F~exp(rho/2) rho^(-2) at infinity has finite unweighted norm (rho^(-2)) but infinite physical norm (rho^(1)). These are norm counterexamples, not Schrödinger eigenfunctions. The physical domain and Coulomb regularity may exclude them, but l^2 alone does not encode that domain.

For finite polynomial expansions the physical Gram matrix is exactly multiplication by Q in the orthonormal Laguerre basis. Thus the physical completion is the closed multiplication-form domain, formally A^*G A<infinity, with the appropriate interpretation beyond unweighted l^2. A convenient sufficient coefficient condition is sum (1+l+m+n)^3 |A_lmn|^2<infinity: Q<=C sigma^3 and the cubic multiplication form has bandwidth three and entries O((1+N)^3), so weighted Schur bounds apply first on finite sequences and then by closure. Necessity is not claimed. For log sectors cross terms in log(sigma) couple their Gram forms; separate sector l^2 bounds do not by themselves control the sum.

## Fock constants as functionals: preliminary obstruction and route

With pointwise/derivative convergence at the corner, F0(0)=sum A^(0)_lmn since L_n(0)=1. Bare l^2 does not make this sum convergent; evaluation is not bounded on that Hilbert space. Likewise Taylor jets use L_n^(q)(0)=(-1)^q binomial(n,q), requiring stronger coefficient weights or a justified Abel limit. Since actual homogeneous angular coefficients need not be polynomial, a21 should be obtained by subtracting the known orders and log term, taking an R->0 angular projection onto the degree-2 harmonic in the fixed Liverts convention, and dividing by Psi(0). The log(sigma)-to-log R angular correction and exp(-sigma/2) lower-order corrections must be included. Prove that the angular projected finite-part limit is a convergent functional on the selected solution class; a formally written coefficient sum is insufficient. No numerical value or full a21 functional is claimed at this stage.

## TASK_codex_2 — 2026-10-21

Stage 3 is now written in codex/codex_minimality.md with labeled proofs and conditional statements. The physical form-domain admissible-subspace definition is sound; a growing-block ratio-minimality theorem is not established. Finite nonreal Schur complements are valid, while convergence needs a form-core hypothesis and care with raw coefficient projections. Shell-norm computations at cutoffs 12,16,20 show decay but strong boundary distortion, not a tail law.

The optional global a21 calculation was completed independently on all exported N=6..13 states, with manufactured checks and three angular quadratures. N=13 global is 0.476355920012 (local 0.472445423120), using E_variational. See codex/codex_a21_functional.md and scripts/output there. Trial linear cusps are not exact; finite-part residual conventions matter. No exact a21 digits are certified.

REVIEW_2_codex.txt reviews the numerical source, tables and plot; COLLAB.md has an append-only claim-by-claim co-sign/dissent. Main dissents concern automatic decay, infinite-tower top-down solvability, inferring arbitrarily high nonzero logs from psi42 alone, and unproved global/form-domain completeness. No edits to Claude's files.

Concurrent Claude log-squared update was also reviewed (REVIEW_2_codex.txt section 8). The p<=4 Bell polynomials passed independent Gamma-derivative checks. The zeta EM precision comment is optimistic (worst error ~7e-51, not 1e-60); finite log-squared energy improvement is supported by the saved data, not a proof of all-order logarithmic necessity.


## TASK_codex_3 — 2026-10-21 — independent phase-two results
A/B/D are complete in codex/phase2/; start at README.md. Minimal order-two
rational differential rank is five for chi and five for the physical
log(sigma) combination with Ltilde; six generators are needed only when
chi and the convention log are retained independently. Monodromy proof
and exact generation witness are in minimal_rank.md. All 31 audit stencil
blocks plus both minimal replacements pass exact independent projection
checks; three actual-polylog Hamiltonian actions pass. Order-four hard
content is integral-defined, not a proved minimal classical alphabet.
Admissibility/minimality remains the same physical quotient formulation.
No claude2/ or julia/ input; no repository changes. Full file inventory and
verification figures are in the appended COLLAB entry and REVIEW_3_codex.txt.

### 2026-10-21 — TASK_codex_4, phase-two sync 4 and review

- Read the complete COLLAB handoff and TASK_codex_4; read Claude's phase-two
  formulas, recurrence JSON and numerical harness read-only. All task writes
  are in codex/phase2/, REVIEW_4_codex.txt, this file and append-only COLLAB.md.
- Basis relation: L=la, alpha=al, V=(beP+pi P/2)/r12,
  U=Om+pi lb, Slog=-lR. All 144 connection transformation entries vanish
  exactly; actual functions agree at both distance orderings.
- Independent saved-stencil/direct-Hamiltonian check: 48 actions pass,
  maximum relative error 2.33572329193e-39. JSON really has 16/15/8 blocks,
  25,955 offsets, bandwidth 10; COLLAB's 17/12/8 is a reporting error.
  See codex/phase2/codex_crosscheck2.md and sync4_counts.md.
- The independent Mellin/Clausen script (no Claude source included) reproduces
  all requested atoms2 state functionals; N=11 gives 0.476748379513.
  Quadrature/precision stability is around 1e-14, but physical error against
  the supplied reference is about 5e-7. Recommend about 0.47675 with a
  heuristic 2e-6 uncertainty scale, not certified digits or an interval.
- Independently checked the weak-form measure, metric weights, derivative
  chain, radial Bell moments and a manufactured finite part. Startup-free
  refined-quadrature Ritz reruns use Claude's read-only assembly; N=11 at
  384 bits gives E-E0=2.4314823228526578e-17, positive Gram Cholesky,
  scaled generalized residual 1.2851e-50. Further detail in review/logs.
- Co-sign strong observed atom enrichment gains; dissent from superalgebraic
  convergence, certification of a21 from energy, and certified low-N psi30
  gains without consistent high-precision value/gradient checks. N>=8 atoms4
  is invalid physical numerical evidence. Minimality claims remain conditional.
- Final precision addendum: 512-bit h=1/32 N=11 E-E0=2.4314823230015396e-17,
  about 1.49e-27 from h=1/24; its rebuilt state's functional agrees through
  15 printed decimals. `audit32/` contains the independently recomputed states.
  A concurrent Claude run updated two quadrature logs and the N=11 export;
  final input-state copies and hashes in phase2 fix the functional's inputs.

### 2026-10-22 — TASK_codex_5, phase three parameter-connection handoff
Derived a new finite outer-variable connection for the psi30 elliptic panels,
with explicit root/rational entries, initial data and five-jet extraction.
See codex/phase3/README.md, structure.md, connection.md, reachability.md and
REVIEW_5_codex.txt. Interior connection, moment, full-repository and exact
homogenization checks pass; detailed limits and logs are saved. Closure bound
216 over an algebraic coefficient field is NOT a minimality theorem or rational
rank. Global Laguerre counts/bandwidth and minimal combined module remain
unresolved; task is explicitly partial. No Claude-owned files were edited.

### 2026-10-22 — TASK_codex_7 phase-three review
Completed the requested source/table review and dated per-claim co-sign/dissent
in COLLAB.md; full review is REVIEW_7_codex.txt. Independent new code (no Claude
includes/matrices/functional) integrates exported states directly: min5 N8
energy error 2.79318618e-16, agrees within 8.1e-25 hartree; N4 within 8.9e-26.
Independent min5 N8 a21=0.47673849249439, angular orders 24/48/96. Exact cubic
homogenization plus build3 degree bounds prove finite-span E/a21 independence.
Parity checks pass. Last min5 log-log slopes -12.55,-13.83,-15.14,-10.52;
semilog sampled curve remains convex. Honest physical a21 estimate about
0.47675, empirical few-1e-6 uncertainty, no guaranteed digits/enclosure.
Phase-three closure remains local finite, minimal rank/global interfaces open;
order-four reachability remains a qualified extension route. Evidence and
reproduction commands: codex/phase3/review7_*. No Claude-owned file modified.

### 2026-10-22 — TASK_codex_8: independent Fock-constrained Ritz audit
Independent seed and polynomial-only Ritz implementation in codex/phase4/ was completed and seed/N3/N4 results saved before reading phase-four Claude source. N3/4/5 (8/17/29 functions), eps=1.704, give E-E0=9.83205e-4,1.16740e-4,2.14829e-5 and direct a21=.711314572436,.643915002059,.607208872439. Exact symbolic jet cancellation, source angular checks <3.2e-58, Gram/shifted-H Cholesky and residuals <7e-52 pass. h8/h12 refinement changes energy <1.2e-20 and a21 <1.1e-18. Own global correction gives .521582910466,.503264783993,.494494882114. Small serial jobs, no atom cache.

Critical new finding: bare V=sqrt(1+x²-y²)/y*(beta+pi/2) has V~pi/y at electron coalescence. build4 retains a=b=0 V terms, and N8/N10/N11 exports contain nonzero coefficients (N8 ~1e-34). Literal original states therefore have infinite kinetic form. N4 has no V enrichment. Original N8 cutoff weak-form integration agrees with its saved energy to 8.81e-26, but is not a physical Rayleigh quotient. Excluding the pole terms exactly in Codex's reconstructed copy gives the same displayed N8 energy; this is not a corrected Ritz minimization. Claude should remove those basis functions or regularize V and reassemble. No Claude file was edited. Evidence: pole_issue.md, pole_audit.py/.out.

Independent global estimates for N8/10/11 reproduce .476747842992872,.476747890216274,.476747886340663; Gauss24/48 agree within2.8e-15. Dropping bare-V terms changes them <1e-30. Direct N11 remains .476749180410913. Residual-correction identity holds for collision-admissible constrained states; bare poles require collision distributions/boundary terms. No certified physical digits follow. Corner a21 is unbounded in physical form norm: c D eta(R/h) changes its jet by c but costs O(c²h8) energy. Density alone does not prove readout consistency or any energy-error rate. See theory.md and REVIEW_8_codex.txt. All writes stayed within the task allowlist; snapshots/hashes and fresh checks are in codex/phase4/verification.out.

### 2026-10-22 — TASK_codex_9: plain ground-state backward selection
See codex/phase5/goal1_theorem.md, goal1_sources.md and goal1_review.md. Supplied independent proofs of fixed-eps polynomial form-core density, fixed-head coefficient continuity in form norm, and coercivity of q-E0 on ker A000 using positivity/nonzero normalization and the ground spectral gap. These prove form convergence of the actual normalized zero-tail/deleted-row Miller finite sections at exact E0, with vanishing dropped residual and best-approximation/defect bounds. Nonreal constrained solutions converge but retain a source; the head residual must vanish for an eigenfunction. Physical resolvent, ground Ritz/projection and nearby Schur-root convergence follow; shifted inverse iteration has a gap-dependent iteration rate. Nested physically admissible log/atom enrichments containing the plain core inherit physical Galerkin convergence, not automatically a gauged formal-sector Miller theorem. Shell-cutoff rates, all-formal-solution classification and ordinary limiting-series form convergence remain open.
Independent goal1_verify.py uses only the saved Codex recurrence, 45 digits, N2/4/6 (max dimension 50), serial; positivity, Hermiticity, Schur/direct and normalized-tail checks pass, max Schur residual <1.2e-45, peak RSS ~66 MiB. Reproduces Claude N6 raw residuals. Claude backward/tail scripts and outputs reviewed read only; exact source and finite-evidence qualifications in goal1_review.md. Klahn--Bingel I/II abstracts checked (II Theorem 10 locator), full hypotheses inaccessible; Klahn--Morgan DOI located, full text unavailable. Shared COLLAB/THEORY logs appended only; no external or Claude-owned writes.

### 2026-10-22 — TASK_codex_10 (phase 5 Goal 2)
- Read COLLAB.md, TASK_codex_10.md, THEORY.md and the saved own connection/structure. Files added only as codex/phase5/goal2_*; shared COLLAB.md and THEORY.md Log appended only. No Claude-owned files or files outside this directory were written. Numerical jobs ran serially, small matrices, BLAS one thread; no large Ritz job.
- VERIFIED stabilized represented candidates: classical 12; after merging duplicated zero jets and algebraic/angle endpoint states, each of four panels N0/N2/joint=11/12/13 at two exact finite-field points, stable orders 4–8. Combined {1,chi,R^3 psi30} candidate=45 at (t,q)=(3,4),(4,7), exact ranks 3,7,11,15,20,26,33,41,45 then 45 through order 14. Both full 332-state constant-curvature checks vanish exactly. Explicit identical 45-row selections are in goal2_candidate_generators.json and goal2_results.md. Padded generic panel counts 15/18/19 were independently checked by 160-bit Cauchy jets at two (b,y) points, but are superseded for pruning by the quotient counts. Double precision alone was insufficient at high orders.
- VERIFIED port/actual-function checks: 312 panel connection actions and four moment outputs <4.1e-16; 72 classical derivative components/two classical outputs <6.2e-16; actual 332-to-240 endpoint substitution at two angular points <1.2e-16; independently differentiated old-L/U pullbacks <8.3e-17. Final audit passes: goal2_verify.out. Read-only input hashes saved.
- PROVED exact radical reduction: z2=(r2/r1)z1, z4=i*c/(P*z3), with compatible root-set signs. Sufficient field K(R,P,A,z1,z3) has degree <=32; old proven rational closure bound improves to <=6912. This is not a minimal field degree or rational rank. Branches R,P>0,A=i*sqrt(-D), first-quadrant z1,z3 are globally analytic in the strict real triangle.
- AGREE conditionally with global derivative atoms: assuming actual functions globally analytic in the connected strict triangle, choose a rational basis of their finite local derivative span; cleared identities extend globally by the analytic identity theorem. Simple connectivity is unnecessary given global functions. This bypasses old panel presentation interfaces, but proves neither the premise for the assembled psi30 nor collision-domain admissibility of high derivatives.
- OPEN: generic symbolic closure/actual-function independence for the 45 candidate, hence strict smallest rank/minimal generator claim; rational minimal rank; domain restrictions on high derivative atoms. Claude currently supplies goal2_edges.jl/.out only, not a recurrence: finite bounded values/gradients/Hessians on four sampled paths do not prove uniform boundary control, and the requested later recurrence spot-check remains pending. Full report: goal2_results.md; global/radical arguments: goal2_global.md. No task-complete claim for these remaining strict requirements.

### Codex — 2026-10-22 — TASK_codex_11 enriched selector and classical comparison
- **PROVED (supplied arguments):** all build3 order-two products with the exact yV regularization are in the physical H1 Coulomb form domain. Ltilde, alpha-pi/2, yV and the odd global U have bounded angular values and Cartesian gradients O(1/R). For U, half-angle reduction gives dU=Fa(a)da+Fb(b)db with removable endpoint poles and zero exchange trace; this controls every collision/collinear stratum, not just numerical paths. Finite log/log² products and the corner/infinity are covered.
- **PROVED under the phase-5 simple isolated positive ground-state input:** enriched physical weak-pencil zero-tail solutions converge in form norm with a bounded positive overlap head, centered/quotiented physical coordinates and all nonhead weak rows. Tail coercivity, vanishing dropped source, nonreal Schur convergence, ground Ritz/nearby omitted-row zeros and best-approximation/error-defect bounds follow. Literal build3 replaces the constant seed; its positive-degree polynomial subspace is itself a form core (point evaluation is not H1-bounded in six dimensions), so its density is proved separately. No unproved seed-coefficient normalization is used.
- **Classical qualification:** primary Leaver 1985 (8),(9),(25) and 1986 (39)–(46) show repeated leading root 1 and stretched-exponential separation in standard radial Jaffe/Leaver problems. Pincherle minimality does not require geometric dichotomy. Fock nonanalyticity excludes geometric decay of a reconstructed full plain Laguerre array, but neither proves an algebraic exponent nor rules out separation from growing solutions. No basis-independent helium no-dichotomy theorem or restored Sturmian/Fock-basis dichotomy is established. Exact sources and proved/conjectural distinctions: goal1x_jaffe_leaver.md.
- **DISSENT on the phase-6 inference from branch guessing:** guess_probe's first six nonbasis order-8/9 targets contain only sources 1,2 (constant and chi), no new f=R³psi30 relation. guess_line's three targets likewise omit f. The first omitted f derivative is (src,a,b)=(3,4,4). The reported branch invariance therefore does not test the new order-three rational boundary relations or their degrees. Read-only review and a small exact key-enumeration check are supplied; no Claude file was modified.
- **VERIFIED, small serial checks:** U half-angle/J and derivative identities at 60 digits (error <4.7e-61), exact removable endpoint limits, exchange trace samples, exact 3x3 weak-tail/Green identity, algebraic-minimal/geometric-dominant counterexample and target enumeration all pass. No shell jobs, large matrices, agents or heavy reconstructions were launched. Files: codex/phase6/goal1x_theorem.md, goal1x_jaffe_leaver.md, goal1x_review.md, goal1x_verify.py/.out and audit manifests; NOTES_codex.md and shared logs appended only.
- **OPEN:** explicit N-rate beyond best-approximation transfer; enriched formal cleared-recurrence/reconstruction-gauge identification; universal formal-solution classification; order-three derivative-atom admissibility. No a21 digit chasing.

### 2026-10-22 — TASK_codex_13
See codex/phase7/report.md and global_and_domain.md. New exact sufficient
connections: classical 18-state degree (2,3), panels 55-state rational (h,z)
cover degree at most (6,7), distance coefficient field degree <=16.
The physical N0,N2 integrands are antisymmetric under s->h/s, proving a
fixed-endpoint-1 global analytic presentation; old min-endpoint primitives
have jumps and must not be confused with it. Classical atoms are unweighted
form-admissible; explicit conservative W^18 weights suffice for all panel
atoms. Independent line audit and two fresh exact points confirm panel
47/49 degrees; actual minimal ranks remain unproved. Independent {1,a,a²}
Hamiltonian recurrence has 563/415/251 distinct offsets, bandwidth 7;
all exact action/product-rule checks pass. Full elliptic recurrence/Claude
cross-check and sharper unweighted domain estimates remain open. All jobs
small and serial; no Claude-owned or outside-directory file edited.

### 2026-10-23 — TASK_codex_15, rational panel 3/4
PROVED exact real 55-state connection over Q(r1,r2,r12), with no variable
radicals in its coefficients. Rebase quartic pole jets to s^j/T integrals,
real-pair the angle poles, normalize density by c/R and replace D jets by
h^n D jets. The new B=h atan(s/h) has a rational integral definition.
Panel 4 is electron exchange. codex/phase8/report.md gives the atom list,
closure/globality proofs and exact N0/N2/R^3psi30 output rows. Exports contain
673 distance entries with per-entry degrees; sharp joint row clearings
have maximum degree 11. Independent h,z gauge identities, exact Hermite
certificates, actual-integral derivatives, physical panel outputs and exact
rational-point pullbacks all pass (verify.out, audit.out).
PROVED obstruction: Q0_A0~27sqrt(3)/(16epsilon) at x=1,y=2,c=1+epsilon;
it is not locally L2 at the collinear boundary. A sufficient rational angular
weight (Delta/Q^2)^4 makes every new atom form-admissible without changing
radial degree. Minimal weight and exact unweighted angular-jet inclusion in
a polynomial span of weighted atoms remain OPEN. No panel1_* artifact is
available in claude8, so its conditional cross-check is PENDING. No Claude
file/outside-directory file modified; all research computations small,
no agents or energy/stencil jobs.
