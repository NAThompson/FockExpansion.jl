# TASK_codex_19 serial work record

Scope: execute the three numbered items in TASK_codex_19.md. Only x19_* artifacts and append-only COLLAB.md / THEORY.md Log writes. No agents, no recurrence construction, no energy solves; single-thread serial jobs, peak RSS below 2 GB. Imports suppress bytecode; Julia disables startup/history/precompilation and uses an in-directory writable depot if needed.

1. Independently transport closed41 coefficients into the Codex five-atom distance connection; verify the homogeneous Laplacian exactly and determine/check the two purity constants independently. Compare held-out spectral data.
2. Evaluate all four x16 23-state families by their defining rational integrals, together with order-two and 18 classical atoms. Fit C1/C2/C3 on disjoint training/held-out points versus degree. Validate quadrature, rank/conditioning and denominator-multiplied controls. Do not infer nonmembership from finite-degree residuals.
3. Derive line-charge/Feynman reduction with explicit radical eliminations, check them symbolically, classify the exhibited geometry separately from essential periods of the physical output.
4. Write statuses and co-sign/dissent in x19_report.md, append shared logs, verify artifacts and ownership.

User's explicit fixed directory/serial restrictions supersede worktree, delegation and generic software implementation workflows. This is a mathematical audit rather than a change to a software product.

Progress: task inputs read; no scientific conclusion yet.

Progress 1: exact six-component covariant PDE identity and homogeneity passed (53.5 s, 172 MB). All integer coefficient relations independently reproduce the pickle. Independent square-domain purity quadrature n=16/24/32/48 converges to ah=0.00082106684957556469415831715509798526044, bh=0.01992734793080898549051527717103411701282. The spec41 coefficient fit includes the nonzero hierarchy harmonic part; it is not the pure inverse. This distinction is to be made explicit in the report.

Quadrature diagnostic: initial 128/256 Float64 panel comparison failed the 1e-10 guard (1.166e-9). No fits consumed those values. Investigate the worst input/state and increase quadrature only after a convergence sequence, rather than loosening the guard.

Final coverage:
- Item 1 complete: exact independent connection identity, relation reconstruction, independent purity and held-out spectral/Feynman checks.
- Item 2 complete as bounded numerical membership tests: all 115 atoms, degree sweeps, held-out points, explicit denominator controls, conditioning controls; added full degree-two singlet high-precision test with complete parity reduction. Exact/unbounded membership remains OPEN, explicitly not inferred from residuals.
- Item 3 complete as explicit corrected integral representation and concrete radical eliminations: auxiliary K3 proved generically; essential-period classification and removal of source-density parameters remain OPEN. No assertion that a requested but unjustified two-parameter full-density formula exists.
- Final self-review (serial user constraint; no delegation): checked charge columns, complete polynomial/singlet spans, harmonic/pure distinction, denominator control divisibility, source-dependent density, discriminant twist and generic minimal fibre counts. Corrected report to distinguish numerical failure from exact finite-degree exclusion. All final script outputs pass; shared logs appended only. Largest peak RSS 1,286,881,280 bytes.
