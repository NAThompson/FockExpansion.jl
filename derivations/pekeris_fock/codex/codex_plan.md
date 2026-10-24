# Independent Pekeris workstream plan

Spec: TASK_codex_1.md. User requests execution now; all writes remain in codex/, codex_*, NOTES_codex.md, and append-only COLLAB.md.

1. Assess every ansatz element; append qualified agreement/dissent and state the required function class.
2. Independently review primary literature, recording theorem/equation precision and verification limits in codex_literature.md.
3. Derive the distance-coordinate Hamiltonian by chain rule, transform and conjugate, and compute Laguerre matrix elements in codex_derive.py. Save equation-index recurrence coefficients in codex_recurrence.json and counts/bandwidth metadata. Verify polynomial actions, boundary indices, log conjugation, and distant-offset cancellation without reading sym/.
4. Only after step 3, read sym/recurrence.json and compare exact symbolic differences; document all discrepancies and machine checks in codex_crosscheck.md.
5. Sketch growing-shell minimality, symmetry dimensions and physical-weight versus coefficient norms in NOTES_codex.md. Append a dated completion entry to COLLAB.md.

Review focus: equation versus source indexing; mixed-derivative factors; negative indices at Laguerre boundaries; derivative tails cancelling; fixed-energy simplification; finite-log conclusions requiring spectral/domain hypotheses; analytic Cartesian versus distance/conormal regularity.
