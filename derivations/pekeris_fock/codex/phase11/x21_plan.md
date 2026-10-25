# TASK_codex_21 research and verification plan

Spec: `TASK_codex_21.md`; read-only context: COLLAB phases 10/11 and x20_report.
Execute inline, serially; all outputs x21_* here, shared logs append-only.
No agents, installations, repository commits, energy jobs or stencil builds.
All explicit writes, including temporary files and Python caches, stay here.

- [x] Read the handoff and actual Green/panel formulas; snapshot shared prefixes
  and relevant read-only input hashes.
- [x] Audit the cyclic D-module exact sequence, solution duality and the
  distinction between a constituent and an endpoint submodule/quotient.
- [x] Prove and independently check the spherical lowering identity for the
  closed K3 measure; test whether P-harmonic sections can exist.
- [x] Compute the local indicial polynomials and convergent analytic solution
  count for the supplied L5. Distinguish one analytic branch from all branches.
- [x] Identify the rigorous bare-surface weight pieces and logarithmic fibre
  weights. Compute the actual density's double residues, with panel
  coefficients, and audit involution/exchange cancellation claims.
- [x] State exactly which full relative-cohomology and physical pairing steps
  remain open; give a bounded next computation and honest cost estimate.
- [x] Run the small exact checker, inspect every result, append COLLAB/THEORY
  entries, verify prefix preservation and input hashes, and report status.

The mathematical result may be OPEN: task item 3 explicitly requests the
cheapest decisive computation if routes 1/2 do not settle Conjecture A.
Neither a nonzero ambient residue nor a failed bounded ansatz is to be
reported as physical survival/non-survival.

Execution note: the task is a research audit, so no feature/TDD, worktree,
commit, or delegated implementation workflow is needed. The user specified
serial execution and the output directory. Scientific identities were proved
and checked inline; the full framed GrW calculation remains explicitly OPEN,
as permitted by task item 3. The optional OS time wrapper failed a sandbox
sysctl query; a clean direct rerun supplies timing/RSS and empty stderr.
