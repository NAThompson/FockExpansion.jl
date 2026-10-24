# TASK_codex_3 execution record

2026-10-21. Scope: only phase2 files, NOTES_codex.md, REVIEW_3_codex.txt,
and append-only COLLAB.md. Repository is read-only. No claude2 or julia reads.

- A: translate source formulas; distinguish closed expressions and numerical inverses.
- B: simplify the order-two dilogarithmic differential module, construct its
  rational connection and atom/log-sector operators, clear denominators, verify.
- D: carry the physical admissibility formulation through the reconstruction map.
- Record independent results before any later comparison. Cross-check deferred
  as TASK_codex_3 requests.

Workflow ruling: this is a prescribed mathematical derivation in an existing
  shared scratch directory, not a feature branch; no worktree, commits, or
  files outside the explicit user allowlist will be created.

Verification must distinguish exact symbolic checks, numerical checks, and
unresolved minimal-rank or flattened-stencil claims.

Completed: A source audit and exact homogenizations; B rational closure,
expanded stencils, direct and exact projection checks; D admissibility note.
Minimality of the atom module is proved separately by analytic monodromy
and the exact generation determinant in minimal_rank.md. The physical
log(sigma) combination uses five atoms; an independently representable
chi plus Slog uses six. All writes stayed in the user allowlist; no
claude2/ or julia/ files were read. Later cross-comparison remains deferred.

Implementation verification notes: the polynomial-ring serializer required
normalizing SymPy symbols with assumptions; the exact checker required a
higher recursion limit for long expanded polynomial ASTs. Root causes
were read from the tracebacks; the final complete projection checks pass.
