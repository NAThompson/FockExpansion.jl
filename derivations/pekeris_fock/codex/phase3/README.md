# TASK_codex_5 handoff — 2026-10-22

**Status: substantial partial result; not all requested claims are established.**

- `structure.md`: general Z,E,a21 distance decomposition of R³ psi30,
  polynomial part, classical alphabet and four panel maps.
- `connection.md`: explicit finite outer-variable Pfaffian connection with
  roots/radicals, all entries, normalization and moment extraction; also the
  exact boundary/interface limitations. No unspecified moment derivatives.
- `panel_connection.py`: executable new connection and independent integral evaluator.
- `classical_reference.py`: attributed copy of the older workspace classical
  formula. Its reuse is deliberate; it is not an independent classical derivation.
- `reachability.md`: reasoned order-four extension route and harmonic data.
- `verify_connection.py/.out`, `connection_checks.json`: 312 component equations
  per step size, 42-digit arithmetic; moment reconstruction ~1e-43.
- `verify_flatness.py/.out`, `flatness_checks.json`: entire padded 65x65
  parameter connection, arbitrary initial vectors; curvature decreases by four
  under halving of a centered-difference step (max 1.12e-16 at step 5e-9).
  Numerical parameter matrix counts 406/333 are NOT Laguerre stencil counts.
- `verify_homogeneous.py/.out`: five exact polynomial homogenization identities.
- `repo_values.jl/.tsv/.err`: read-only repository BigFloat/Float64 references.
- `verify_repo.py/.out`, `repo_checks.json`: all contributions agree at three
  general parameter sets, both distance orderings and Z=0, to about 48 digits.
  Maximum observed Float64 discrepancy 6.06e-15. Repository invocation emits
  one SpecialFunctions/OpenSpecFun binding warning, but returns all values.

Run from the task directory, with bytecode generation disabled:

    PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase3/verify_connection.py
    PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase3/verify_flatness.py
    PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase3/verify_homogeneous.py
    JULIA_PKG_PRECOMPILE_AUTO=0 julia --startup-file=no --history-file=no --compiled-modules=no --project=/Users/nathompson7/FockExpansion.jl codex/phase3/repo_values.jl
    PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase3/verify_repo.py

The Julia command's stdout must replace only `codex/phase3/repo_values.tsv`
before running the final Python comparison. Other checks write JSON beside
themselves. All output files belong to Codex. No Claude source is imported.

## Unresolved task requirements

1. Minimal combined order-two/order-three atom set and a proof of minimality.
   We have a chartwise **upper bound** of 216 over an explicit algebraic field.
2. Rational-field restriction, globally compatible atom branches/interfaces,
   row-cleared global Laguerre stencils, actual offset counts and shell bandwidth.
   Finite chartwise differential closure alone is insufficient to supply these numbers.
3. Order-four explicit minimal connections are a reasoned extension route,
   not computed matrices or a proof of membership in the order-three module.

The task is not labeled complete. The new usable result is the parameter
connection, which the earlier finite word evaluator did not itself provide.
