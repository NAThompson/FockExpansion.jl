# Phase 4 independent constrained Ritz audit

All explicit writes are confined to the TASK 8 allowlist. No Claude implementation is included or executed. The main independent code was written and the seed checked before reading phase-4 source. Sources and qualifications: plan.md and theory.md.

Run in the task directory, serially (one numerical Julia process at a time):

```
venv/bin/python codex/phase4/verify_jet.py > codex/phase4/verify_jet.out
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --startup-file=no --compiled-modules=no --project=/Users/nathompson7/FockExpansion.jl codex/phase4/check_seed.jl > codex/phase4/check_seed.out 2>&1
NS=3,4,5 HINV=8 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --startup-file=no --compiled-modules=no --project=/Users/nathompson7/FockExpansion.jl codex/phase4/constrained.jl > codex/phase4/run_h8.out 2>&1
NS=3,4,5 HINV=12 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --startup-file=no --compiled-modules=no --project=/Users/nathompson7/FockExpansion.jl codex/phase4/constrained.jl > codex/phase4/run_h12.out 2>&1
```

constrained.jl assembles tiny dense matrices with exact radial moments and BigFloat angular quadrature. Only 8,17,29 basis functions at N=3,4,5. Scaled Gram and shifted Hamiltonian Cholesky must succeed. Shifted inverse iteration uses E0 below the finite Ritz energy and prints the generalized residual. No cache is written outside this folder; Julia startup and compiled-module caches are disabled. Dependencies/project in the repository are read only. Refinement differences are empirical quadrature checks, not interval error bounds. The omitted endpoint region is controlled numerically, not by a rigorous remainder estimate.

review_state.jl evaluates Claude's exported trial states and their jets with independent actual-function formulas, derivatives and radial moments. It does not rebuild their variational matrices or independently minimize their richer spaces. corrected.jl implements the radial Mellin finite-part residual correction in theory.md; its Gauss nodes originate in a Float64 eigenproblem, so BigFloat arithmetic does not certify all displayed digits.

Additional serial verification commands:

```
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --startup-file=no --compiled-modules=no --project=/Users/nathompson7/FockExpansion.jl codex/phase4/review_state.jl > codex/phase4/review_state.out 2>&1
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --startup-file=no --compiled-modules=no --project=/Users/nathompson7/FockExpansion.jl codex/phase4/corrected.jl codex/phase4/poly_N3_h12.txt codex/phase4/poly_N4_h12.txt codex/phase4/poly_N5_h12.txt > codex/phase4/corrected.out 2>&1
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --startup-file=no --compiled-modules=no --project=/Users/nathompson7/FockExpansion.jl codex/phase4/global_state.jl > codex/phase4/global_state.out 2>&1
```

Claude state inputs are frozen copies in input_states/ with SHA256 values in input_hashes.json. The final source library removes unused phase-3 driver functions; only letters, Gauss nodes and exact radial moment formulas remain.

The first state check used /usr/bin/time -l: all numerical results completed, but its final sysctl resource query was denied by the sandbox, so the wrapper exited 1. A clean driver rerun without this resource wrapper is recorded separately. No measured peak-memory claim is made.

Fresh clean cross-check (after reducing the shared own-code library), omitting the unavailable resource query:

```
HS=12 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --startup-file=no --compiled-modules=no --project=/Users/nathompson7/FockExpansion.jl codex/phase4/review_state.jl > codex/phase4/review_state_clean.out 2>&1
```

Collision audit and independent pole-removal diagnostic:

```
venv/bin/python codex/phase4/pole_audit.py > codex/phase4/pole_audit.out
DROP_POLES=1 HS=12 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --startup-file=no --compiled-modules=no --project=/Users/nathompson7/FockExpansion.jl codex/phase4/review_state.jl > codex/phase4/review_state_pole_removed.out 2>&1
```

Read pole_issue.md before interpreting the original min5 state integrals: nonzero bare V coefficients cause a 1/r12 pole and invalidate literal physical Rayleigh-quotient claims for the original exports. Their finite cutoff quadratures remain reproducible numerical evidence. No Claude source or export was changed.

```
DROP_POLES=1 NS=8,10,11 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --startup-file=no --compiled-modules=no --project=/Users/nathompson7/FockExpansion.jl codex/phase4/global_state.jl > codex/phase4/global_state_pole_removed.out 2>&1
```
