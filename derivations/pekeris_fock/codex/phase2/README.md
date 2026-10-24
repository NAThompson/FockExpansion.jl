# TASK_codex_3 independent results

Read `codex_atoms.md` for A, `codex_order2_recurrence.md` and
`minimal_rank.md` for B, and `codex_minimality2.md` for D.
No Claude phase-two or numerical harness files were read.

The smallest order-two rational differential module has five generators.
For chi use (1,L,alpha,V,U); for the physical log(sigma) coefficient use
(1,Ltilde,alpha,V,U). A six-generator audit system also retains chi and
the log-convention correction separately. This distinction is explained
in the derivation and prevents omission of the angular log correction.

`codex_recurrence.json` is the complete six-atom audit stencil (large,
about 79 MiB because every coefficient is fully expanded).
`minimal_sigma_corrections.json` plus the explicit delete/rename/pi-scale
prescription in the derivation gives the entire five-atom physical
recurrence without a second copy of all unchanged polynomials.
`counts.md` and `counts_minimal_sigma.md` list every block count and range.

Reproduce from this working directory, with bytecode writes disabled:

```sh
PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase2/order2.py
PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase2/recurrence.py
PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase2/minimal_sigma.py
PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase2/verify_polynomials.py
PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase2/verify_atoms.py
PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase2/verify_operator.py
PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase2/verify_stencils.py
PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase2/verify_stencils.py --minimal
```

The saved `.out` files record the successful runs. No proof of finite
classical differential closure through order four, no new variational
state calculation, and no cross-comparison is asserted. The requested
later task can perform sync points 3/4 against the saved independent data.
