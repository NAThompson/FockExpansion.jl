Direct Fock recurrence certificate

With uv installed, run from the repository root:

    ./verification/direct_recurrence/verify_all.py

Or from any directory:

    uv run --script /path/to/FockExpansion.jl/verification/direct_recurrence/verify_all.py

uv selects a compatible Python and installs the pinned inline dependencies in
an isolated environment. No manual virtual environment or pip step is needed.
The first run requires access to download any missing Python/dependencies.
All subprocesses use the same uv-managed Python interpreter.

Expected final messages:
    ✓ Certificate verified

The summary lists all eight recurrence sectors and the coefficient-table links.

The runner writes certificate.json and diagnostic logs alongside itself.
A failed check exits nonzero. This is the algebraic certificate; the analytic
regularity and uniqueness arguments are in the accompanying text.

Without uv: install requirements.txt, then run python verify_all.py.
The runner regenerates every residual and checks all coefficients exactly
over Q(sqrt(2),i)(T,U,pi,G,log(2)), using i only to rationalize cotangents.
No numerical samples or third-order series formula enter the zero tests.
Constants pi,G,log(2) are treated as algebraically independent.

Inputs are snapshots of the compact classical expression, polynomial
integral weights, proved classical conic evaluations, and lower-order chi.
sym_H includes unused historical helper functions; only chi() is used.
No third-order auxiliary assembly is called.

See REGULARITY.txt and the manuscript appendix for analytic justification.
Generated .pkl files and logs are outputs, never trusted by the runner.

Integration identities: integration_proofs contains the exact conic-gradient,
vertex, rational-integral, and arcsine-derivative checks used before the final
residual reduction. Run conic_proof.py inf 1 -1 from that directory for all six
conic gradient identities, then conic_vertex.py and rational_terms_proof.py.
SUPPLEMENTARY_DERIVATIONS.txt retains the detailed recurrence argument.


Native Julia counterpart
------------------------
With Julia installed, from the repository root run:

    ./verification/direct_recurrence/verify_all.jl

Or: julia --startup-file=no verification/direct_recurrence/verify_all.jl
The entry point activates and instantiates the separate julia/ environment.
It does not modify the main package's dependencies. The first run downloads
missing dependencies and precompiles them. The final message is:
    ✓ PASS  Certificate verified · 8/8 sectors
    Report: certificate_julia.json

The Julia code rebuilds the panel differential action, verifies the integral
primitive and rational remainder identities, checks the 18 classical table
links, differentiates the undifferentiated formula, and checks all eight
recurrence sectors using exact arithmetic. Deliberate changes to the formula
and source are required to produce nonzero residuals as negative controls.

Scope of independence: Julia does not run Python or read raw_*.pkl,
flint_*.pkl, panels.pkl, or either prior certificate report. It reads the
portable expression DAG julia/formula_inputs.json, containing the formula,
lower-order source, explicit integral primitives and conic evaluations, and
both sides of coefficient-table identities. This is shared mathematical
input, exported from the Python formula definitions, not an independent
transcription from the typeset paper. Differentiation, atom collection,
Clausen duplication and rational reduction are implemented in Julia.
Nemo uses FLINT, so the two implementations share an arithmetic backend.
Neither runner replaces the analytic regularity/uniqueness argument.

Maintenance only: julia/export_inputs.py regenerates the portable inputs
using the Python formula definitions after verify_all.py has generated the
panel primitives. End users need neither Python nor generated pickle files
to run the Julia certificate. The JSON contains no differentiated residuals.
The pinned Julia environment and certificate hashes record the tested inputs.

Terminal presentation: the Python runner uses Rich; the Julia runner uses
native styled text without an extra display dependency. Redirected logs are
plain text. Set NO_COLOR=1 to disable color, or FOCK_CERT_VERBOSE=1 to show
every individual Julia identity. Failures exit nonzero and retain diagnostics.
