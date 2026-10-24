# Controlled comparison (Codex review 1, point 8): same channel interior and exterior,
# Codex's constrained least-squares fit versus the eigenvector of the mismatch.
module CodexRef
include("/Users/nathompson7/Documents/Codex/2026-09-05/how-much-compu/outputs/boundary_matching/fixed_energy_matching.jl")
end
using Printf
E = CodexRef.E_REFERENCE
println("d,channels,radius,a21_leastsquares,a21_eigenvector,lambda_eig")
for d in (4, 8, 12, 16), rm in (0.75, 1.0)
    b = CodexRef.cuspmodel(d; nquad = 100)
    z = CodexRef.fixedfit(b, E; rm)
    w = CodexRef.cuspmismatch(b, E; rm, details = true)
    @printf("%d,%d,%.2f,%.12f,%.12f,%.3e\n", d, length(b.K), rm, z.a21, w.a21, w.match)
    flush(stdout)
end
