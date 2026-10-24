include("steklov.jl")
d = 8
model = channel_model(d)
setup = interior_setup(model, E_HE, 1.0, d)
for sets in [(:outer,), (:gauge,)], omega in (8, 10, 12, 14), tol in (1e-12, 1e-14, 1e-16)
    t = @elapsed r = run_match(model, E_HE, 1.0, d; omega, setup, lift = :gauge, sets, tol)
    pc = physical_constants(r.setup, r.a)
    @printf("%-12s om=%2d tol=%.0e: lambda=%+.3e a21=%.12f kept=%d/%d lmin=%.1e (%.1fs)\n",
            string(sets), omega, tol, r.lambda, pc.a21, r.ext.kept, r.ext.total, r.ext.lmin, t)
end
