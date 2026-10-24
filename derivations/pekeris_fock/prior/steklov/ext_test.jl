include("steklov.jl")
d = 8
model = channel_model(d)
setup = interior_setup(model, E_HE, 1.0, d)
for (lift, sets, omega, tol) in [
        (:gauge, (:outer, :compact), 12, 1e-14), (:gauge, (:outer, :compact), 12, 1e-12),
        (:gauge, (:outer, :compact), 12, 1e-16),
        (:poly, (:outer, :compact), 10, 1e-14), (:poly, (:outer, :compact), 12, 1e-14),
        (:poly, (:outer, :gauge), 10, 1e-14), (:poly, (:outer, :gauge), 12, 1e-14),
        (:poly, (:gauge,), 12, 1e-14), (:poly, (:outer,), 12, 1e-14)]
    t = @elapsed r = run_match(model, E_HE, 1.0, d; omega, setup, lift, sets, tol)
    pc = physical_constants(r.setup, r.a)
    @printf("%-6s %-20s om=%2d tol=%.0e: lambda=%+.3e a21=%.12f kept=%d/%d (%.1fs)\n",
            lift, string(sets), omega, tol, r.lambda, pc.a21, r.ext.kept, r.ext.total, t)
end
