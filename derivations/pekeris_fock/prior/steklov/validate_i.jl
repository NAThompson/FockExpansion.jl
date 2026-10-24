include("steklov.jl")
# Validation (i): no e-e repulsion, Z = 2, E = -4, Psi = exp(-2(r1+r2)), Phi = 1.
for (d, L, Rs, omega, lift) in [(2, 2, 1.0, 6, :frozen), (4, 4, 1.0, 8, :frozen),
                                (6, 6, 1.0, 10, :frozen), (6, 6, 0.6, 10, :frozen),
                                (6, 4, 1.0, 10, :frozen), (6, 6, 1.0, 6, :gauge)]
    model = channel_model(d; Z = 2.0, g = 0.0)
    t = @elapsed r = run_match(model, -4.0, Rs, L; omega, lift, sets = (:outer,), kappa = 2.0)
    pc = physical_constants(r.setup, r.a)
    @printf("d=%d L=%d R*=%.2f omega=%d lift=%s: lambda=%.3e lambda2=%.3e max|a[2:end]|=%.2e a21=%.2e asym=%.1e ode=%.1e kept=%d/%d (%.1fs)\n",
            d, L, Rs, omega, lift, r.lambda, r.lambda2, maximum(abs, r.a[2:end]), pc.a21,
            r.setup.asym, r.setup.ode, r.ext.kept, r.ext.total, t)
end
