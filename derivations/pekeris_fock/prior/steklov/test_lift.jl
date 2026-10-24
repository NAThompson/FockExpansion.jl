include("steklov.jl")
include("exterior_dd.jl")
model = channel_model(8)
setup = interior_setup(model, E_HE, 1.0, 8)
for lift in (:decay, :decayK), omega in (8, 10, 12)
    t = @elapsed r = run_match_dd(model, E_HE, 1.0, 8; omega, setup, lift)
    pc = physical_constants(r.setup, r.a)
    @printf("lift=%s om=%d: lambda=%+.4e a21=%.12f %s %d/%d (%.1fs)\n", lift, omega,
            r.lambda, pc.a21, r.ext.method, r.ext.kept, r.ext.total, t)
end
