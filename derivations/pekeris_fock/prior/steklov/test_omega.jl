include("steklov.jl")
include("exterior_dd.jl")
L = parse(Int, ARGS[1])
model = channel_model(L)
setup = interior_setup(model, E_HE, 1.0, L)
for omega in parse.(Int, ARGS[2:end])
    t = @elapsed r = run_match_dd(model, E_HE, 1.0, L; omega, setup, verbose = true)
    pc = physical_constants(r.setup, r.a)
    @printf("L=%d om=%d n=%d: lambda=%+.4e a21=%.12f %s %d/%d (%.1fs)\n", L, omega, r.nterms,
            r.lambda, pc.a21, r.ext.method, r.ext.kept, r.ext.total, t)
    flush(stdout)
end
