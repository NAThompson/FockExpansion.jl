include("steklov.jl")
include("exterior_dd.jl")
L = 4
model = channel_model(8)
setup = interior_setup(model, E_HE, 1.0, L)
k = sqrt(-2(E_HE + 2))
zsets = [(2.0, k, 10), (4.0, 4.0, 6), (8.0, 8.0, 6)]
for fast in (true, false)
    t = @elapsed r = run_match_dd(model, E_HE, 1.0, L; zsets, nchan = 4, lchan = 4, setup, verbose = true, fast)
    pc = physical_constants(r.setup, r.a)
    @printf("fast=%s lambda=%+.10e a21=%.14f (%.1fs)\n", fast, r.lambda, pc.a21, t)
end
