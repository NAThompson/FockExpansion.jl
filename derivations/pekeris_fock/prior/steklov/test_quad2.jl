include("steklov.jl")
include("exterior_dd.jl")
L = 4
model = channel_model(8)
setup = interior_setup(model, E_HE, 1.0, L)
k = sqrt(-2(E_HE + 2))
zsets = [(2.0, k, 10), (4.0, 4.0, 6), (8.0, 8.0, 6)]
for (nR, edges, npanel, beta) in [(40, [0.25, 0.75, 2.0, 5.0], 12, 2.5), (40, [0.25, 0.75, 2.0, 5.0], 16, 2.5),
                                  (50, [0.25, 0.75, 2.0, 5.0], 20, 2.5), (40, [0.15, 0.4, 1.0, 2.5, 6.0], 16, 2.5),
                                  (60, [0.15, 0.4, 1.0, 2.5, 6.0], 20, 2.0)]
    t = @elapsed r = run_match_dd(model, E_HE, 1.0, L; zsets, nchan = 4, lchan = 4, setup, nR, nphi = 24, nu = 20, beta, edges, npanel)
    pc = physical_constants(r.setup, r.a)
    @printf("nR=%d edges=%s np=%d: lambda=%+.10e a21=%.14f (%.1fs)\n", nR, string(edges), npanel, r.lambda, pc.a21, t)
    flush(stdout)
end
