include("steklov.jl")
include("exterior_dd.jl")
L = 4
model = channel_model(8)
setup = interior_setup(model, E_HE, 1.0, L)
k = sqrt(-2(E_HE + 2))
zsets = [(2.0, k, 10), (4.0, 4.0, 6), (8.0, 8.0, 6)]
for (nR, nphi, nu, beta) in [(60, 40, 30, 2.5), (80, 48, 36, 2.5), (60, 40, 30, 3.5), (40, 30, 24, 2.5), (40, 24, 20, 2.5), (50, 32, 24, 3.0)]
    t = @elapsed r = run_match_dd(model, E_HE, 1.0, L; zsets, nchan = 4, lchan = 4, setup, nR, nphi, nu, beta)
    pc = physical_constants(r.setup, r.a)
    @printf("nR=%d nphi=%d nu=%d beta=%.1f: lambda=%+.10e a21=%.14f (%.1fs)\n", nR, nphi, nu, beta, r.lambda, pc.a21, t)
    flush(stdout)
end
