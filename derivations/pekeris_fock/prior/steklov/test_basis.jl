include("steklov.jl")
include("exterior_dd.jl")
L = 4
model = channel_model(L)
setup = interior_setup(model, E_HE, 1.0, L)
k = sqrt(-2(E_HE + 2))
cases = [
    ("outer8", [(2.0, k, 8)], 0),
    ("outer8+chan3", [(2.0, k, 8)], 3),
    ("outer8+chan6", [(2.0, k, 8)], 6),
    ("outer8+z4", [(2.0, k, 8), (4.0, 4.0, 6)], 0),
    ("outer8+z4+z8", [(2.0, k, 8), (4.0, 4.0, 6), (8.0, 8.0, 6)], 0),
    ("outer10+chan6", [(2.0, k, 10)], 6),
    ("outer10+z4+z8+chan4", [(2.0, k, 10), (4.0, 4.0, 6), (8.0, 8.0, 6)], 4),
]
for (name, zsets, nchan) in cases
    t = @elapsed r = run_match_dd(model, E_HE, 1.0, L; zsets, nchan, setup, verbose = true)
    pc = physical_constants(r.setup, r.a)
    @printf("%-22s n=%4d: lambda=%+.4e a21=%.12f (%.1fs)\n", name, r.nterms, r.lambda, pc.a21, t)
    flush(stdout)
end
