include("steklov.jl")
include("exterior_dd.jl")
L = 4
model = channel_model(8)        # interior/trace harmonics up to degree 8 available for channels
setup = interior_setup(model, E_HE, 1.0, L)
k = sqrt(-2(E_HE + 2))
H1 = [(2.0, k, 10), (4.0, 4.0, 6), (8.0, 8.0, 6)]
H2 = [(2.0, k, 12), (4.0, 4.0, 8), (8.0, 8.0, 6)]
H3 = [(2.0, k, 10), (1.2, 1.2, 6), (4.0, 4.0, 6), (8.0, 8.0, 6), (16.0, 16.0, 4)]
cases = [("H1 c4 P4", H1, 4, 4), ("H1 c4 P8", H1, 8, 4), ("H1 c6 P6", H1, 6, 6),
         ("H1 c8 P6", H1, 6, 8), ("H2 c4 P4", H2, 4, 4), ("H3 c4 P4", H3, 4, 4),
         ("H3 c6 P8", H3, 8, 6)]
for (name, zsets, nchan, lchan) in cases
    t = @elapsed r = run_match_dd(model, E_HE, 1.0, L; zsets, nchan, lchan, setup, verbose = true)
    pc = physical_constants(r.setup, r.a)
    @printf("%-10s n=%4d: lambda=%+.6e a21=%.12f a3=%.10f a4=%.10f (%.1fs)\n", name, r.nterms,
            r.lambda, pc.a21, r.a[3], r.a[4], t)
    flush(stdout)
end
