# Exterior basis convergence at fixed interior: for each L and basis, a21 at d_int = 24.
include("production.jl")
Rs = parse(Float64, ARGS[1])
L = parse(Int, ARGS[2])
for name in ARGS[3:end]
    ext = get_ext(L, Rs, E_HE, name)
    for d in (16, 24)
        p = solve_point(d, L, Rs, E_HE, ext)
        @printf("R*=%.2f L=%2d %-7s nb=%4d d=%2d: lambda=%+.6e a21=%.13f a3=%+.11f a4=%+.11f\n",
                Rs, L, name, ext.nb, d, p.lambda, p.a21, p.a[3], p.a[4])
    end
    flush(stdout)
end
