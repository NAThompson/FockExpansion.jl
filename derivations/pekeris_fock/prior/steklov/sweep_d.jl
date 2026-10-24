# Interior-degree sweep with cached exterior forms.
include("production.jl")
Rs = parse(Float64, ARGS[1])
name = ARGS[2]
Ls = parse.(Int, split(ARGS[3], ","))
ds = parse.(Int, split(ARGS[4], ","))
for L in Ls
    ext = get_ext(L, Rs, E_HE, name)
    for d in ds
        d < L && continue
        t = @elapsed p = solve_point(d, L, Rs, E_HE, ext)
        @printf("R*=%.2f L=%2d %s d=%2d: lambda=%+.6e a21=%.13f a3=%+.12f a4=%+.12f asym=%.1e ode=%.1e (%.1fs)\n",
                Rs, L, name, d, p.lambda, p.a21, p.a[3], p.a[4], p.asym, p.ode, t)
        flush(stdout)
    end
end
