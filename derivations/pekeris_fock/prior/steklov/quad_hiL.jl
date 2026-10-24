# Angular quadrature refinement for the production exterior forms (Codex review 2, point 1).
include("production.jl")
Rs = parse(Float64, ARGS[1])
name = ARGS[2]
L = parse(Int, ARGS[3])
quad = (parse(Int, ARGS[4]), parse(Int, ARGS[5]))
ds = parse.(Int, split(ARGS[6], ","))
ext = get_ext(L, Rs, E_HE, name; quad)
for d in ds
    p = solve_point(d, L, Rs, E_HE, ext)
    @printf("R*=%.2f L=%2d %s q=%dx%d kept=%d/%d d=%2d: lambda=%+.6e a21=%.13f a3=%+.12f a4=%+.12f a5=%+.12f\n",
            Rs, L, name, quad..., ext.kept, ext.total, d, p.lambda, p.a21, p.a[3], p.a[4], p.a[5])
    flush(stdout)
end
