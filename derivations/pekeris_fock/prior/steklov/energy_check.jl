# Smallest matching eigenvalue at E0 and E0 +- 1e-6 (L = 10, B10P6, R* = 1, d = 32).
include("production.jl")
L, Rs, name, d = 10, 1.0, "B10P6", 32
lam = Float64[]
for dE in (-1e-6, 0.0, 1e-6)
    E = E_HE + dE
    ext = get_ext(L, Rs, E, name)
    p = solve_point(d, L, Rs, E, ext)
    push!(lam, p.lambda)
    @printf("E = E0 %+.0e: lambda_min = %+.6e  lambda_2 = %.4e  a21 = %.13f\n", dE, p.lambda, p.lambda2, p.a21)
end
slope = (lam[3] - lam[1]) / 2e-6
@printf("d lambda/dE = %.6e;  zero crossing at E0 %+.3e\n", slope, -lam[2] / slope)
