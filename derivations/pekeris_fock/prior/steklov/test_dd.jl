include("steklov.jl")
include("exterior_dd.jl")
# Cross-check analytic-r12 block against 3D quadrature in Float64.
terms = hylleraas_terms(4, 2.0, 1.3444)
Rs = 1.0
r2 = exterior_rule_2d(Rs, Float64)
A2 = correction_block(terms, r2, Rs, E_HE, 2.0, 1.0, Float64)
r3 = exterior_rule(Rs; nR = 60, nphi = 40, nu = 40)
A3 = exterior_matrix([], terms, r3, Rs, E_HE; Z = 2.0, g = 1.0)
println("corr block analytic vs 3D: ", norm(A2 - A3, Inf) / norm(A3, Inf))
# Double64 vs Float64 analytic
r2d = exterior_rule_2d(Rs, Double64)
A2d = correction_block(terms, r2d, Rs, E_HE, 2.0, 1.0, Double64)
println("Double64 vs Float64: ", Float64(norm(A2d - A2, Inf) / norm(A2d, Inf)))
# validation (i) in DD
model = channel_model(6; Z = 2.0, g = 0.0)
t = @elapsed r = run_match_dd(model, -4.0, 1.0, 6; omega = 6, kappa = 2.0, verbose = true)
@printf("val(i) DD gauge lift: lambda=%.3e max|a|=%.2e (%.1fs)\n", r.lambda, maximum(abs, r.a[2:end]), t)
# helium d=8
model = channel_model(8)
setup = interior_setup(model, E_HE, 1.0, 8)
for omega in (8, 10, 12)
    t = @elapsed r = run_match_dd(model, E_HE, 1.0, 8; omega, setup, verbose = true)
    pc = physical_constants(r.setup, r.a)
    @printf("He d=L=8 om=%d n=%d: lambda=%+.4e a21=%.12f %s %d/%d (%.1fs)\n", omega, r.nterms,
            r.lambda, pc.a21, r.ext.method, r.ext.kept, r.ext.total, t)
end
