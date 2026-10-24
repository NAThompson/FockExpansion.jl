# Validation (i): g = 0, Z = 2, E = -4, Psi = exp(-2(r1+r2)), Phi = 1.
# Exact: every cusp seed beyond a_1 is 0, a21 = 0, and the physical order-4 projections are
# <(-2 sigma)^4/24, Y> (sigma = cos phi + sin phi), computed by quadrature below.
include("production.jl")
Z = 2.0
E = -4.0
kap = 2.0
rule = sphere_rule(80, 80)
sig = rule.c .+ rule.s
for (lift, zsets, nchan, L, d) in [
        (:gauge, [(Z, kap, 6)], 0, 4, 8),
        (:frozen, [(Z, kap, 8), (4.0, 4.0, 4)], 2, 4, 8),
        (:frozen, [(Z, kap, 10), (4.0, 4.0, 6), (8.0, 8.0, 6)], 4, 4, 8),
        (:frozen, [(Z, kap, 12), (4.0, 4.0, 8), (8.0, 8.0, 6)], 6, 4, 8),
        (:frozen, [(Z, kap, 12), (4.0, 4.0, 8), (8.0, 8.0, 6)], 6, 6, 12)]
    model = channel_model(d; Z, g = 0.0)
    setup = interior_setup(model, E, 1.0, L)
    t = @elapsed ext = exterior_form_dd(model, E, 1.0, setup.nL; zsets, nchan, lchan = L,
                                        lift, gamma = 3.0)
    r = match(setup, Float64.(ext.S))
    pc = physical_constants(setup, r.a)
    i4 = findall(==(4), model.basis.K)
    A = harmonic_values(model.basis, sphere_q(rule), sphere_s(rule))
    exact40 = A[:, i4]' * (rule.w .* (2 .* sig) .^ 4 ./ 24)
    @printf("lift=%-6s L=%d d=%2d nb=%4d: lambda=%+.2e lambda2=%.2e max|a_2..|=%.2e a21=%+.2e psi40 err=%.2e (exact %s) (%.0fs)\n",
            lift, L, d, ext.nb, r.lambda, r.lambda2, maximum(abs, r.a[2:end]), pc.a21,
            maximum(abs, pc.psi40 .- exact40), string(round.(exact40; digits = 10)), t)
    flush(stdout)
end
