# Split the 2s x 2s interior error: interior channel degree d versus trace degree L.
include("steklov.jl")
Z = 2.0
E = -Z^2 / 4
phi_exact(r1, r2) = (1 - Z * r1 / 2) * exp(Z * r1 / 2) * (1 - Z * r2 / 2) * exp(Z * r2 / 2)
Rs = 1.0
println("d,L,err_a20")
for d in (8, 12, 16, 20, 24, 28)
    model = channel_model(d; Z, g = 0.0)
    ser = fock_series(model, E, max(44, 2d + 28), [Rs])
    rule = sphere_rule(2d + 40, 2d + 40)
    A = harmonic_values(model.basis, sphere_q(rule), sphere_s(rule))
    y = A' * (rule.w .* phi_exact.(Rs .* rule.c, Rs .* rule.s))
    deg = harmonic_degree(model.basis)
    i20 = findfirst(==((2, 0)), model.basis.labels)
    for L in 4:4:d
        yL = copy(y)
        yL[deg .> L] .= 0
        a = ser.F[1] \ yL
        a ./= a[1]
        @printf("%d,%d,%+.3e\n", d, L, a[i20] + 1 / 32)
    end
end
