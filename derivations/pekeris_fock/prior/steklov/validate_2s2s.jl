# Interior validation with a nonzero cusp-gauge seed (suggested by Codex review 1):
# g = 0, Z = 2, E = -Z^2/4, Psi = (1 - Z r1/2)(1 - Z r2/2) exp(-Z(r1+r2)/2).
# With S = -Z(r1+r2): Phi = prod (1 - x)exp(x), x = Z r/2, so
#   Phi_20 = -Z^2/8 (a_2 = 0),  Phi_40 = -Z^4 q^2/128,  a_(2,0) = -Z^4/512 = -1/32, a_(0,2) = 0.
# E = -1 lies in the continuum (threshold -2), so only the interior is tested here:
# project the exact trace onto the channels and recover a = F^{-1} y.
include("steklov.jl")
Z = 2.0
E = -Z^2 / 4
phi_exact(r1, r2) = (1 - Z * r1 / 2) * exp(Z * r1 / 2) * (1 - Z * r2 / 2) * exp(Z * r2 / 2)
for Rs in (0.5, 1.0), d in (4, 8, 12, 16, 20, 24)
    model = channel_model(d; Z, g = 0.0)
    ser = fock_series(model, E, max(44, 2d + 28), [Rs])
    rule = sphere_rule(2d + 40, 2d + 40)
    A = harmonic_values(model.basis, sphere_q(rule), sphere_s(rule))
    y = A' * (rule.w .* phi_exact.(Rs .* rule.c, Rs .* rule.s))
    a = ser.F[1] \ y
    a ./= a[1]
    i20 = findfirst(==((2, 0)), model.basis.labels)
    i02 = findfirst(==((0, 2)), model.basis.labels)
    @printf("R*=%.1f d=%2d m=%3d: a_2=%+.3e  a_(2,0)+1/32=%+.3e  a_(0,2)=%+.3e  cond(F)=%.1e\n",
            Rs, d, length(a), a[2], a[i20] + 1 / 32, a[i02], cond(ser.F[1]))
end
