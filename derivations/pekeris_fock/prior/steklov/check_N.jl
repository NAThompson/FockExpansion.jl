# Interior series-order check: N = 2d + 28 (default) versus larger N, at several R*.
include("production.jl")
for Rs in (0.5, 0.75, 1.0, 1.5, 2.0), d in (16, 24)
    model = channel_model(d)
    vals = Float64[]
    for N in (2d + 28, 2d + 44, 2d + 64)
        ser = fock_series(model, E_HE, N, [Rs])
        push!(vals, norm(ser.F[1]) + norm(ser.Fp[1]))
        @printf("R*=%.2f d=%d N=%3d ode=%.1e |F|+|F'|=%.15e cond(F)=%.2e\n", Rs, d, N,
                ode_residual(model, E_HE, Rs, ser.F[1], ser.Fp[1], ser.Fpp[1]), vals[end], cond(ser.F[1]))
    end
end
