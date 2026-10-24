using FockExpansion, ForwardDiff, Printf
E0 = big"-2.903724377034119598311159245194404446696925310"; A21 = big"0.47674787900"
f(z, rt) = psi30(z[1], z[2]; Z = big(2), E = E0, a21 = A21, rtol = rt)
for (bits, rt) in ((128, big"1e-25"), (192, big"1e-35"), (256, big"1e-45"))
    setprecision(BigFloat, bits) do
        z = [big"0.8", big"1.2"]
        f(z, rt); ForwardDiff.gradient(w -> f(w, rt), z)  # compile
        t = @elapsed v = f(z, BigFloat(rt))
        t2 = @elapsed g = ForwardDiff.gradient(w -> f(w, BigFloat(rt)), z)
        vf = psi30(0.8, 1.2; Z = 2.0, E = Float64(E0), a21 = Float64(A21), rtol = 1e-13)
        @printf("%d bits: value %.25e (%.3fs)  grad (%.3fs)  Float64 diff %.2e\n", bits, v, t, t2, Float64(v) - vf)
        println("   grad ", g)
    end
end
