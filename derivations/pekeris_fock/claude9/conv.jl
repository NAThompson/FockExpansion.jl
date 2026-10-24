include("psi41z2_big.jl")
setprecision(BigFloat, 256) do
    for (α, θ) in ((0.7, 1.1), (1.55, 0.05), (0.02, 3.1))
        a, t = BigFloat(α), BigFloat(θ)
        v = [psi41z2_raw(a, t; h = BigFloat(1) / k) for k in (32, 64, 128)]
        @time psi41z2_raw(a, t; h = BigFloat(1) / 64)
        println((α, θ), "  ", Float64(v[1] - v[2]), "  ", Float64(v[2] - v[3]))
    end
end
