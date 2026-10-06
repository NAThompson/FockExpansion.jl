include("fock_spectral.jl")
using FockExpansion
setprecision(BigFloat, 256)
Z, E, a21 = big"2.0", big"-2.9", big"0.4"
for n in (32, 48, 64)
    t=@elapsed begin
        g=SGrid(BigFloat, n)
        ψ=sfock(g; Z, E, kmax=4, free=Dict((2,1)=>a21))
    end
    # compare with the BigFloat closed form of ψ₃₀ at three grid nodes
    errs=BigFloat[]
    for (i, j) in ((n÷4, n÷3), (n÷2, n÷2), (3n÷4, n÷5))
        α=g.α[i, j]; β=g.β[i, j]; θ=acos(cos(β)/sin(α))
        ref=psi30(α, θ; Z, E, a21, rtol=big(1e-40))
        push!(errs, abs(ψ[(3, 0)][i, j]-ref))
    end
    @printf("n=%d (256-bit, chain to k=4 in %.0f s): max |ψ₃₀ - closed form| = %.1e\n", n, t, Float64(maximum(errs)))
end
