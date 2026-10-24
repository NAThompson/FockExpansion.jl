include("psi40_big.jl")
setprecision(BigFloat, 256)
E0 = big"-2.903724377034119598311"; a21 = big"0.47674788"
g = SGrid(BigFloat, 80); ψ = sfock(g; Z = 2, E = E0, kmax = 4, free = Dict((2, 1) => a21))
errs = Tuple{Float64,Float64,Float64}[]
for a in range(0.05, π - 0.05, length = 15), t in range(0.05, π - 0.05, length = 15)
    r = FockExpansion.psi40(a, t; Z = 2.0, E = Float64(E0), a21 = Float64(a21))
    push!(errs, (abs(r - Float64(interp(g, ψ[(4, 0)], big(a), big(t)))), a, t))
end
sort!(errs, rev = true)
println("max errors: ", errs[1:6]); println("median: ", errs[end÷2][1])
