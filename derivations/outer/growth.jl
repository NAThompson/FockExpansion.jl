include("outer_series.jl")
E=big"-2.903724377034119598311159245194404446696905"
t=@elapsed o=outer(E, 40)
@printf("κ = %.12f, σ = %.12f  (%.1f s)\n", Float64(o.κ), Float64(o.σ), t)
for k in 0:2:40
    v=maximum(abs(evalg(o.g[(k, l)], big(1)/2)) for l in 0:k if haskey(o.g, (k, l)))
    @printf("k=%2d  max_l |g_kl(r₁=0.5)| = %.3e   k!/(κ₂-κ₁)^k = %.3e\n", k, Float64(v), Float64(factorial(big(k))/(sqrt(2(-big(1)/2-E))-o.κ)^k))
end
