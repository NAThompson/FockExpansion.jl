# The truncated outer series is a formal solution: (H - E)ψ_K / ψ_K = O(R^{-K-1}).
# For L = 0, ∇₁² + ∇₂² = ∂²₁ + (2/r₁)∂₁ + ∂²₂ + (2/r₂)∂₂ + (1/r₁² + 1/r₂²)(∂²_θ + cot θ ∂_θ),
# applied by central differences in 512-bit arithmetic.
include("outer_series.jl")
E=big"-2.903724377034119598311159245194404446696905"
o=outer(E, 24)
P(l, x)=l==0 ? one(x) : l==1 ? x : ((2l-1)*x*P(l-1, x)-(l-1)*P(l-2, x))/l
ψout(r1, R, θ, K)=exp(-o.κ*R)*R^(o.σ-1)*sum(R^(-k)*evalg(o.g[(k, l)], r1)*P(l, cos(θ)) for k in 0:K for l in 0:k if haskey(o.g, (k, l)))
function residual(r1, R, θ, K; h=big(10)^-30)
    f(a, b, c)=ψout(a, b, c, K)
    f0=f(r1, R, θ)
    d2(g, x)=(g(x+h)-2g(x)+g(x-h))/h^2
    d1(g, x)=(g(x+h)-g(x-h))/(2h)
    lap=d2(x->f(x, R, θ), r1)+2/r1*d1(x->f(x, R, θ), r1)+d2(x->f(r1, x, θ), R)+2/R*d1(x->f(r1, x, θ), R)+
        (1/r1^2+1/R^2)*(d2(x->f(r1, R, x), θ)+cot(θ)*d1(x->f(r1, R, x), θ))
    r12=sqrt(r1^2+R^2-2r1*R*cos(θ))
    (-lap/2+(-Z/r1-Z/R+1/r12-E)*f0)/f0
end
r1, θ=big"0.7", big"1.1"
@printf("relative residual (H-E)ψ_K/ψ_K at r₁=0.7, θ=1.1\n")
for R in (big"10.0", big"20.0", big"40.0")
    @printf("R=%2d:", Int(R))
    for K in (0, 2, 4, 8, 12); @printf("  K=%-2d %.1e", K, Float64(abs(residual(r1, R, θ, K)))); end
    println()
end
