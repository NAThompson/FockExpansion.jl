# Does adding R⁴(ψ₄₂ log²R + ψ₄₁ log R) to Ψ⁽³⁾ reduce (H-E)Ψ?  ψ₄₀ is set to 0.
# The residual of a truncated Fock series is exact and finite (see README of this dir):
#   Ψ⁽³⁾:        R²[s₀ + s₁L] - E R³(ψ₃₀ + ψ₃₁L)
#   +ψ₄₂:        R²[(s₀-ψ₄₂) + (s₁-12ψ₄₂)L] + R³[-Eψ₃₀ - Eψ₃₁L + Vψ₄₂L²] - E R⁴ψ₄₂L²
#   +ψ₄₂,ψ₄₁:    R²[s₀-6ψ₄₁-ψ₄₂] + R³[-Eψ₃₀ + (Vψ₄₁-Eψ₃₁)L + Vψ₄₂L²] - E R⁴(ψ₄₁L + ψ₄₂L²)
#   full k=4:    R³[(Vψ₄₀-Eψ₃₀) + (Vψ₄₁-Eψ₃₁)L + Vψ₄₂L²] - E R⁴(ψ₄₀ + ψ₄₁L + ψ₄₂L²)
# with L = log R, s₀ = Vψ₃₀ - Eψ₂₀, s₁ = Vψ₃₁ - Eψ₂₁. ψ₄₀ uses a₄₀ = a₄₂ = 0, which
# affects only the R³ and R⁴ terms.
using FockExpansion, Printf, Base.Threads
const F=FockExpansion
const Z, E, a21 = 2.0, -2.9037243770341195983, 0.47674787900
function node(α, θ; n=8, levels=8)
    ξ=F.xi_stable(α, θ); V=1/ξ-Z*(1/sin(α/2)+1/cos(α/2))
    p30=psi30(α, θ; Z, E, a21, rtol=1e-12); p31=psi31(α, θ; Z)
    s0=V*p30-E*psi20(α, θ; Z, E, a21); s1=V*p31-E*psi21(α, θ; Z)
    p42=psi42(α, θ; Z); p41=psi41(α, θ; Z, E, a21)
    p40=psi40(α, θ; Z, E, a21, n, levels)
    (; V, p30, p31, s0, s1, p42, p41, p40)
end
res3(c, R)=(L=log(R); R^2*(c.s0+c.s1*L)-E*R^3*(c.p30+c.p31*L))
res42(c, R)=(L=log(R); R^2*((c.s0-c.p42)+(c.s1-12c.p42)*L)+R^3*(-E*c.p30-E*c.p31*L+c.V*c.p42*L^2)-E*R^4*c.p42*L^2)
res4(c, R)=(L=log(R); R^2*(c.s0-6c.p41-c.p42)+R^3*(-E*c.p30+(c.V*c.p41-E*c.p31)*L+c.V*c.p42*L^2)-E*R^4*(c.p41*L+c.p42*L^2))
res40(c, R)=(L=log(R); R^3*((c.V*c.p40-E*c.p30)+(c.V*c.p41-E*c.p31)*L+c.V*c.p42*L^2)-E*R^4*(c.p40+c.p41*L+c.p42*L^2))
const Rs=(1e-6, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2, 0.3, 0.5)

println("Pointwise |(H-E)Ψ|/R² at the angles of the ψ₃₀ paper (Z=2, ground-state E, a₂₁):")
for (α, θ) in ((1.1, 1.2), (0.4, 0.7), (2.2, 2.4))
    c=node(α, θ; n=10, levels=10)
    @printf("\n(α,θ)=(%.1f,%.1f): s₀=%.4f s₁=%.4f  ψ₄₂=%.5f ψ₄₁=%.5f  s₀'=s₀-6ψ₄₁-ψ₄₂=%.4f\n", α, θ, c.s0, c.s1, c.p42, c.p41, c.s0-6c.p41-c.p42)
    @printf("%8s %14s %14s %14s %14s\n", "R", "Ψ⁽³⁾", "+ψ₄₂", "+ψ₄₂+ψ₄₁", "+ψ₄₀ (all k=4)")
    for R in Rs
        @printf("%8.0e %14.6f %14.6f %14.6f %14.6f\n", R, abs(res3(c,R))/R^2, abs(res42(c,R))/R^2, abs(res4(c,R))/R^2, abs(res40(c,R))/R^2)
    end
end

# L² norm over the angles (singlet: α ∈ (0, π/2] suffices), coarse graded grid.
rule=F.gauss(6)
A=F.graded_points(0.0, π/2, true, true; levels=4)
T=F.graded_points(0.0, Float64(π), true, false; levels=4)
nodes=[(m+h*x, w*h) for (lo,hi) in zip(A[1:end-1],A[2:end]) for (x,w) in zip(rule...) for (m,h) in (((lo+hi)/2,(hi-lo)/2),)]
tnodes=[(m+h*x, w*h) for (lo,hi) in zip(T[1:end-1],T[2:end]) for (x,w) in zip(rule...) for (m,h) in (((lo+hi)/2,(hi-lo)/2),)]
grid=[(a, t, wa*wt*sin(a)^2*sin(t)) for (a,wa) in nodes for (t,wt) in tnodes]
println("\nL² grid: ", length(grid), " nodes, ", nthreads(), " threads")
cs=Vector{Any}(undef, length(grid))
t=@elapsed @threads for i in eachindex(grid)
    a, θ, _=grid[i]; cs[i]=node(a, θ; n=6, levels=6)
end
@printf("grid evaluated in %.0f s\n", t)
nrm(f)=sqrt(sum(g[3]*f(c)^2 for (g,c) in zip(grid,cs)))
@printf("\n‖(H-E)Ψ‖/R² over the angles\n%8s %14s %14s %14s %14s\n", "R", "Ψ⁽³⁾", "+ψ₄₂", "+ψ₄₂+ψ₄₁", "+ψ₄₀ (all k=4)")
for R in Rs
    n3=nrm(c->res3(c,R))/R^2; n42=nrm(c->res42(c,R))/R^2; n4=nrm(c->res4(c,R))/R^2; n40=nrm(c->res40(c,R))/R^2
    @printf("%8.0e %14.6f %14.6f %14.6f %14.6f\n", R, n3, n42, n4, n40)
end
open(joinpath(@__DIR__, "residual_grid.csv"), "w") do io
    println(io, "alpha,theta,weight,V,psi30,psi31,s0,s1,psi42,psi41,psi40")
    for (g,c) in zip(grid,cs)
        println(io, join((g..., c.V, c.p30, c.p31, c.s0, c.s1, c.p42, c.p41), ","))
    end
end
