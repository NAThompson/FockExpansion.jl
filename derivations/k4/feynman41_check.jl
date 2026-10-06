include("feynman41.jl")
using FockExpansion
const F=FockExpansion
# zonal solutions satisfy the ODE? check (Λ²-32)U = f on a test zonal function by FD in γ
for kind in (:inv, :log, :xlog), ρ in (0.8,)
    f(c)=kind==:inv ? 1/(1-ρ*c) : kind==:log ? log(1-ρ*c) : (1-ρ*c)*log(1-ρ*c)
    γ=1.1; h=1e-3
    U(g)=zonal(kind, ρ, 1-ρ^2, g, 1-ρ*cos(g))
    lap=(U(γ+h)-2U(γ)+U(γ-h))/h^2+2cot(γ)*(U(γ+h)-U(γ-h))/(2h)
    # (Λ²-32)U = -4(Δ+8)U should equal f - κ'U₂ (harmonic projection); compare up to that:
    lhs=-4*(lap+8U(γ))
    # expect (Λ²-32)U = f - 4κU₂; compare the ratio of the difference to U₂ at two angles
    γ2=0.6; lhs2=-4*((U(γ2+h)-2U(γ2)+U(γ2-h))/h^2+2cot(γ2)*(U(γ2+h)-U(γ2-h))/(2h)+8U(γ2))
    U2(g)=4cos(g)^2-1
    @printf("%-5s ρ=%.1f: (lhs-f)/U₂ at γ=1.1: %.7f, at γ=0.6: %.7f (must agree)\n", kind, ρ, (lhs-f(cos(γ)))/U2(γ), (lhs2-f(cos(γ2)))/U2(γ2))
end
pts=[(0.7,1.1),(1.5,0.3),(1.2,2.0),(0.3,0.5),(2.6,2.5)]
d=Float64[]; Ys=Vector{Float64}[]
for (a,t) in pts
    v=psi41_z2_feynman(a,t); g=F.psi41_z2(a,t)
    push!(d, v-g); push!(Ys, [F.harmonic40(a,t), F.harmonic42(a,t)])
    @printf("(%.1f,%.1f) feynman %.12f green %.12f diff %.6e\n", a, t, v, g, v-g)
end
M=reduce(vcat, permutedims.(Ys)); b=M[1:2,:]\d[1:2]
println("harmonic fit from 2 points: ", b, "  residual at others: ", d[3:end]-M[3:end,:]*b)
