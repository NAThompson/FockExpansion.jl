include("outer_series.jl")
E=big"-2.903724377034119598311159245194404446696905"
o=outer(E, 40)
κ2=sqrt(2(-big(1)/2-E))
for r in (big"0.25", big"1.0", big"3.0")
    s=[evalg(o.g[(k, 0)], r) for k in 0:40]
    @printf("r₁=%.2f signs of g_k0 for k=20..40: %s\n", Float64(r), join([x>0 ? "+" : "-" for x in s[21:41]]))
    # fit growth g_k ~ C Γ(k+β)/A^k from three consecutive even orders
    k=36; q1=s[k+3]/s[k+1]; q2=s[k+1]/s[k-1]
    @printf("     g_{k+2}/g_k at k=36: %.4f, k=38: %.4f   (κ₂-κ₁)⁻² k² ≈ %.1f\n", Float64(q2), Float64(q1), Float64((k+1)^2/(κ2-o.κ)^2))
end
