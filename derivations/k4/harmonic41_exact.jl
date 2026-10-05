# Pure-Z harmonic coefficients of ψ₄₁ to ~30 digits:
#   c_l(Z) = ⟨Y_l, Vψ₃₀ - ψ₄₂⟩/(6⟨Y_l,Y_l⟩)  at E = a₂₁ = 0,  = c₁Z + c₂Z² + c₃Z³.
# Nested adaptive BigFloat quadrature over α ∈ (0, π/2] (doubled by exchange symmetry), θ ∈ (0, π).
using FockExpansion, QuadGK, Printf, Base.Threads
setprecision(BigFloat, 128)
const TOL=big(1e-31)
Y40(a,t)=4cos(a)^2-1; Y42(a,t)=sin(a)^2*(3cos(t)^2-1)/2
function proj(Z)
    Z=big(Z)
    function f(a, t)
        ξ=sqrt(2sin(big(π)/4-a/2)^2+2sin(a)*sin(t/2)^2)
        V=1/ξ-Z*(1/sin(a/2)+1/cos(a/2))
        p42=(big(π)-2)*(5big(π)-14)/(180big(π)^2)*Z^2*(1-2sin(a)^2*sin(t)^2)
        g=V*psi30(a, t; Z, E=big(0), a21=big(0), rtol=TOL)-p42
        w=sin(a)^2*sin(t)
        [Y40(a,t)*g*w, Y42(a,t)*g*w, Y40(a,t)^2*w, Y42(a,t)^2*w]
    end
    inner(a)=quadgk(t->f(a, t), big(0), big(π)/2, big(π); rtol=TOL, atol=TOL)[1]
    v=quadgk(inner, big(0), big(π)/4, big(π)/2; rtol=TOL, atol=TOL)[1]
    v[1]/(6v[3]), v[2]/(6v[4])
end
res=Vector{Any}(undef, 3)
@threads for i in 1:3
    res[i]=proj(i)
    println("Z=$i done: ", Float64.(res[i])); flush(stdout)
end
# c(Z) = c₁Z + c₂Z² + c₃Z³ from Z = 1, 2, 3
A=BigFloat[1 1 1; 2 4 8; 3 9 27]
for (l, k) in ((0, 1), (2, 2))
    c=A\[res[i][k] for i in 1:3]
    println("l=$l: c₁ = ", round(c[1], sigdigits=30), "\n     c₂ = ", round(c[2], sigdigits=30), "\n     c₃ = ", round(c[3], sigdigits=30))
end
