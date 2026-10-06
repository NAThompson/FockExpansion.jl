# Pure-Z harmonic coefficients of ψ₄₁ at high precision:
#   c_l(Z) = ⟨Y_l, Vψ₃₀ - ψ₄₂⟩/(6⟨Y_l,Y_l⟩)  at E = a₂₁ = 0,  = c₁Z + c₂Z² + c₃Z³.
# Fixed Gauss–Legendre tensor rules, graded toward α = 0, α = π/2 and θ = 0 (the Coulomb
# singularities; α ∈ (0, π/2] suffices by exchange symmetry), at two resolutions.
using FockExpansion, QuadGK, Printf, Base.Threads
setprecision(BigFloat, 128)
const TOL=big(1e-28)
Y40(a,t)=4cos(a)^2-1; Y42(a,t)=sin(a)^2*(3cos(t)^2-1)/2
function graded(a, b, toa, tob; levels, q=big(0.3))
    pts=[a, b]
    toa && append!(pts, [a+(b-a)/2*q^k for k in 0:levels])
    tob && append!(pts, [b-(b-a)/2*q^k for k in 0:levels])
    sort!(unique!(pts))
end
function rule(pts, n)
    x, w=gauss(BigFloat, n)
    xs=BigFloat[]; ws=BigFloat[]
    for i in 1:length(pts)-1
        h=(pts[i+1]-pts[i])/2; m=(pts[i+1]+pts[i])/2
        append!(xs, m .+ h .* x); append!(ws, h .* w)
    end
    xs, ws
end
function proj(Z, n, levels)
    Z=big(Z); P=big(π)
    A=rule(graded(big(0), P/2, true, true; levels), n)
    T=rule(graded(big(0), P, true, false; levels), n)
    acc=zeros(BigFloat, 4)
    for (a, wa) in zip(A...)
        for (t, wt) in zip(T...)
            ξ=sqrt(2sin(P/4-a/2)^2+2sin(a)*sin(t/2)^2)
            V=1/ξ-Z*(1/sin(a/2)+1/cos(a/2))
            p42=(P-2)*(5P-14)/(180P^2)*Z^2*(1-2sin(a)^2*sin(t)^2)
            g=V*psi30(a, t; Z, E=big(0), a21=big(0), rtol=TOL)-p42
            w=wa*wt*sin(a)^2*sin(t)
            y0, y2=Y40(a, t), Y42(a, t)
            acc.+=w .* [y0*g, y2*g, y0^2, y2^2]
        end
    end
    acc[1]/(6acc[3]), acc[2]/(6acc[4])
end
n, levels=parse(Int, ARGS[1]), parse(Int, ARGS[2])
res=Vector{Any}(undef, 3)
@threads for i in 1:3
    res[i]=proj(i, n, levels)
end
A=BigFloat[1 1 1; 2 4 8; 3 9 27]
for (l, k) in ((0, 1), (2, 2))
    c=A\[res[i][k] for i in 1:3]
    println("n=$n levels=$levels l=$l: ", join([string(round(x, sigdigits=28)) for x in c], "  "))
end
