# Fraction of the helium ground-state probability inside the hypersphere ρ ≤ ρ₀ (Hylleraas).
# dτ ∝ ρ⁵ sin²α sin θ dρ dα dθ.
include("hylleraas.jl")
using Serialization
setprecision(BigFloat, 128)
h=deserialize(joinpath(tempdir(), "hylleraas_12.jls"))
h=(; E=h.E, c=BigFloat.(h.c), idx=h.idx, ζs=BigFloat.(h.ζs))
function gauss_legendre(m)
    x=zeros(m); w=zeros(m)
    for i in 1:m
        z=cos(π*(i-0.25)/(m+0.5))
        for _ in 1:100
            p0, p1=1.0, z
            for k in 2:m; p0, p1=p1, ((2k-1)*z*p1-(k-1)*p0)/k; end
            dp=m*(z*p1-p0)/(z^2-1); z-=p1/dp
        end
        p0, p1=1.0, z
        for k in 2:m; p0, p1=p1, ((2k-1)*z*p1-(k-1)*p0)/k; end
        x[i]=z; w[i]=2/((1-z^2)*(m*(z*p1-p0)/(z^2-1))^2)
    end
    x, w
end
x, w=gauss_legendre(24)
ang=[(π/2*(xa+1), π/2*(xt+1), wa*wt*π^2/4) for (xa, wa) in zip(x, w) for (xt, wt) in zip(x, w)]
shell(ρ)=sum(wt*sin(a)^2*sin(t)*Float64(hyl_eval(h, big(ρ)*cos(big(a)/2), big(ρ)*sin(big(a)/2), cos(big(t))))^2 for (a, t, wt) in ang)*ρ^5
# ρ ∈ [0, 12] in panels
edges=[0.0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 6, 9, 12]
cum=[0.0]
for i in 1:length(edges)-1
    a, b=edges[i], edges[i+1]
    push!(cum, cum[end]+sum(wi*(b-a)/2*shell((a+b)/2+(b-a)/2*xi) for (xi, wi) in zip(x, w)))
end
for (e, c) in zip(edges[2:end], cum[2:end])
    @printf("P(ρ ≤ %4.1f) = %.6f\n", e, c/cum[end])
end
