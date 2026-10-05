# Fast azimuthal kernel ∫₀^{2π} G dφ for the k = 4 Green's function. With
# 1 - cos γ = δ + 2d sin²(φ/2), the only difficulty is the peak at φ = 0 when δ ≪ d
# (nearby points). On φ ∈ [0, π/2] use sin(φ/2) = σ sinh v, σ = √(δ/(2d)), so that
# 1 - cos γ = δ cosh²v and the integrand is smooth in v; plain Gauss elsewhere.
using FockExpansion, QuadGK, Printf
const F=FockExpansion
function kernel_fast(δ, d, n)
    x, w=gauss(n)
    G=F.green_k4
    s=0.0
    if δ<d
        σ=sqrt(δ/(2d)); vmax=asinh(sin(π/4)/σ)
        for (xi, wi) in zip(x, w)
            v=vmax*(xi+1)/2; sh=sinh(v)
            sφ=σ*sh; cφ=sqrt(1-sφ^2)                # sin(φ/2), cos(φ/2)
            s+=wi*vmax/2*G(δ*cosh(v)^2)*2σ*cosh(v)/cφ
        end
        lo=π/2
    else
        lo=0.0
    end
    for (xi, wi) in zip(x, w)
        φ=lo+(π-lo)*(xi+1)/2
        s+=wi*(π-lo)/2*G(δ+2d*sin(φ/2)^2)
    end
    2s
end
ref(δ, d)=2quadgk(φ->F.green_k4(δ+2d*sin(φ/2)^2), 0, 1e-8, 1e-6, 1e-4, 1e-2, 0.1, 1, π; rtol=1e-14)[1]
for (δ, d) in ((0.3, 0.2), (1e-2, 0.5), (1e-5, 0.5), (1e-10, 0.3), (1e-14, 0.8), (0.5, 0.0))
    r=ref(δ, d)
    @printf("δ=%-6g d=%-4g ref % .14f  n=12 %.1e  n=20 %.1e  n=32 %.1e\n", δ, d, r, kernel_fast(δ,d,12)-r, kernel_fast(δ,d,20)-r, kernel_fast(δ,d,32)-r)
end
