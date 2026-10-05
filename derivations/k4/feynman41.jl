# ψ₄₁⁽²⁾ as a one-dimensional integral of elementary functions.
#
# Source (Λ²-32)ψ = B Σ± Q±(x)/(d± ξ), with x ∈ S³, u = x₁ = cos α, w = x₂ = sin α cos θ,
#   d± = √(2(1∓u)), ξ = √(1-w),  Q± = -1/3 - (7/3)w + (5/3)w² ± uw.
# Feynman: 1/(d ξ) = (1/(√2π)) ∫₀¹ dt/(√(t(1-t)) ℓ),  ℓ = 1 - p·x,  p± = (±t, 1-t, 0, 0).
# Polynomial factors via p-derivatives: x_i/ℓ = -∂ᵢ log ℓ,  x_i x_j/ℓ = ∂ᵢ∂ⱼ(ℓ log ℓ).
# Each of 1/ℓ, log ℓ, ℓ log ℓ is zonal about p̂; its regular solution is elementary
# (variation of parameters with W = U sin γ, W'' + 9W = -sin γ f/4 + κ sin 3γ).
using ForwardDiff, QuadGK, Printf, LinearAlgebra

# ∫₀^γ cos^m τ dτ, m = 0..6
function Cpow(γ, m)
    s, c=sincos(γ)
    m==0 && return γ
    m==1 && return s
    c^(m-1)*s/m+(m-1)/m*Cpow(γ, m-2)
end
# K_n = ∫₀^γ cos^n τ/(1-ρ cos τ) dτ via K_n = (K_{n-1} - C_{n-1})/ρ, K₀ = A(γ);
# ω = 1-ρ² is passed separately to avoid cancellation as ρ → 1.
function Ktab(γ, ρ, ω, N)
    A=γ==π ? oftype(γ, π)/sqrt(ω) : 2/sqrt(ω)*atan((1+ρ)/sqrt(ω)*tan(γ/2))
    K=[A]
    for n = 1:N
        push!(K, (K[end]-Cpow(γ, n-1))/ρ)
    end
    K  # K[n+1] = K_n
end
# J_n(c) = ∫ c^n/(1-ρc) dc (antiderivative) given lg = log(1-ρc); H_n = ∫ c^n log(1-ρc) dc
function Jtab(c, ρ, lg, N)
    J=[-lg/ρ]
    for n = 1:N
        push!(J, (J[end]-c^n/n)/ρ)
    end
    J
end
Hn(c, ρ, lg, n, J)=c^(n+1)*lg/(n+1)+ρ*J[n+2]/(n+1)

# U = W/sin γ for f ∈ (:inv, :log, :xlog), given ρ, ω = 1-ρ², γ and ℓ = 1-ρ cos γ.
function zonal(kind, ρ, ω, γ, ℓ)
    I1, I2=zonal_I(kind, ρ, ω, γ, ℓ)
    _, I2π=zonal_I(kind, ρ, ω, oftype(γ, π), 1+ρ; atpi = true)
    κ=-2I2π/oftype(γ, π)
    I1+=κ*sin(3γ)^2/6
    I2+=κ*(γ/2-sin(6γ)/12)
    (sin(3γ)*I1-cos(3γ)*I2)/(3sin(γ))
end
# Stable form: every coefficient that vanishes as ρ → 1 is written with explicit
# powers of ω = 1-ρ² (derivations/k4/zonal_coeffs.py), so nothing cancels.
#   I₁·(-4) = a1 log(1-ρ) + ac log ℓ + poly,    I₂·(-4) = (P log ℓ -) ρ^m [R A + trig],
#   Σₙ aₙ Kₙ = A·R - Σₙ aₙ Σ_{k<n} C_k ρ^{k-n},  R = Σ aₙ ρ⁻ⁿ (factored).
function trigsum(a, γ, ρ)
    s=zero(γ*ρ)
    for (n, an) in a, k = 0:(n-1)
        s-=an*Cpow(γ, k)*ρ^(k-n)
    end
    s
end
function zonal_I(kind, ρ, ω, γ, ℓ; atpi = false)
    c=cos(γ)
    omc=2sin(γ/2)^2                          # 1 - c
    lg=log(ℓ)
    lg1=log(ω/(1+ρ))                         # log(1-ρ)
    sω=sqrt(ω)
    at=atpi ? oftype(γ, π)/2 : atan((1+ρ)*tan(γ/2)/sω)   # A = 2at/√ω
    if kind==:inv
        a1=(3ρ^2-4)/ρ^4
        poly=-omc*(4c^2*ρ^2+4c*ρ^2+6c*ρ-5ρ^2+6ρ+12)/(3ρ^3)
        I1=a1*(lg1-lg)+poly
        a=((4, -4.0), (2, 5.0), (0, -1.0))   # S(c) = sin3τ sinτ
        I2=(ρ^2-4)*sω/ρ^4*2at+trigsum(a, γ, ρ)
    elseif kind==:log
        a1=ω*(ρ^2-2)/(2ρ^4)
        ac=ℓ*(1+ρ*c)*(2c^2*ρ^2-3ρ^2+2)/(2ρ^4)
        poly=-omc*(3c^3*ρ^3+3c^2*ρ^3+4c^2*ρ^2-6c*ρ^3+4c*ρ^2+6c*ρ-6ρ^3-14ρ^2+6ρ+12)/(12ρ^3)
        I1=a1*lg1+ac*lg+poly
        P=sin(2γ)/4-sin(4γ)/8                # ∫₀^γ S(cos τ) dτ;  P sin τ = c(1-c²)²
        a=((1, 1.0), (3, -2.0), (5, 1.0))
        I2=P*lg-ρ*(ω*sω/ρ^5*2at+trigsum(a, γ, ρ))
    else # :xlog
        a1=(ω/(1+ρ))^2*(2ρ^3-ρ^2-4ρ-2)/(10ρ^4)
        ac=ℓ^2*(8c^3*ρ^3+6c^2*ρ^2-10c*ρ^3+4c*ρ-5ρ^2+2)/(10ρ^4)
        poly=omc*(48c^4*ρ^4+48c^3*ρ^4-15c^3*ρ^3-52c^2*ρ^4-15c^2*ρ^3-20c^2*ρ^2-52c*ρ^4+
                  60c*ρ^3-20c*ρ^2-30c*ρ-52ρ^4+60ρ^3+130ρ^2-30ρ-60)/(300ρ^3)
        I1=a1*lg1+ac*lg+poly
        P3=sin(2γ)/4-sin(4γ)/8-ρ*(sin(γ)/4-sin(5γ)/20)
        a=((1, 1.0), (3, -2.0), (5, 1.0), (0, -ρ/5), (2, -2ρ/5), (4, 7ρ/5), (6, -4ρ/5))
        I2=P3*lg-ρ*(ω^2*sω/(5ρ^5)*2at+trigsum(a, γ, ρ))
    end
    -I1/4, -I2/4
end
# U_f at x for p = (p₁, p₂, 0, 0) inside the ball, with ω = 1-|p|² supplied accurately.
function Uf(kind, p1, p2, ω, x)
    ρ=sqrt(p1^2+p2^2)
    d2=(x[1]-p1)^2+(x[2]-p2)^2+x[3]^2+x[4]^2   # |x-p|²
    ℓ=(d2+ω)/2                                 # 1 - p·x
    sγ2=(d2-(1-ρ)^2)/(4ρ)                      # sin²(γ/2) = (1-c)/2, c = p̂·x
    γ=2asin(sqrt(clamp(sγ2, zero(sγ2), one(sγ2))))
    zonal(kind, ρ, ω, γ, ℓ)
end

const B41=(π-2)/(3π)
# Integrand at p = (σt, 1-t). Near t = 0 use variables (p₁, q = 1-p₂), near t = 1
# (r = 1-|p₁|, p₂), so that ω = 1-|p|² = 2t(1-t) never cancels.
# Integrand at p = (σt, 1-t), given t and 1-t separately. Near t = 0 use variables
# (p₁, q = 1-p₂), near t = 1 (r = 1-|p₁|, p₂), so that ω = 1-|p|² = 2t(1-t) never cancels.
function integrand(t, omt, x, σ)
    Flo(v, kind)=Uf(kind, v[1], 1-v[2], 2v[2]-v[2]^2-v[1]^2, x)        # v = (p₁, 1-p₂)
    Fhi(v, kind)=Uf(kind, σ*(1-v[1]), v[2], 2v[1]-v[1]^2-v[2]^2, x)    # v = (1-|p₁|, p₂)
    F, v0, J=t<=1/2 ? (Flo, [σ*t, t], [1 0; 0 -1]) : (Fhi, [omt, omt], [-σ 0; 0 1])
    J=convert(Matrix{typeof(t)}, J)
    d1=J*ForwardDiff.gradient(v->F(v, :log), v0)
    H=J*ForwardDiff.hessian(v->F(v, :xlog), v0)*J
    -F(v0, :inv)/3+(7/3)*d1[2]+(5/3)*H[2, 2]+σ*H[1, 2]
end
integrand(t, x, σ)=integrand(t, 1-t, x, σ)

# Fixed rule in φ (t = sin²φ) on panels graded geometrically toward φ = 0 and π/2,
# where the integrand grows like log φ.
function feynman_rule(; n=12, levels=14, q=0.25)
    pts=sort(unique(vcat([π/4*q^k for k in 0:levels], [π/2-π/4*q^k for k in 0:levels], [0.0, π/4, π/2])))
    x, w=gauss(n)
    φs=Float64[]; ws=Float64[]
    for i in 1:length(pts)-1
        h=(pts[i+1]-pts[i])/2; m=(pts[i+1]+pts[i])/2
        append!(φs, m.+h.*x); append!(ws, h.*w)
    end
    φs, ws
end
const FEYNMAN_RULE=feynman_rule()
function psi41_z2_feynman(α, θ; rule=FEYNMAN_RULE)
    x=[cos(α), sin(α)*cos(θ), sin(α)*sin(θ), 0.0]
    φs, ws=rule
    s=0.0
    for (φ, w) in zip(φs, ws)
        t, omt=sin(φ)^2, cos(φ)^2
        s+=2w*(integrand(t, omt, x, 1)+integrand(t, omt, x, -1))
    end
    B41/(sqrt(2)*π)*s
end
