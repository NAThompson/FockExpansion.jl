module FockExpansion
using QuadGK, SpecialFunctions, ForwardDiff
primal(x) = x
primal(x::ForwardDiff.Dual) = primal(ForwardDiff.value(x))
adnorm(x::Real) = abs(x)
adnorm(x::ForwardDiff.Dual) =
    max(adnorm(ForwardDiff.value(x)), maximum(adnorm, ForwardDiff.partials(x)))
typedpi(x) = oftype(primal(float(x)), Base.MathConstants.pi)
export psi00, psi10, psi20, psi21, psi31, psi30, psi30_parts, clausen2

"""Real Clausen function Cl₂(x), evaluated by a convergent series after periodic reduction."""
function clausen2(x::Real)
    x = float(x)
    π = typedpi(x)
    u = x - round(x/(2π))*(2π)
    iszero(u) && return zero(u)
    q = (u/(2π))^2
    ans = u*(1-log(abs(u)))
    power = u*q
    for n = 1:max(100, precision(x))
        term = zeta(oftype(x, 2n))*power/(n*(2n+1))
        ans += term
        abs(term)<eps(x)*max(abs(ans), one(x))/4 && return ans
        power *= q
    end
    error("Clausen series failed to converge")
end
# Fixed coefficients remove repeated zeta calls in the Float64 hot path.
const CLAUSEN64 = ntuple(n -> zeta(Float64(2n))/(n*(2n+1)), 32)
function clausen2(x::Float64)
    u=rem2pi(x, RoundNearest)
    iszero(u) && return zero(u)
    q=(u/(2π))^2
    u*(1-log(abs(u)))+u*q*evalpoly(q, CLAUSEN64)
end
L(u) = iszero(u) ? zero(u) : u*log(abs(2sin(u)))+clausen2(2u)/2
T(u) = -u*log(abs(2cos(u)))+clausen2(typedpi(u)-2u)/2
# Analytic scalar derivatives of the elementary/Clausen primitives. Nested
# duals propagate these rules to second derivatives without log(0) cancellation.
function clausen2(x::ForwardDiff.Dual{Tag}) where {Tag}
    u=ForwardDiff.value(x)
    ForwardDiff.Dual{Tag}(clausen2(u), -log(abs(2sin(u/2)))*ForwardDiff.partials(x))
end
function L(x::ForwardDiff.Dual{Tag}) where {Tag}
    u=ForwardDiff.value(x)
    derivative=iszero(primal(u)) ? one(u)-u*u/3 : u/tan(u)
    ForwardDiff.Dual{Tag}(L(u), derivative*ForwardDiff.partials(x))
end
function T(x::ForwardDiff.Dual{Tag}) where {Tag}
    u=ForwardDiff.value(x)
    ForwardDiff.Dual{Tag}(T(u), u*tan(u)*ForwardDiff.partials(x))
end
function geometry(α, θ)
    α, θ=promote(float(α), float(θ))
    π=typedpi(α)
    (0<α<π && 0<θ<π) || throw(
        DomainError(
            (α, θ),
            "This evaluator requires 0 < α, θ < π; boundary limits are not implemented.",
        ),
    )
    r1, r2=cos(α/2), sin(α/2)
    v=sin(α)*cos(θ)
    (; α, θ, r1, r2, σ = r1+r2, v, ξ = sqrt(1-v))
end
psi00(α, θ) = (geometry(α, θ); 1.0)
function psi10(α, θ; Z)
    g=geometry(α, θ)
    g.ξ/2-Z*g.σ
end
# All dilogarithm arguments in χ lie on the unit circle.
function unit_dilog(z)
    π=typedpi(real(z))
    u=mod(angle(z), 2π)
    complex(π^2/6-u*(2π-u)/4, clausen2(u))
end
function chi(α, μ)
    π=typedpi(α)
    x=min(α, π-α)
    s=sin(x)
    c=cos(x)
    v=s*μ
    σ=sqrt(1+s)
    ξ=sqrt(1-v)
    g=sqrt(1-v^2)
    h=-(3π+10-16clausen2(π/2))/(24π)
    if abs(c)<8eps(primal(float(α)))
        return (g+1-2sqrt(oftype(primal(α), 2))*ξ)/6-μ*log(sqrt(oftype(primal(α), 2))+ξ)/3+g*asin(
            μ,
        )/(3π)+h*μ
    end
    β=abs(μ)==1 ? μ*x : asin(v)
    dl=clausen2(x-β)+clausen2(π+β-x)-clausen2(π-x-β)-clausen2(x+β)
    logs=β*(2log(s)-2log(g+c))
    logs+=μ == -1 ? zero(β) : (β+x)*log1p(μ)
    logs+=μ == 1 ? zero(β) : (β-x)*log1p(-μ)
    (g+s-2σ*ξ)/6+(-2π*v*log(σ+ξ)+π*c*log((c+σ*ξ)^2/(σ^2*(g+c)))+2g*β+c*(logs-dl))/(6π)+h*v

end
"""Closed classical expression for ψ₂₀, with total-projection a₂₁."""
function psi20(α, θ; Z, E, a21)
    π=typedpi(α)
    g=geometry(α, θ)
    (1-2E)/12+Z*(chi(g.α, cos(g.θ))-g.σ*g.ξ/3)+Z^2*(one(g.α)/3+sin(g.α)/2)+(
        a21-Z*(π+4)/(9π)
    )*g.v
end
function psi21(α, θ; Z)
    π=typedpi(α)
    g=geometry(α, θ)
    -Z*(π-2)*g.v/(3π)
end
function psi31(α, θ; Z)
    π=typedpi(α)
    g=geometry(α, θ)
    Z*(π-2)*g.ξ*(5g.ξ^2-6)/(36π)+Z^2*(π-2)*g.σ*g.v/(6π)
end
function K(m, d, s, x, h, Z)
    π=typedpi(d)
    for _ = 1:mod(m, 4)
        d, s=-s, d
    end
    Z*(d*(-9s^2+6s*x+8x^2-8)-h*(5-2x^2))/(144π)-Z^2*5*(x^2-1)*(d+h)/(32π)
end
function Ncoef(m, d, s, x, h)
    π=typedpi(d)
    ν = m in (0, 2) ? x^2-s^2 : m==1 ? 9s^2+x^2-2-8s*h : 9s^2+x^2-2+8s*h
    d*ν/(36π)
end
elam(d, s, x, h, Z) =
    Z*(7s*d*h+9s^3+8s*x^2-26s+10x^3-12x)/72+Z^2*s*(137d*h+64s^2-154x^2+26)/288
function classical(α, p0, Z)
    π=typedpi(α)
    r1, r2=cos(α/2), sin(α/2)
    s, d=r1+r2, r1-r2
    x, h0=sqrt(oftype(primal(α), 2))*sin(p0/2), sqrt(oftype(primal(α), 2))*cos(p0/2)
    a, p=α-π/2, p0-π
    G=clausen2(π/2)
    l2=log(oftype(primal(α), 2))
    v=0.0
    for sign in (1, -1)
        w=(a+sign*p)/4
        h=sign*h0
        kk(m) = K(m, d, s, x, h, Z)+Z^2*Ncoef(m, d, s, x, h)
        v+=sum(4kk(m)*L(w+m*π/4) for m in (-1, 0, 1))-4kk(2)*T(w)+elam(d, s, x, h, Z)*log(
            abs(2cos(w)),
        )
    end
    c9=-h0*(Z*(2x^2-5)/(9π)-5Z^2*(x^2-1)/(2π))
    v+=c9*L(p/2)/2
    v+=2Z^2/(9π)*(r2*(8r2^2+x^2-5)*L(α/2)+r1*(8r1^2+x^2-5)*L((π-α)/2))
    qa2=Z*x*(2x^2+1)/(144π)-5Z^2*x*(x^2-1)/(32π)
    qa=7Z*x*d*s/144+Z^2*(x*d*s/64+d*(x^2-1)/(9π))
    qp2=Z^2*s*(4s^2+x^2-9)/(36π)
    qp=-Z*h0*(5x^2-4)/(36π)+Z^2*s*x*h0/(18π)
    v+=qa2*a^2+qa*a+qp2*p^2+qp*p
    e0=x*(3-2x^2)/72
    e1=-x*(5x^2-6)*(l2+4G/π)/72-π*x*(1-x^2)/36+x*(91x^2-106)/(144π)+(
        4s^3+4s^2*x-32s*x^2+4s-27x^3+26x
    )/288
    e2=-s*(1-x^2)*(l2+4G/π)/12+53s*(1-x^2)/(72π)+π*(16s^3-4s*x^2-28s-45x^3+45x)/576-(
        503s^3-375s^2*x+384s*x^2-1212s-656x^3+960x
    )/864
    e3=-s*(5s^2-3)/36
    v+e0+Z*e1+Z^2*e2+Z^3*e3
end
function panel(b, y, c0, c2; rtol)
    π=typedpi(b)
    h=sqrt((1+b)/(1-b))
    # Chart decisions use only primal values, never dual partials at a tie.
    upper=primal(h)<=1 ? h : one(h)
    function integrand(q)
        q == 0 && return (2/π)*(c0*π^2/4+c2)/sqrt(1+b)
        a=atan(1/q)
        d=atan(q/h)
        numerator=c0*(a-d)*(a+d)+c2*((h*d/q)^2-(q*a)^2)
        # Positive form of the transformed quartic avoids cancellation near y=0.
        T=(1-q^2)*((1+b)-(1-b)*q^2)+4*y^2*q^2
        (2/π)*numerator/sqrt(T)
    end
    # Fixed real endpoints let duals carry the moving-upper-limit derivative.
    lo=zero(primal(h))
    hi=one(lo)
    quadgk(
        t->upper*integrand(upper*t),
        lo,
        hi/2,
        hi;
        rtol,
        atol = rtol/100,
        norm = adnorm,
        order = primal(h) isa BigFloat ? 21 : 7,
    )
end
"""Decompose ψ₃₀ into its classical, E/a₂₁, and elliptic contributions.
`quadrature_error` is an adaptive quadrature estimate, not a rigorous error bound.
Inputs are dimensionless and a₂₁ uses the total-projection convention.
"""
function psi30_parts(α, θ; Z, E, a21, rtol = 1e-11)
    π=typedpi(α)
    g=geometry(α, θ)
    (; r1, r2, σ, v, ξ)=g
    α=primal(g.α)<=π/2 ? g.α : π-g.α
    cl=classical(α, acos(v), Z)
    state=E*(Z*σ*(2+r1*r2)/18-(6-ξ^2)*ξ/72)-a21*(Z*σ*v/2-(6-5ξ^2)*ξ/12)
    ell=0.0
    err=0.0
    if Z!=0
        for b in (cos(α), -cos(α))
            y=-ξ/sqrt(oftype(primal(α), 2))
            c0=-sqrt(oftype(primal(α), 2))*Z/576*(9Z*b^2-18Z-8b^2+8)
            c2=sqrt(oftype(primal(α), 2))*Z/288*(b-1)*(90Z*y^2-45Z-18b-32y^2+16)
            val, e=panel(b, y, c0, c2; rtol)
            ell+=val
            err+=e
        end
        for y in (-r1, -r2)
            b=-v
            val, e=panel(b, y, Z^2*(1-b^2)/18, 2Z^2*(b-1)*(2y^2-1)/9; rtol)
            ell+=val
            err+=e
        end
    end
    (; classical = cl, state, elliptic = ell, value = cl+state+ell, quadrature_error = err)
end
"""Numerical ψ₃₀ at an interior angle. Uses four one-dimensional elliptic quadratures."""
psi30(α, θ; kwargs...) = psi30_parts(α, θ; kwargs...).value
include("green.jl")
include("langner.jl")
export LangnerTable, psi30_langner
export psi30_green
end
