module FockExpansion
using QuadGK, SpecialFunctions, ForwardDiff
using ClausenFunctions: cl2
primal(x) = x
primal(x::ForwardDiff.Dual) = primal(ForwardDiff.value(x))
adnorm(x::Real) = abs(x)
adnorm(x::ForwardDiff.Dual) =
    max(adnorm(ForwardDiff.value(x)), maximum(adnorm, ForwardDiff.partials(x)))
typedpi(x) = oftype(primal(float(x)), Base.MathConstants.pi)
export psi00, psi10, psi20, psi21, psi31, psi30, psi30_parts, psi41, psi42

# L(u) = u log|2 sin u| + Cl₂(2u)/2 and T(u) = -u log|2 cos u| + Cl₂(π-2u)/2,
# given s = sin u and c = cos u so that callers can share them.
L(u, s, c) = iszero(u) ? zero(u) : u*log(abs(2s))+cl2(2u)/2
T(u, s, c) = -u*log(abs(2c))+cl2(typedpi(u)-2u)/2
L(u) = L(u, sincos(u)...)
T(u) = T(u, sincos(u)...)
# Analytic scalar derivatives L′ = u cot u and T′ = u tan u. Nested
# duals propagate these rules to second derivatives without log(0) cancellation.
function L(x::ForwardDiff.Dual{Tag}, s, c) where {Tag}
    u, s, c=ForwardDiff.value(x), ForwardDiff.value(s), ForwardDiff.value(c)
    derivative=iszero(primal(u)) ? one(u)-u*u/3 : u*c/s
    ForwardDiff.Dual{Tag}(L(u, s, c), derivative*ForwardDiff.partials(x))
end
function T(x::ForwardDiff.Dual{Tag}, s, c) where {Tag}
    u, s, c=ForwardDiff.value(x), ForwardDiff.value(s), ForwardDiff.value(c)
    ForwardDiff.Dual{Tag}(T(u, s, c), u*s/c*ForwardDiff.partials(x))
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
    complex(π^2/6-u*(2π-u)/4, cl2(u))
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
    h=-(3π+10-16cl2(π/2))/(24π)
    if abs(c)<8eps(primal(float(α)))
        return (g+1-2sqrt(oftype(primal(α), 2))*ξ)/6-μ*log(sqrt(oftype(primal(α), 2))+ξ)/3+g*asin(
            μ,
        )/(3π)+h*μ
    end
    β=abs(μ)==1 ? μ*x : asin(v)
    dl=cl2(x-β)+cl2(π+β-x)-cl2(π-x-β)-cl2(x+β)
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
    rt=1/sqrt(oftype(primal(α), 2))
    sp0, cp0=sincos(p0/2)
    x, h0=sp0/rt, cp0/rt
    a, p=α-π/2, p0-π
    G=oftype(primal(α), Base.MathConstants.catalan) # Cl₂(π/2)
    l2=log(oftype(primal(α), 2))
    v=0.0
    for sign in (1, -1)
        w=(a+sign*p)/4
        h=sign*h0
        sw, cw=sincos(w)
        # sin and cos of w ∓ π/4
        sm, cm=(sw-cw)*rt, (cw+sw)*rt
        sp, cp=(sw+cw)*rt, (cw-sw)*rt
        kk(m) = K(m, d, s, x, h, Z)+Z^2*Ncoef(m, d, s, x, h)
        v+=4kk(-1)*L(w-π/4, sm, cm)+4kk(0)*L(w, sw, cw)+4kk(1)*L(w+π/4, sp, cp)
        v+=-4kk(2)*T(w, sw, cw)+elam(d, s, x, h, Z)*log(abs(2cw))
    end
    # sin(p/2) = -cos(p0/2), cos(p/2) = sin(p0/2), sin(α/2) = r2, cos(α/2) = r1
    c9=-h0*(Z*(2x^2-5)/(9π)-5Z^2*(x^2-1)/(2π))
    v+=c9*L(p/2, -cp0, sp0)/2
    v+=2Z^2/(9π)*(r2*(8r2^2+x^2-5)*L(α/2, r2, r1)+r1*(8r1^2+x^2-5)*L((π-α)/2, r1, r2))
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
# The elliptic moments are integrals over 0 < q < upper = min(h, 1), h = √((1+b)/(1-b)).
# Write q = upper*t. The quartic under the square root vanishes at t = 1; as y → 0
# (always together with b → 0) the integrand develops a peak of width σ ~ |y| there.

# Quartic (1-q²)((1+b)-(1-b)q²) at q = upper*t, with m = 1-t² computed as (1-t)(1+t).
# Written as a sum of nonnegative terms, so the zero at t = 1 costs no relative precision.
# The caller states which limit applies, so primal ties never switch formulas under duals.
quartic(b, h, m, upper_is_h) = upper_is_h ? m*h^2*((1+b)*m-2b) : m*((1-b)*m+2b)

# Moment integrand c0*N₀ + c2*N₂ at q, given a = atan(1/q), d = atan(q/h) and the quartic P.
function moment_integrand(q, a, d, h, b, y, c0, c2, P)
    π=typedpi(q)
    iszero(q) && return (2/π)*(c0*π^2/4+c2)/sqrt(1+b)
    numerator=c0*(a-d)*(a+d)+c2*((h*d/q)^2-(q*a)^2)
    (2/π)*numerator/sqrt(P+4*y^2*q^2)
end

# ∫₀¹ f(t, 1-t) dt for an integrand with a peak of width σ at t = 1, as one quadrature
# over s ∈ [0, 2] with a break at s = 1. The bulk t ∈ [0, 1/2] stays linear in s, since a
# stretching map there would crowd the branch points at t = ±i. A narrow peak gets
# x = 1-t = σ sinh(τ) on the tail, which turns it into a smooth ramp. A single quadgk
# call keeps one compiled specialization per integrand, which matters for nested duals.
# Fixed real endpoints let duals carry the derivatives of σ and of the moving upper limit.
function endpoint_quadgk(f, σ, lo; rtol)
    hi=one(lo)
    narrow=primal(σ)<1/4
    τmax=narrow ? asinh((hi/2)/σ) : zero(σ)
    function integrand(s)
        s<=hi && return f(s/2, 1-s/2)/2
        u=s-hi
        narrow || return f(1-u/2, u/2)/2
        τ=τmax*u
        x=σ*sinh(τ)
        σ*τmax*cosh(τ)*f(1-x, x)
    end
    order=primal(lo) isa BigFloat ? 21 : 7
    quadgk(integrand, lo, hi, 2hi; rtol, atol = rtol/100, norm = adnorm, order)
end

"""One elliptic moment combination c0*N₀(b, y) + c2*N₂(b, y), with its quadrature error estimate,
by adaptive quadrature."""
function panel_adaptive(b, y, c0, c2; rtol)
    π=typedpi(b)
    h=sqrt((1+b)/(1-b))
    # Chart decisions use only primal values, never dual partials at a tie.
    upper_is_h=primal(h)<=1
    upper=upper_is_h ? h : one(h)
    function integrand(t, x)
        q=upper*t
        P=quartic(b, h, x*(1+t), upper_is_h)
        upper*moment_integrand(q, π/2-atan(q), atan(q/h), h, b, y, c0, c2, P)
    end
    endpoint_quadgk(integrand, abs(y), zero(primal(h)); rtol)
end

# The four moment combinations of ψ₃₀ as two quadratures in t, one per pair of panels
# that share both their arctangents and their peak width. The pair b = ±cos(α) has
# h₂ = 1/h₁, so both panels need only atan(t) and atan(t/h₁); the pair with b = -v
# shares q, a and d outright. Each pair costs two arctangents per node instead of four.
function elliptic(α, g, Z; rtol)
    π=typedpi(α)
    r=sqrt(oftype(primal(α), 2))
    lo=zero(primal(α))
    # Panels 1 and 2: b = ±c with c = cos(α) ≥ 0 because α ≤ π/2, so h₁ ≥ 1 ≥ h₂.
    c=cos(α)
    y=-g.ξ/r
    h1=sqrt((1+c)/(1-c))
    h2=inv(h1)
    c0₁₂=-r*Z/576*(9Z*c^2-18Z-8c^2+8)
    c2₁=r*Z/288*(c-1)*(90Z*y^2-45Z-18c-32y^2+16)
    c2₂=r*Z/288*(-c-1)*(90Z*y^2-45Z+18c-32y^2+16)
    function pair12(t, x)
        m=x*(1+t)
        A=atan(t)
        D=atan(t/h1)
        moment_integrand(t, π/2-A, D, h1, c, y, c0₁₂, c2₁, quartic(c, h1, m, false))+
        h2*moment_integrand(h2*t, π/2-D, A, h2, -c, y, c0₁₂, c2₂, quartic(-c, h2, m, true))
    end
    # Panels 3 and 4: b = -v, y = -r1 and y = -r2.
    b=-g.v
    h=sqrt((1+b)/(1-b))
    upper_is_h=primal(h)<=1
    upper=upper_is_h ? h : one(h)
    c0₃₄=Z^2*(1-b^2)/18
    c2₃=2Z^2*(b-1)*(2g.r1^2-1)/9
    c2₄=2Z^2*(b-1)*(2g.r2^2-1)/9
    function pair34(t, x)
        q=upper*t
        a=π/2-atan(q)
        d=atan(q/h)
        P=quartic(b, h, x*(1+t), upper_is_h)
        upper*(
            moment_integrand(q, a, d, h, b, -g.r1, c0₃₄, c2₃, P)+
            moment_integrand(q, a, d, h, b, -g.r2, c0₃₄, c2₄, P)
        )
    end
    v12, e12=endpoint_quadgk(pair12, -y, lo; rtol)
    v34, e34=endpoint_quadgk(pair34, min(g.r1, g.r2), lo; rtol)
    v12+v34, e12+e34
end
# Fixed-rule evaluation of the panel integral ∫₀ᵘ (2/π) n(q)/√T(q) dq.
# The quartic factors exactly as T = (1-b)(q²-s₁)(q²-s₂) = Q(q)F(q), where the
# real quadratic Q = (q-m)² ± t² holds the two roots nearest the endpoint u
# and F > 0 on [0,u]. As y → 0 those roots pinch onto u and the integrand
# becomes nearly singular. The substitution q = m + t sinh z (or m - t cosh z)
# gives dq/√Q = dz exactly, so the z-integrand n(q)/√F(q) is smooth and a
# fixed Gauss–Legendre rule converges, with a node count growing like
# log(1/|y|). Differences that vanish as y → 0 use T(u) = 4y²u², i.e.
# (u²-s₁)(u²-s₂) = 4y²u²/(1-b), so none is formed by cancellation.
struct PanelMap{T}
    h::T
    za::T
    zb::T
    m::T
    t::T
    cosh_branch::Bool # q = m - t cosh z, else q = m + t sinh z
    f1::T             # F(q) = (1-b)(q² + f1 q + f0)
    f0::T
end
acosh1p(δ) = log1p(δ+sqrt(δ*(2+δ))) # acosh(1+δ), accurate for small δ
function panelmap(b, y)
    h2=(1+b)/(1-b)
    h=sqrt(h2)
    u=primal(h)<=1 ? h : one(h)
    u2=u*u
    # B-u², with B=(s₁+s₂)/2; fma forms b∓2y² with a single rounding.
    Bu=primal(h)<=1 ? -fma(2y, y, b)/(1-b) : fma(-2y, y, b)/(1-b)
    P=4y^2*u2/(1-b) # (u²-s₁)(u²-s₂)
    disc=fma(Bu, Bu, -P) # B²-s₁s₂
    if primal(disc)<0
        # Complex pair r, r̄ with r=√s₁: Q=(q-m)²+t², m=Re r, t=Im r.
        B=Bu+u2
        sq=sqrt(-disc)
        a=hypot(B, sq) # |s₁|
        if primal(B)>=0
            rr=sqrt((a+B)/2)
            ri=sq/(2rr)
        else
            ri=sqrt((a-B)/2)
            rr=sq/(2ri)
        end
        # ρ=r-u=(s₁-u²)/(r+u), with s₁-u²=Bu+i sq
        den=(rr+u)^2+ri^2
        ρr=(Bu*(rr+u)+sq*ri)/den
        t=ri
        m=u+ρr
        PanelMap(h, asinh(-m/t), asinh(-ρr/t), m, t, false, 2rr, a)
    elseif primal(Bu+u2)>0
        # Two roots u² < s₁ ≤ s₂: Q=(q-r₁)(q-r₂)=(q-m)²-t², r₁=√s₁ ≥ u.
        sd=sqrt(disc)
        δ2=-(Bu+sd) # u²-s₂
        δ1=P/δ2     # u²-s₁ without cancellation
        r1=sqrt(u2-δ1)
        r2=sqrt(u2-δ2)
        t=sd/(r1+r2) # (r₂-r₁)/2
        ε1=-δ1/(r1+u) # r₁-u
        PanelMap(h, acosh1p(ε1/t), acosh1p(r1/t), r1+t, t, true, r1+r2, r1*r2)
    else
        # Two negative roots: Q=q²-s_near, F=(1-b)(q²-s_far).
        sd=sqrt(disc)
        sfar=Bu+u2-sd
        t=sqrt(-h2/sfar) # √(-s_near), s_near s_far=h²
        PanelMap(h, zero(t), asinh(u/t), zero(t), t, false, zero(t), -sfar)
    end
end
# Gauss–Legendre rules on [-1,1] with 1, …, 96 nodes, padded with zero-weight
# nodes at 0 to a multiple of 16 entries: the vectorized node loop handles 16
# nodes per iteration and falls back to scalar code for any remainder.
function padded_gauss(n)
    x, w=gauss(n)
    k=16cld(n, 16)-n
    vcat(x, zeros(k)), vcat(w, zeros(k))
end
const PANEL_RULES=[padded_gauss(n) for n = 1:96]
# Node count for z-range L and d=-log₁₀(rtol) requested digits: an upper
# envelope of the nodes needed against 200-bit references over a 31×32 grid of
# b∈[-0.9999,0.9999], y∈[-1,-1e-10]. Relative accuracy is limited to about
# 1e-14 (worst near b=±0.9999), so d is capped at 14.
function panel_nodes(L, rtol)
    d=clamp(-log10(rtol), 4, 14)
    n=(0.83d+4.3)+(0.18d-0.35)*L
    clamp(4ceil(Int, n/4), 8, 96)
end
const PANEL_ERROR=1e-14
include("panel_node_table.jl")
# Nodes for a panel from the calibrated table over (atanh(b), L), or 0 outside it.
# The table resolves interior points, where the limiting singularities are the
# arctangent branch points rather than the pinching roots, far more tightly than
# the envelope in `panel_nodes`; the tolerance level used is the first tabulated
# one at or below rtol.
function panel_nodes_table(b, L, rtol)
    β=atanh(b)
    β0, β1, nβ=PANEL_NODE_BETA
    L0, L1, nL=PANEL_NODE_L
    (β0<=β<=β1 && L<=L1) || return 0
    k=something(findfirst(<=(rtol), PANEL_NODE_LEVELS), length(PANEL_NODE_LEVELS))
    i=clamp(floor(Int, (β-β0)/(β1-β0)*nβ)+1, 1, nβ)
    j=clamp(floor(Int, (L-L0)/(L1-L0)*nL)+1, 1, nL)
    Int(PANEL_NODE_TABLE[k, i, j])
end
# Branch-free kernels for the Float64 node loop, so that it vectorizes.
# @horner(x, c₀, c₁, …) expands c₀ + x(c₁ + x(…)) inline (evalpoly is a call).
macro horner(x, cs...)
    ex=esc(cs[end])
    for c in reverse(cs[1:(end-1)])
        ex=:(muladd(t, $ex, $(esc(c))))
    end
    :(let t=$(esc(x)); $ex; end)
end
# atan01 is fdlibm's atan (as in Base) restricted to 0 ≤ x ≤ 1, with the
# reduction intervals [0,7/16), [7/16,11/16), [11/16,1] selected by ifelse.
@inline function atan01(x::Float64)
    lo=x<7/16
    mid=!lo&(x<11/16)
    num=ifelse(lo, x, ifelse(mid, 2x-1, x-1))
    den=ifelse(lo, 1.0, ifelse(mid, 2+x, x+1))
    hi=ifelse(lo, 0.0, ifelse(mid, 4.63647609000806093515e-01, 7.85398163397448278999e-01))
    lo_=ifelse(lo, 0.0, ifelse(mid, 2.26987774529616870924e-17, 3.06161699786838301793e-17))
    t=num/den
    t2=t*t
    t4=t2*t2
    p=t2*@horner(t4, 3.33333333333329318027e-01, 1.42857142725034663711e-01,
        9.09088713343650656196e-02, 6.66107313738753120669e-02,
        4.97687799461593236017e-02, 1.62858201153657823623e-02)
    q=t4*@horner(t4, -1.99999999998764832476e-01, -1.11111104054623557880e-01,
        -7.69187620504482999495e-02, -5.83357013379057348645e-02,
        -3.65315727442169155270e-02)
    hi-((t*(p+q)-lo_)-t)
end
# exp(z) for |z| ≤ 700 without branches: z = k log 2 + r, |r| ≤ log(2)/2,
# and the Taylor polynomial of degree 13 for exp(r) (error < 4e-18 relative).
@inline function exp_nobranch(z::Float64)
    k=round(z*1.4426950408889634)
    r=fma(-k, 6.93147180369123816490e-01, z)
    r=fma(-k, 1.90821492927058770002e-10, r)
    e=@horner(r, 1.0, 1.0, 1/2, 1/6, 1/24, 1/120, 1/720, 1/5040, 1/40320, 1/362880,
        1/3628800, 1/39916800, 1/479001600, 1/6227020800)
    e*reinterpret(Float64, (unsafe_trunc(Int64, k)+1023)<<52)
end
# Float64 node loop: q = m + t sinh z (or m - t cosh z) via one exp, and
# atan(1/q) = π/2 - atan(q) for 0 < q ≤ 1.
function panel_fixed(b::Float64, y::Float64, c0::Float64, c2::Float64, map::PanelMap{Float64}, n)
    (; h, za, zb, m, t, cosh_branch, f1, f0)=map
    x, w=PANEL_RULES[n]
    half=(zb-za)/2
    mid=(zb+za)/2
    ih=1/h
    st=cosh_branch ? -t/2 : t/2      # q = m + st*(e ± 1/e)
    sg=cosh_branch ? 1.0 : -1.0
    acc=0.0
    @inbounds @simd for i in eachindex(x)
        e=exp_nobranch(mid+half*x[i])
        q=m+st*(e+sg/e)
        iq=1/q
        a=1.5707963267948966-atan01(q)
        d=atan01(q*ih)
        numerator=c0*(a-d)*(a+d)+c2*((h*d*iq)^2-(q*a)^2)
        # F(q)/(1-b) > 0; abs lets the compiler drop sqrt's domain check
        acc+=w[i]*numerator/sqrt(abs(q*(q+f1)+f0))
    end
    (2/pi)*half*acc/sqrt(1-b)
end
function panel_fixed(b, y, c0, c2, map::PanelMap, n)
    π=typedpi(b)
    (; h, za, zb, m, t, cosh_branch, f1, f0)=map
    x, w=PANEL_RULES[n]
    half=(zb-za)/2
    mid=(zb+za)/2
    ih=1/h
    acc=zero(half*c0)
    @inbounds for i in eachindex(x)
        z=mid+half*x[i]
        q=cosh_branch ? m-t*cosh(z) : m+t*sinh(z)
        iq=1/q
        a=atan(iq)
        d=atan(q*ih)
        numerator=c0*(a-d)*(a+d)+c2*((h*d*iq)^2-(q*a)^2)
        acc+=w[i]*numerator/sqrt(q*(q+f1)+f0)
    end
    (2/π)*half*acc/sqrt(1-b)
end
"""Panel integral ∫₀ᵘ (2/π)(c₀n₀(q)+c₂n₂(q))/√T(q) dq for one pair (b,y).
Returns (value, error estimate). Float64 inputs (including ForwardDiff duals)
use a fixed Gauss–Legendre rule after removing the endpoint near-singularity;
other number types use adaptive quadrature with tolerance `rtol`.
"""
function panel(b, y, c0, c2; rtol)
    b, y=promote(float(b), float(y))
    if primal(b) isa Float64
        map=panelmap(b, y)
        L=primal(map.zb-map.za)
        # The double-root limit t → 0 (measure zero) and extreme ranges fall back.
        if isfinite(L) && L<=60
            # The table is calibrated on values; derivatives need more nodes, so
            # duals keep the envelope.
            n=b isa Float64 ? panel_nodes_table(b, L, rtol) : 0
            iszero(n) && (n=panel_nodes(L, rtol))
            val=panel_fixed(b, y, c0, c2, map, n)
            return val, max(rtol, PANEL_ERROR)*abs(primal(val))
        end
    end
    panel_adaptive(b, y, c0, c2; rtol)
end
# The four moment combinations of ψ₃₀ for Float64 inputs, one fixed-rule panel each.
# Each panel has its own map z → q, so unlike `elliptic` the pairs cannot share
# arctangents; removing the near-singularity analytically more than pays for that
# except at interior points with loose rtol.
function elliptic_fixed(α, g, Z; rtol)
    r=sqrt(oftype(primal(α), 2))
    y=-g.ξ/r
    ell=zero(α*Z)
    err=zero(primal(ell))
    for b in (cos(α), -cos(α))
        c0=-r*Z/576*(9Z*b^2-18Z-8b^2+8)
        c2=r*Z/288*(b-1)*(90Z*y^2-45Z-18b-32y^2+16)
        val, e=panel(b, y, c0, c2; rtol)
        ell+=val
        err+=e
    end
    b=-g.v
    for yr in (-g.r1, -g.r2)
        val, e=panel(b, yr, Z^2*(1-b^2)/18, 2Z^2*(b-1)*(2yr^2-1)/9; rtol)
        ell+=val
        err+=e
    end
    ell, err
end
"""Decompose ψ₃₀ into its classical, E/a₂₁, and elliptic contributions.
`quadrature_error` is an estimate, not a rigorous error bound: for Float64 inputs it is
max(rtol, 1e-14) times each panel (from the calibration of the fixed rule), otherwise
the adaptive quadrature estimate.
Inputs are dimensionless and a₂₁ uses the total-projection convention.
"""
function psi30_parts(α, θ; Z, E, a21, rtol = 1e-11)
    π=typedpi(α)
    g=geometry(α, θ)
    (; r1, r2, σ, v, ξ)=g
    α=primal(g.α)<=π/2 ? g.α : π-g.α
    cl=classical(α, acos(v), Z)
    state=E*(Z*σ*(2+r1*r2)/18-(6-ξ^2)*ξ/72)-a21*(Z*σ*v/2-(6-5ξ^2)*ξ/12)
    ell, err=if Z==0
        (zero(cl), zero(primal(cl)))
    elseif primal(α) isa Float64
        elliptic_fixed(α, g, Z; rtol)
    else
        elliptic(α, g, Z; rtol)
    end
    (; classical = cl, state, elliptic = ell, value = cl+state+ell, quadrature_error = err)
end
"""Numerical ψ₃₀ at an interior angle. Uses four one-dimensional elliptic quadratures."""
psi30(α, θ; kwargs...) = psi30_parts(α, θ; kwargs...).value
include("green.jl")
include("green_k4.jl")
include("fourth_order.jl")
include("langner.jl")
export LangnerTable, psi30_langner
export psi30_green
end
