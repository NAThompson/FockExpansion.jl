# First piece of the ψ₄₀ core: L₄⁺[-2V₀ F_a], V₀ = 1/ξ, with F_a the b = ±cos α panels of ψ₃₀.
# With atan(a/b)/a = ∫₀¹ b du/(b² + a²u²), each panel is
#   c₀N₀ + c₂N₂ = (2/π) ∫₀¹dt ∫₀¹du L(t) (c₀ + c₂t²) (-2yt) / D,   D = 4y²t² + P u²,
# L(t) = log((1+t)/(1-t)), y = -ξ/√2 (so -2yt = √2 ξ t), and D is linear in x:
#   D = C (1 - q·x),  C = 2t² + u²(1+t⁴),  q = (∓u²(1-t⁴), 2t²(1-u²))/C   for b = ±x₁.
# The 1/ξ cancels the ξ, leaving Q(x)/(1 - q·x) with Q quadratic: zonal solutions about q.
include("feynman41.jl")
const SQ2=sqrt(2.0)
# coefficient polynomials of F_a in x (b = s·x₁, y² = (1-x₂)/2)
function ca(x1, x2, s, Z)
    b=s*x1; y2=(1-x2)/2
    c0=-SQ2*Z*(9Z*b^2-18Z-8b^2+8)/576
    c2=SQ2*Z*(b-1)*(90Z*y2-45Z-18b-32y2+16)/288
    c0, c2
end
# L₄⁺[Q(x)/(1-q·x)] for Q given as a function of (x₁,x₂): expand Q to its monomials by
# exact differentiation at 0, then use x_i/ℓ = -∂ᵢlog ℓ, x_i x_j/ℓ = ∂ᵢ∂ⱼ(ℓ log ℓ).
function Lplus_Q_over_l(Qcoef, q1, q2, ω, x)
    # Qcoef = (q0, q1x, q2x, q11, q12, q22) for q0 + a x₁ + b x₂ + c x₁² + d x₁x₂ + e x₂²
    a0, a1, a2, a11, a12, a22=Qcoef
    # differentiate in δ = q - q₀ with ω = ω₀ - 2q₀·δ - |δ|², ω₀ exact: no cancellation near |q| = 1
    F(v, kind)=Uf(kind, q1+v[1], q2+v[2], ω-2(q1*v[1]+q2*v[2])-v[1]^2-v[2]^2, x)
    v0=[0.0, 0.0]
    g=ForwardDiff.gradient(v->F(v, :log), v0)
    H=ForwardDiff.hessian(v->F(v, :xlog), v0)
    a0*F(v0, :inv)-a1*g[1]-a2*g[2]+a11*H[1, 1]+a12*H[1, 2]+a22*H[2, 2]
end
function graded_rule(n, levels; q=0.25)
    pts=sort(unique(vcat([0.5*q^k for k in 0:levels], [1-0.5*q^k for k in 0:levels], [0.0, 0.5, 1.0])))
    x, w=gauss(n); xs=Float64[]; ws=Float64[]
    for i in 1:length(pts)-1
        h=(pts[i+1]-pts[i])/2; m=(pts[i+1]+pts[i])/2
        append!(xs, m.+h.*x); append!(ws, h.*w)
    end
    xs, ws
end
function core_a0(α, θ, Z; nt=24, nu=24, levels=0)
    x=[cos(α), sin(α)*cos(θ), sin(α)*sin(θ), 0.0]
    tt, tw=levels==0 ? gauss(nt, 0.0, 1.0) : graded_rule(nt, levels)
    uu, uw=levels==0 ? gauss(nu, 0.0, 1.0) : graded_rule(nu, levels)
    total=0.0
    for (t, wt) in zip(tt, tw), (u, wu) in zip(uu, uw), s in (1, -1)
        C=2t^2+u^2*(1+t^4)
        q1=-s*u^2*(1-t^4)/C; q2=-2t^2*(u^2-1)/C
        ω=4t^2*u^2*(1+t^2)^2/C^2
        # monomial coefficients of (c₀ + c₂t²) as a polynomial in x
        f(x1, x2)=(c=ca(x1, x2, s, Z); c[1]+c[2]*t^2)
        a0=f(0.0, 0.0)
        a1=ForwardDiff.derivative(z->f(z, 0.0), 0.0); a2=ForwardDiff.derivative(z->f(0.0, z), 0.0)
        a11=ForwardDiff.derivative(z->ForwardDiff.derivative(w->f(w, 0.0), z), 0.0)/2
        a22=ForwardDiff.derivative(z->ForwardDiff.derivative(w->f(0.0, w), z), 0.0)/2
        a12=ForwardDiff.derivative(z->ForwardDiff.derivative(w->f(z, w), 0.0), 0.0)
        val=Lplus_Q_over_l((a0, a1, a2, a11, a12, a22), q1, q2, ω, x)/C
        total+=wt*wu*log((1+t)/(1-t))*SQ2*t*val
    end
    -2*(2/π)*total
end
