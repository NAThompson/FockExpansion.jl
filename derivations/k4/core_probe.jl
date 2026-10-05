include("core_a0.jl")
PURE[]=false
α, θ, Z = 0.7, 1.1, 2.0
x=[cos(α), sin(α)*cos(θ), sin(α)*sin(θ), 0.0]
function point(t, u, s)
    C=2t^2+u^2*(1+t^4); q1=-s*u^2*(1-t^4)/C; q2=-2t^2*(u^2-1)/C; ω=4t^2*u^2*(1+t^2)^2/C^2
    f(x1, x2)=(c=ca(x1, x2, s, Z); c[1]+c[2]*t^2)
    a0=f(0.0, 0.0); a1=ForwardDiff.derivative(z->f(z, 0.0), 0.0); a2=ForwardDiff.derivative(z->f(0.0, z), 0.0)
    a11=ForwardDiff.derivative(z->ForwardDiff.derivative(w->f(w, 0.0), z), 0.0)/2
    a22=ForwardDiff.derivative(z->ForwardDiff.derivative(w->f(0.0, w), z), 0.0)/2
    a12=ForwardDiff.derivative(z->ForwardDiff.derivative(w->f(z, w), 0.0), 0.0)
    Lplus_Q_over_l((a0, a1, a2, a11, a12, a22), q1, q2, ω, x)/C*t
end
for (t, u) in ((0.5,0.5),(0.5,1e-2),(0.5,1e-3),(0.5,1e-4),(1e-2,0.5),(1e-3,0.5),(1e-4,0.5),(0.999,0.5),(0.5,0.999))
    @printf("t=%-7g u=%-7g  %+.6e  %+.6e\n", t, u, point(t,u,1), point(t,u,-1))
end
tt, tw=gauss(48, 0.0, 1.0); uu, uw=gauss(48, 0.0, 1.0)
vals=[(abs(point(t,u,s)), t, u, s, point(t,u,s)) for t in tt for u in uu for s in (1,-1)]
sort!(vals, rev=true)
for v in vals[1:6]; @printf("t=%.6f u=%.6f s=%+d  value %+.4e\n", v[2], v[3], v[4], v[5]); end
