# Azimuthal average of the k = 4 Green's function in one variable s ∈ (2, ∞):
# with C = a + d cos φ, (π-γ)/sin γ = ∫₀^∞ dt/(t² - 2tC + 1), the φ-average of
# -(π-γ)T₃(C)/sin γ becomes, after s = t + 1/t,
#   -2∫₂^∞ ds/√(s²-4) [ (s³-3s)/(2√((s-2c₁)(s-2c₂))) - (4a²+2d²-3 + 2as + s²)/2 ],  c₁,₂ = a ± d.
using QuadGK, Printf
G(γ)=((γ-π)*cos(3γ)+sin(3γ)/6)/(4π^2*sin(γ))
for (a, d) in ((0.2, 0.3), (-0.5, 0.4), (0.6, 0.35))
    direct=quadgk(φ->G(acos(a+d*cos(φ))), 0, π; rtol=1e-13)[1]/π      # φ-average
    c1, c2=a+d, a-d
    # the two terms cancel to O(1/s) at large s: evaluate in BigFloat
    J=setprecision(256) do
        A, D=big(a), big(d); C1, C2=A+D, A-D
        f(s)=(s^3-3s)/(2sqrt((s-2C1)*(s-2C2)))-(4A^2+2D^2-3+2A*s+s^2)/2
        Float64(-2quadgk(v->f(2cosh(v)), big(0), big(60); rtol=big(1e-20))[1])   # s = 2cosh v
    end
    U2avg=4(a^2+d^2/2)-1
    @printf("a=%.1f d=%.2f: direct %.13f   s-form %.13f\n", a, d, direct, (J+U2avg/6)/(4π^2))
end
