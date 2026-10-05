using Test, FockExpansion
const F4 = FockExpansion

# ∫ f sin²α sinθ dα dθ, graded toward the Coulomb singularities at α = 0, π/2, π and θ = 0.
function sphere_integral(f; n = 16, levels = 8)
    rule = F4.gauss(n)
    A = sort(unique(vcat(F4.graded_points(0.0, π/2, true, true; levels),
                         F4.graded_points(π/2, Float64(π), true, true; levels))))
    T = F4.graded_points(0.0, Float64(π), true, false; levels)
    F4.rule_sum(a -> sin(a)^2*F4.rule_sum(t -> f(a, t)*sin(t), T, rule), A, rule)
end
function lambda2(f, a, t; h = 2e-3)
    f0 = f(a, t)
    faa = (f(a+h, t)-2*f0+f(a-h, t))/h^2
    fa = (f(a+h, t)-f(a-h, t))/(2h)
    ftt = (f(a, t+h)-2*f0+f(a, t-h))/h^2
    ft = (f(a, t+h)-f(a, t-h))/(2h)
    -4*(faa+2cot(a)*fa+(ftt+cot(t)*ft)/sin(a)^2)
end

@testset "Fourth-order logarithmic coefficients" begin
    B = (π-2)/(3π)
    ξ(a, t) = F4.xi_stable(a, t)
    E = -2.9
    # Green's function solver against the closed-form Z¹ and Z³ components.
    h1(a, t) = (π-2)/(18π)*(6*(1-2E)+(12*E-5)*ξ(a, t)^2)
    p1(a, t) = (π-2)/(2880π)*(3*(32*E-15)-8*(12*E-5)*ξ(a, t)^2)
    h3(a, t) = B*(2+tan(a/2)+cot(a/2))*sin(a)*cos(t)
    p3(a, t) = -(π-2)/(120π)*(4+5sin(a))*sin(a)*cos(t)
    for (a, t) in ((0.7, 1.1), (2.6, 2.5))
        @test F4.solve_k4(h1, a, t; n = 8, levels = 8) ≈ p1(a, t) atol=1e-13
        @test F4.solve_k4(h3, a, t; n = 8, levels = 8) ≈ p3(a, t) atol=1e-13
    end
    # ψ₄₂ is a k = 4 harmonic.
    @test lambda2((a, t) -> psi42(a, t; Z = 2.0), 0.9, 1.3)-32psi42(0.9, 1.3; Z = 2.0) ≈ 0 atol=1e-6
    # Z² component: recurrence by finite differences, and exchange symmetry.
    function h2(a, t)
        ς = cos(a/2)+sin(a/2)
        V1 = 1/sin(a/2)+1/cos(a/2)
        24psi42(a, t; Z = 1.0)-2*B*ς*sin(a)*cos(t)/2/ξ(a, t)+2V1*B*ξ(a, t)*(5ξ(a, t)^2-6)/12
    end
    z2(a, t) = F4.psi41_z2(a, t)
    @test lambda2(z2, 0.7, 1.1)-32z2(0.7, 1.1) ≈ h2(0.7, 1.1) atol=1e-6
    @test F4.psi41_z2(π-0.7, 1.1) ≈ F4.psi41_z2(0.7, 1.1) atol=1e-15
    # The one-dimensional representation against the Green's function.
    for (a, t) in ((0.7, 1.1), (1.5, 0.3), (0.05, 1.5))
        @test F4.psi41_z2(a, t) ≈ F4.psi41_z2_green(a, t; n = 10, levels = 10) atol=1e-13
    end
    # Harmonic part: solvability of the ψ₄₀ equation, ⟨Y₄ₗ, 12ψ₄₁+2ψ₄₂-2Vψ₃₀+2Eψ₂₀⟩ = 0.
    # The pure parts of ψ₄₁ are orthogonal to Y₄ₗ, so only its harmonic part enters.
    Z, E, a21 = 1.7, -2.2, 0.31
    for Y in (F4.harmonic40, F4.harmonic42)
        f(a, t) = Y(a, t)*(12F4.psi41_harmonic(a, t; Z, E, a21)+2psi42(a, t; Z)-
                           2*(1/ξ(a, t)-Z*(1/sin(a/2)+1/cos(a/2)))*psi30(a, t; Z, E, a21, rtol = 1e-12)+
                           2E*psi20(a, t; Z, E, a21))
        @test sphere_integral(f) ≈ 0 atol=1e-10
    end
    @test psi41(π-0.7, 1.1; Z, E, a21) ≈ psi41(0.7, 1.1; Z, E, a21) atol=1e-14
    # psi41 assembles its components with the closed forms above.
    a, t = 0.7, 1.1
    z1 = (π-2)/(2880π)*(3*(32*E-15)-8*(12*E-5)*ξ(a, t)^2)
    @test psi41(a, t; Z, E, a21)-Z^2*F4.psi41_z2(a, t)-F4.psi41_harmonic(a, t; Z, E, a21) ≈
          Z*z1+Z^3*p3(a, t) atol=1e-15
    # ψ₄₀ at Z = 0 against its closed form, and exchange symmetry.
    let E = -2.9, a21 = 0.4
        p40_0(a, t) = (12*E^2-11*E+6*a21+1)/1152-(72*E*a21+E-30*a21-2)*sin(a)*cos(t)/720
        @test psi40(0.7, 1.1; Z = 0.0, E, a21) ≈ p40_0(0.7, 1.1) atol=1e-13
        @test psi40(π-0.7, 1.1; Z = 2.0, E, a21, step = 0.125) ≈ psi40(0.7, 1.1; Z = 2.0, E, a21, step = 0.125) atol=1e-14
    end
end
