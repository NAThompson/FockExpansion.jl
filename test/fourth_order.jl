using Test, FockExpansion, ForwardDiff
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
    # Closed form against the Feynman integral, which is accurate away from the boundaries;
    # nearer the coalescences, Float64 against the same closed form in BigFloat.
    for a in (0.1, 0.7, 1.2, π/2, 2.3, 2.9), t in (0.05, 0.8, 1.6, 2.4, 2.8)
        @test F4.psi41_z2(a, t) ≈ F4.psi41_z2_feynman(a, t) atol=2e-14
    end
    setprecision(256) do
        for (a, t) in ((1e-8, 1.0), (0.3, 1e-9), (2.5, π-1e-7), (π/2+1e-7, 1e-7), (π/2-1e-9, 0.3), (π-1e-6, 2.0))
            @test F4.psi41_z2(a, t) ≈ F4.psi41_z2(big(a), big(t)) atol=2e-16
        end
        # 60-digit spectral values (derivations/pekeris_fock/claude9/spec41.txt), which include
        # the Z² part of the hierarchy's Y₄ₗ component; inputs and values truncated to 44 digits.
        c40, c42 = big"0.00305932040176846773106267150356500648940567", big"-0.00536795642440787076597402151228072313907757"
        for (a, t, r) in ((big"1.82372502624377676315790053928930750097518991", big"2.24353309293961854956350058611345611763909834",
                           big"-0.016804267094946716128345686750419657949132468"),
                          (big"0.08367293768469548707429164363234652793861687", big"1.15685392202119356462696844442813816961611090",
                           big"0.0097401818358765451959522972472719742714383846"),
                          (big"0.50332423416928334559178192157255312865667850", big"1.89434614490375545989111083635481249789945774",
                           big"0.0044515806202731535659140238620186541218321096"))
            @test F4.psi41_z2(a, t)+c40*F4.harmonic40(a, t)+c42*F4.harmonic42(a, t) ≈ r atol=1e-40
        end
    end
    # The recurrence with ForwardDiff derivatives instead of finite differences.
    let a = 0.7, t = 1.1, f = x -> F4.psi41_z2(x[1], x[2])
        g, H = ForwardDiff.gradient(f, [a, t]), ForwardDiff.hessian(f, [a, t])
        Λ2 = -4*(H[1, 1]+2cot(a)*g[1]+(H[2, 2]+cot(t)*g[2])/sin(a)^2)
        @test Λ2-32f([a, t]) ≈ h2(a, t) atol=1e-13
        @test ForwardDiff.derivative(b -> psi41(b, t; Z = 2.0, E, a21 = 0.4), a) ≈
              (psi41(a+1e-5, t; Z = 2.0, E, a21 = 0.4)-psi41(a-1e-5, t; Z = 2.0, E, a21 = 0.4))/2e-5 atol=1e-9
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
    # Independent reference values from the 256-bit spectral solution of the whole hierarchy
    # (derivations/k4/fock_spectral_ref.jl; grids n = 48 and 64 agree to 1e-38).
    let Z = 1.7, E = -2.2, a21 = 0.31
        for (a, t, r41, r40) in ((0.7, 1.1, -0.01214924667877228035760333, 0.6827950979687855037644383),
                                 (1.5, 0.3, -0.1603145467166121681419629, 0.4233893486868100167537526),
                                 (2.3, 2.0, 0.04158027503238696413555123, 0.6286370389039812304249224))
            @test psi41(a, t; Z, E, a21) ≈ r41 atol=2e-14
            @test psi40(a, t; Z, E, a21, a40 = 0.1, a42 = -0.2) ≈ r40 atol=5e-14
        end
    end
    # Near and at the poles θ = π, 0, against the spectral solution (128-bit, n = 48, 64).
    let Z = 1.7, E = -2.2, a21 = 0.31
        for (t, r40) in ((π-0.01, 0.4479160029353247366910543), (π, 0.4479070511588533572141819),
                         (0.0, 0.5968282234720722454807952))
            @test psi40(1.0, t; Z, E, a21) ≈ r40 atol=5e-14
        end
    end
    # Near α = π/2, where an 8-point G₂ azimuthal rule gave errors up to 8e-13 (60-digit spectral
    # data at E = a₂₁ = 0, derivations/pekeris_fock/claude9/spec40.txt).
    for (a, t, r40) in ((1.5363787940726878, 1.8064754809386592, 0.3722263813907205728692752),
                        (1.7299206379108274, 1.7460469161933476, 0.3622273973978488858445726),
                        (1.6301394819622352, 1.6946884488181788, 0.3527984914123717304411289))
        @test psi40(a, t; Z = 2.0, E = 0.0, a21 = 0.0) ≈ r40 atol=3e-14
    end
    # Float64-only functions refuse wider inputs instead of dropping precision.
    @test_throws ArgumentError psi41(big"0.7", 1.1; Z = 2.0, E, a21 = 0.4)
    @test_throws ArgumentError psi40(0.7, 1.1; Z = big"2.0", E, a21 = 0.4)
    @test_throws ArgumentError F4.psi41_z2_feynman(big"0.7", 1.1)
end
