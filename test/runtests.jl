using Test, FockExpansion
@testset "Independent controls" begin
    α, θ = acos(0.3), acos(0.19/sqrt(0.91))
    @test psi30(α, θ; Z = 2, E = -2.9, a21 = 0.123) ≈ -1.9731181939987159649 atol=2e-11
    @test psi30(π-α, θ; Z = 2, E = -2.9, a21 = 0.123) ≈
          psi30(α, θ; Z = 2, E = -2.9, a21 = 0.123) atol=2e-11
    for α in (0.2, 0.8, 1.4), θ in (0.1, 1.2, 2.9)
        ξ = sqrt(1-sin(α)*cos(θ))
        E=-1.2
        a21=0.17
        exact = ξ*((1-2E)/24+a21/2)+ξ^3*((E-2)/72-5a21/12)
        @test psi30(α, θ; Z = 0, E, a21) ≈ exact atol=2e-13
        @test psi10(α, θ; Z = 2) ≈ ξ/2-2*(sin(α/2)+cos(α/2))
    end
    @test_throws DomainError psi30(0.0, 0.3; Z = 2, E = -2.9, a21 = 0.0)
end
@testset "Second order independent reference" begin
    for (α, θ, ref) in (
        (0.2, 0.1, 1.4354605825548435011),
        (0.8, 1.2, 2.2590497743025423171),
        (1.4, 2.9, 2.7382876747800926447),
    )
        @test psi20(α, θ; Z = 2, E = -2.9, a21 = 0.123) ≈ ref atol=3e-14
        @test psi20(π-α, θ; Z = 2, E = -2.9, a21 = 0.123) ≈ ref atol=3e-14
    end
end
@testset "Additional third order references" begin
    for (a, d, ref) in
        ((-0.27, 1.13, -1.1121070428799242549), (0.2, 1.0, -1.1432807675114344605))
        @test psi30(acos(a), acos((1-d^2)/sqrt(1-a^2)); Z = 1.7, E = -2.2, a21 = 0.31) ≈ ref atol=3e-11
    end
    for α in (0.02, 0.3, π/2, 2.2, π-0.02), θ in (0.02, 0.7, π/2, 2.6, π-0.02)
        p=psi30_parts(α, θ; Z = 2, E = -2.9, a21 = 0.123)
        @test isfinite(p.value)
        @test p.value ≈ psi30(π-α, θ; Z = 2, E = -2.9, a21 = 0.123) atol=3e-10
        @test p.value ≈ psi30(α, θ; Z = 2, E = -2.9, a21 = 0.123, rtol = 1e-12) atol=3e-10
    end
end

@testset "Elementary moment reduction: independent high-precision integrals" begin
    @test first(FockExpansion.panel(0.3, -0.6, 1, 0; rtol = 1e-12)) ≈
          0.6732045211243516227815778 atol=3e-12
    @test first(FockExpansion.panel(0.3, -0.6, 0, 1; rtol = 1e-12)) ≈
          0.3431919685359679475151188 atol=3e-12
    @test first(FockExpansion.panel(-0.4, -0.7, 1, 0; rtol = 1e-12)) ≈
          0.7305274265354698536868275 atol=3e-12
    @test first(FockExpansion.panel(-0.4, -0.7, 0, 1; rtol = 1e-12)) ≈
          0.3504151804356251522926664 atol=3e-12
    @test first(FockExpansion.panel(0, -0.5, 1, 0; rtol = 1e-12)) ≈
          0.7298385186530032081928597 atol=3e-12
    @test first(FockExpansion.panel(0, -0.5, 0, 1; rtol = 1e-12)) ≈
          0.3659786890044870941375693 atol=3e-12
    @test first(FockExpansion.panel(0.99, -0.7, 1, 0; rtol = 1e-12)) ≈
          0.6024826219923981854038521 atol=3e-12
    @test first(FockExpansion.panel(0.99, -0.7, 0, 1; rtol = 1e-12)) ≈
          0.3161481578438969642349222 atol=3e-12
    @test first(FockExpansion.panel(-0.99, -0.7, 1, 0; rtol = 1e-12)) ≈
          0.9558201890751653119529847 atol=3e-12
    @test first(FockExpansion.panel(-0.99, -0.7, 0, 1; rtol = 1e-12)) ≈
          0.3798778625047066710252676 atol=3e-12
    @test first(FockExpansion.panel(0, -0.01, 1, 0; rtol = 1e-12)) ≈
          0.911000717823391833817099 atol=3e-12
    @test first(FockExpansion.panel(0, -0.01, 0, 1; rtol = 1e-12)) ≈
          0.4674961501884395226930881 atol=3e-12
    @test first(FockExpansion.panel(0, -0.999, 1, 0; rtol = 1e-12)) ≈
          0.6170337679986385553232287 atol=3e-12
    @test first(FockExpansion.panel(0, -0.999, 0, 1; rtol = 1e-12)) ≈
          0.3046470886379322523119919 atol=3e-12
end

include("forwarddiff.jl")

include("panels.jl")

include("alternative_representations.jl")
include("fourth_order.jl")
