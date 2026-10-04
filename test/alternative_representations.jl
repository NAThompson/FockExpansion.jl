using Test, FockExpansion
@testset "Full Green recurrence inversion" begin
    @test psi30_green(
        acos(0.3),
        acos(0.19/sqrt(0.91));
        Z = 2.0,
        E = -2.9,
        a21 = 0.123,
        rtol = 1e-7,
    ) ≈ -1.973118193998716 atol=2e-7
    @test psi30_green(0.8, 1.2; Z = 0.0, E = -1.2, a21 = 0.17, rtol = 1e-7) ≈
          psi30(0.8, 1.2; Z = 0.0, E = -1.2, a21 = 0.17) atol=2e-7
end
@testset "Independent corrected Langner series" begin
    for (Z, E, a21, a, d, ref) in (
        (2.0, -2.9, 0.123, 0.3, 0.9, -1.9731181939987159649),
        (1.7, -2.2, 0.31, -0.27, 1.13, -1.1121070428799242549),
        (1.7, -2.2, 0.31, 0.2, 1.0, -1.1432807675114344605),
    )
        table=LangnerTable(; Z)
        α=acos(a)
        θ=acos((1-d*d)/sqrt(1-a*a))
        for method in (:powers, :horner, :direct)
            @test psi30_langner(α, θ, table; E, a21, method) ≈ ref atol=3e-13
        end
        @test psi30_green(α, θ; Z, E, a21, rtol = 1e-7) ≈ ref atol=2e-7
    end
    @test_throws ArgumentError psi30_langner(
        1.0,
        1.0,
        LangnerTable();
        E = -2.9,
        a21 = 0.1,
        N = 61,
    )
end
@testset "Float64 Clausen against arbitrary precision" begin
    for x in range(-4π, 4π, length = 201)
        ref=setprecision(192) do ;
            Float64(clausen2(BigFloat(x)))
        end
        @test clausen2(x) ≈ ref atol=3e-15
    end
end
@testset "Endpoint source limits" begin
    for x in (0.3, 1.2), μ in (-1.0, 1.0)
        @test isfinite(FockExpansion.chi(x, μ))
        @test FockExpansion.chi(x, μ) ≈
              FockExpansion.chi(x, μ*(1-1e-12)) atol=3e-10
    end
end
