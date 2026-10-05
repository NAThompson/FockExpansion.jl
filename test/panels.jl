using Test, ForwardDiff, FockExpansion
@testset "Fixed-rule panels against adaptive BigFloat quadrature" begin
    # includes y → 0 (coalescence), where the endpoint is nearly singular,
    # and b near ±1
    for (b, y) in ((0.3, -0.6), (0.0, -0.5), (0.3, -1e-6), (0.0, -1e-6), (-0.5, -1e-8),
                   (0.999, -0.3), (-0.999, -0.01), (0.5, -1.0)),
        (c0, c2) in ((1, 0), (0, 1))
        ref = setprecision(BigFloat, 128) do
            Float64(first(FockExpansion.panel_adaptive(big(b), big(y), big(c0), big(c2);
                                                      rtol = big(1e-24))))
        end
        @test first(FockExpansion.panel(b, y, c0, c2; rtol = 1e-14)) ≈ ref rtol=2e-14
        @test first(FockExpansion.panel(b, y, c0, c2; rtol = 1e-8)) ≈ ref rtol=1e-8
    end
end
@testset "Fixed-rule panel derivatives" begin
    for (b, y) in ((0.3, -0.6), (-0.4, -0.05), (0.2, -1e-4))
        fast(v) = first(FockExpansion.panel(v[1], v[2], 0.7, -0.3; rtol = 1e-14))
        slow(v) = first(FockExpansion.panel_adaptive(v[1], v[2], 0.7, -0.3; rtol = 1e-13))
        @test ForwardDiff.gradient(fast, [b, y]) ≈ ForwardDiff.gradient(slow, [b, y]) rtol=1e-9
    end
end
