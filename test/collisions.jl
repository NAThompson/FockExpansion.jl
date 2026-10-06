using Test, FockExpansion
@testset "psi30 near the collisions" begin
    # 192-bit references. Near α = 0, π the Float64 panels degenerate (cos α rounds to
    # ±1) and psi30 extrapolates in α; near the e-e coalescence b = -v → ∓1.
    for (α, θ, ref) in (
        (1.0e-30, 0.7, -0.34315540728835314326),
        (1.0e-9, 1.0e-9, -0.34315540796464227580),
        (1.0e-6, 3.0, -0.34315602283751334072),
        (3.141592653588793, 2.0, -0.34315540728898634136),
        (1.5707963267948966, 1.0e-9, -1.9972689383294891190),
        (1.5707963257948965, 1.0e-9, -1.9972689379313956659),
        (1.5707963267948966, 3.141592652589793, -1.1948607235122435170),
        (1.2, 1.0e-12, -1.5958232974746062353),
    )
        @test psi30(α, θ; Z = 1.7, E = -2.2, a21 = 0.31) ≈ ref atol=1e-11
    end
end
