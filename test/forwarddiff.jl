using Test, ForwardDiff, FockExpansion
@testset "ForwardDiff recurrence and precision" begin
    function residual(x; rtol)
        α, θ=x
        Z=oftype(α, 2)
        E=parse(typeof(α), "-2.9037243770341195983111592451944044467")
        a21=parse(typeof(α), "0.47674787900")
        f(t) = psi30(t[1], t[2]; Z, E, a21, rtol)
        v=f(x)
        g=ForwardDiff.gradient(f, x)
        h=ForwardDiff.hessian(f, x)
        lhs=-4*(h[1, 1]+2cos(α)/sin(α)*g[1]+(h[2, 2]+cos(θ)/sin(θ)*g[2])/sin(α)^2)-21v
        V=1/sqrt(1-sin(α)*cos(θ))-Z*(1/cos(α/2)+1/sin(α/2))
        rhs=10psi31(α, θ; Z)-2V*psi20(α, θ; Z, E, a21)+2E*psi10(α, θ; Z)
        lhs-rhs
    end
    @test abs(residual([1.1, 1.2]; rtol = 1e-11)) < 1e-8
    setprecision(128) do
        @test psi30(big"1.1", big"1.2"; Z = big"2", E = big"-2.9", a21 = big".123") isa
              BigFloat
        @test abs(residual(BigFloat[big"1.1", big"1.2"]; rtol = big"1e-24")) < big"1e-22"
        @test abs(residual(BigFloat[BigFloat(π)/2, big"1.2"]; rtol = big"1e-24")) <
              big"1e-22"
        @test abs(residual(BigFloat[BigFloat(π)/2, BigFloat(π)/2]; rtol = big"1e-24")) <
              big"1e-22"
    end
end
@testset "Removable primitive derivatives" begin
    @test ForwardDiff.derivative(FockExpansion.L, 0.0) == 1
    @test ForwardDiff.derivative(x->ForwardDiff.derivative(FockExpansion.L, x), 0.0) == 0
    @test ForwardDiff.derivative(FockExpansion.T, 0.0) == 0
end
