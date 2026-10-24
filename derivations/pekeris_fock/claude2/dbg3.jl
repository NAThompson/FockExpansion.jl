include("atoms_ritz.jl")
setprecision(BigFloat, 128)
function psi30h_big(x, y)
    α, μ = angles(x, y)
    hx = sqrt((1 + x^2) / 2)
    hx^3 * psi30(α, acos(μ); Z = big(2.0), E = E0, a21 = A21, rtol = big(1e-20))
end
x, y = big"0.3", big"0.3" + big"1e-12"
@time g = ForwardDiff.gradient(z -> psi30h_big(z[1], z[2]), [x, y])
println(Float64.(g))
@time g = ForwardDiff.gradient(z -> psi30h_big(z[1], z[2]), [x, y])
