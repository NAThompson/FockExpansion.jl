setprecision(BigFloat, 256)
include("atoms7.jl")
PSI30_DIGITS[] = 45
for (x, y) in ((big"0.3", big"0.6"), (big"0.1", big"0.9"), (big"0.45", big"0.5"))
    p = pieces_xy(x, y)
    st = begin
        α, μ = angles(x, y); θ = theta(x, y); g = FockExpansion.geometry(α, θ)
        sqrt((1 + x^2) / 2)^3 * FockExpansion.psi30_parts(α, θ; Z = big(2), E = E0, a21 = A21, rtol = big(10)^-45).state
    end
    println("sum of pieces + state - psi30 atom = ", Float64(sum(p) + st - psi30h_big(x, y)))
end
