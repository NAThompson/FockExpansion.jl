setprecision(BigFloat, 256)
include("atoms3.jl")
using Random
Random.seed!(1)
worst = Dict{Symbol,Float64}()
for i = 1:30
    y = BigFloat(rand()) * big"0.98" + big"0.01"; x = y * BigFloat(rand()) * big"0.98"
    worst[:chi] = max(get(worst, :chi, 0.0), Float64(abs(chi_codex(x, y) - chi_xy(x, y))))
    # psi30: 256-bit/45 digits vs 192-bit/35 digits, and vs Float64
    PSI30_DIGITS[] = 45; v1 = psi30h_big(x, y); g1 = ForwardDiff.gradient(z -> psi30h_big(z[1], z[2]), [x, y])
    v2 = setprecision(BigFloat, 192) do
        PSI30_DIGITS[] = 35; psi30h_big(BigFloat(x), BigFloat(y))
    end
    worst[:psi_prec] = max(get(worst, :psi_prec, 0.0), Float64(abs(v1 - v2)))
    worst[:psi_f64] = max(get(worst, :psi_f64, 0.0), Float64(abs(v1 - psi30h_xy(Float64(x), Float64(y)))))
    h = big"1e-30"
    fd = (psi30h_big(x + h, y) - psi30h_big(x - h, y)) / 2h
    worst[:psi_grad] = max(get(worst, :psi_grad, 0.0), Float64(abs(fd - g1[1])))
end
println(worst)
