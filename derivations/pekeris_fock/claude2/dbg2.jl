include("atoms_ritz.jl")
x, y = 0.3, 0.6
println(psi30h_xy(x, y), " ", ForwardDiff.gradient(z -> psi30h_xy(z[1], z[2]), [x, y]))
h = 1e-6
println(((psi30h_xy(x + h, y) - psi30h_xy(x - h, y)) / 2h), " ", (psi30h_xy(x, y + h) - psi30h_xy(x, y - h)) / 2h)
@time for i = 1:100; ForwardDiff.gradient(z -> psi30h_xy(z[1], z[2]), [x, y + i * 1e-4]); end
