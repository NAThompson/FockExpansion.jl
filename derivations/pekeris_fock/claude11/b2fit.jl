# Fit the double variation B2(alpha) = c0 + a cos(2 alpha) + b sin(2 alpha) on arc a (alpha = 3pi/4 + pi ty/4).
include(joinpath(@__DIR__, "logcoef.jl"))
setprecision(BigFloat, 256); T = BigFloat
L = deserialize("chain_sectors.jls")
c = L[(1, (4, 0), 2)]
P = T(π)
tys = [T(-0.95) + T(1.9) * (i - 1) / 39 for i = 1:40]
A = [f(3P / 4 + P * ty / 4) for ty in tys, f in (x -> one(T), x -> cos(2x), x -> sin(2x))]
b = [chebeval(c, ty) for ty in tys]
x = A \ b
println("c0, a, b = ", x)
println("max residual = ", Float64(maximum(abs.(A * x - b))))
