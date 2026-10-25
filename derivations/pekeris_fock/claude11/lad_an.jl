include(joinpath(@__DIR__, "logcoef.jl"))
setprecision(BigFloat, 256); T = BigFloat
R80 = deserialize("ladder_80.jls"); R64 = deserialize("ladder_64.jls")
b80 = R80.out[(0, 0)][(4, 0)]; b64 = R64.out[(0, 0)][(4, 0)]
c80 = chebcoef(b80); c64 = chebcoef(b64)
println("conv: ", Float64(maximum(abs(chebeval(c80, T(x)) - chebeval(c64, T(x))) for x in (-0.9, -0.3, 0.2, 0.7))), " tail ", Float64(abs(c80[end])))
serialize("ladder_sectors.jls", Dict((0, (4, 0), 0) => c80))
# trig fit on the far edge: alpha = 5pi/4 + pi ty/4
P = T(π)
tys = [T(-0.95) + T(1.9) * (i - 1) / 59 for i = 1:60]
for K in (2, 4, 6, 8)
    fs = vcat([x -> one(T)], [x -> cos(m * x) for m = 1:K], [x -> sin(m * x) for m = 1:K])
    A = [f(5P / 4 + P * ty / 4) for ty in tys, f in fs]
    b = [chebeval(c80, ty) for ty in tys]
    x = A \ b
    println("K=$K trig fit residual ", Float64(maximum(abs.(A * x - b))))
end
fs = [x -> one(T), x -> cos(2x), x -> sin(2x), x -> cos(x), x -> sin(x)]
A = [f(5P / 4 + P * ty / 4) for ty in tys, f in fs]
b = [chebeval(c80, ty) for ty in tys]
x = A \ b
println("c0, cos2a, sin2a, cos a, sin a: ", [string(round(xx, sigdigits = 40)) for xx in x])
println("times pi: ", [string(round(xx * P, sigdigits = 30)) for xx in x[1:3]], " times pi^2: ", [string(round(xx * P^2, sigdigits = 30)) for xx in x[1:3]])
