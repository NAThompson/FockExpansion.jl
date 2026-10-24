# Grid convergence of high-order ψ_kp: interpolated values at off-grid points for several n.
include("../k4/fock_spectral.jl")
kmax=parse(Int, get(ARGS, 1, "24")); ns=parse.(Int, split(get(ARGS, 2, "48,64,80,96"), ","))
pts=((0.7, 1.1), (1.5, 0.3), (2.3, 2.0), (1.2, 2.9))
vals=Dict{Int,Any}()
for n in ns
    g=SGrid(Float64, n)
    t=@elapsed ψ=sfock(g; Z=2, E=-2.903724377034119598311159, kmax)
    vals[n]=Dict(k=>[interp(g, ψ[(k, p)], a, θ) for p in 0:k÷2, (a, θ) in pts] for k in 0:kmax)
    @printf("n=%d done in %.0f s\n", n, t)
end
nref=ns[end]
for k in 0:2:kmax
    sc=maximum(abs, vals[nref][k])
    @printf("k=%2d scale %.2e  rel err:", k, sc)
    for n in ns[1:end-1]; @printf("  n=%d %.1e", n, maximum(abs, vals[n][k]-vals[nref][k])/sc); end
    println()
end
