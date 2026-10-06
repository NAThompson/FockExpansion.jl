# Size of the Fock coefficients ψ_kp as k grows (helium, E = exact, free constants 0).
include("../k4/fock_spectral.jl")
T=length(ARGS)>0 ? BigFloat : Float64
setprecision(BigFloat, 256)
n=parse(Int, get(ARGS, 2, "48")); kmax=parse(Int, get(ARGS, 3, "24"))
g=SGrid(T, n)
t=@elapsed ψ=sfock(g; Z=2, E=T(-2.903724377034119598311159), kmax)
@printf("n=%d %s, k ≤ %d in %.1f s\n", n, T, kmax, t)

for k in 0:kmax
    m=[Float64(maximum(abs, ψ[(k, p)])) for p in 0:k÷2]
    @printf("k=%2d  max_p |ψ_kp| = %.3e (p=%d)   |ψ_k0| = %.3e\n", k, maximum(m), argmax(m)-1, m[1])
end
