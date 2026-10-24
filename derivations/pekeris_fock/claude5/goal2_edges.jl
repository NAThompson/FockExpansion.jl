# Goal 2 admissibility probe: behaviour of the psi30 atom (h_x^3 psi30) and its derivatives near the corners/edges of the
# angular triangle 0 <= |x| <= y <= 1 (corner (0,0): r12 -> 0; corner (1,1): r2 -> 0; edges y=x and y=1: collinear faces).
setprecision(BigFloat, 192)
include("../claude3/atoms3.jl")
PSI30_DIGITS[] = 35
using Printf
f(z) = psi30h_big(z[1], z[2])
H(z) = ForwardDiff.hessian(f, z)
for (name, path) in (("e-e corner (0,0) along x=0.5y", t -> [t / 2, t]),
                     ("nucleus corner (1,1) along x = 1-2t, y = 1-t", t -> [1 - 2t, 1 - t]),
                     ("collinear edge y=x, at x=0.4: y = 0.4+t", t -> [big"0.4", big"0.4" + t]),
                     ("collinear edge y=1, at x=0.3: y = 1-t", t -> [big"0.3", 1 - t]))
    println(name)
    for t in (big"1e-2", big"1e-3", big"1e-4", big"1e-5")
        z = path(t)
        v = f(z); g = ForwardDiff.gradient(f, z); h = H(z)
        @printf("   t=%.0e  value=% .6e  |grad|=%.4e  |hess|=%.4e\n", Float64(t), Float64(v), Float64(norm(g)), Float64(norm(h)))
    end
end
