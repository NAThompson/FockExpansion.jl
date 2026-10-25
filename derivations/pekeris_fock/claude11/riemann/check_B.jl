# Check B_{A+B}(Y) = int int_{z0>1} T3(z0)/(2pi sqrt(z0^2-1)) R(X,eta) dX deta against the chain's far-edge log coefficients.
include(joinpath(@__DIR__, "..", "chain.jl"))   # Elem, lobatto, SGrid ... (main guarded)
using QuadGK, Serialization
setprecision(BigFloat, 256)
n = 64; N = 48; P = Float64(π); a = P / 2; b = 3P / 2; Xm = P
T = BigFloat
eA = Elem(T, T(π) / 2, T(π), N); eB = Elem(T, 3 * T(π) / 2, T(π), N)
st = Dict(L => deserialize(joinpath(@__DIR__, "..", "ck", "lad$(L)64_step_Z0_s0.jls")) for L in ("A", "B"))
VS(kp) = Float64.(st["A"][1][kp] + st["B"][1][kp])
VF(kp) = (fa = st["A"][2][kp]; fb = st["B"][2][kp]; [Float64.((m <= length(fa) ? fa[m] : zero(fa[1])) + (m <= length(fb) ? fb[m] : zero(fb[1]))) for m = 1:max(length(fa), length(fb))])
Bchain = Float64.(st["A"][3][(4, 0)] + st["B"][3][(4, 0)])
xg = [cos(P * (j + 0.5) / n) for j = 0:n-1]; bwg = [(-1)^j * sin(P * (j + 0.5) / n) for j = 0:n-1]
Ygrid = P .+ P / 2 .* xg
function ywts(η); x = 2(η - P) / P; w = bwg ./ (x .- xg); w ./ sum(w); end
function lobw(e, X)   # barycentric weights on Lobatto nodes of element e at X
    Xn = Float64.(e.X); Nn = length(Xn)
    c = [(j == 1 || j == Nn ? 0.5 : 1.0) * (-1)^(j - 1) for j = 1:Nn]
    i = findfirst(==(X), Xn); i === nothing || return [k == i ? 1.0 : 0.0 for k = 1:Nn]
    w = c ./ (X .- Xn); w ./ sum(w)
end
S3 = VS((3, 0)); S41 = VS((4, 1)); F3 = VF((3, 0)); F41 = VF((4, 1))
Lc(X) = -log(-cos(X))   # chain's log variable
function fval(X, η, Smat, Flist)
    wy = ywts(η)
    if X <= Xm
        return dot(lobw(eA, X), Smat * wy)
    else
        wx = lobw(eB, X); L = Lc(X)
        return sum(dot(wx, Flist[m] * wy) * L^(m - 1) for m in eachindex(Flist))
    end
end
function R(X, η)
    α, β = (X + η) / 2, (η - X) / 2
    w = 2sin(α) * sin(β); wv0 = 2sqrt(2) * sin(α) * cos(β / 2)
    -(3 / 4) * w * fval(X, η, S41, F41) + (wv0 / 8) * fval(X, η, S3, F3)
end
z0(Y, X, η) = cos((b + Y) / 2) * cos((X + η) / 2) + cos((Y - b) / 2) * cos((η - X) / 2)
function inner(Y, X)
    es = range(a, b; length = 801); f(e) = z0(Y, X, e) - 1
    cuts = Float64[]; f(es[1]) > 0 && push!(cuts, a)
    for i = 1:800
        f(es[i]) * f(es[i+1]) < 0 || continue
        lo, hi = es[i], es[i+1]
        for _ = 1:80; m = (lo + hi) / 2; f(lo) * f(m) <= 0 ? (hi = m) : (lo = m); end
        push!(cuts, (lo + hi) / 2)
    end
    f(es[end]) > 0 && push!(cuts, b)
    tot = 0.0
    for k = 1:2:length(cuts)-1
        e1, e2 = cuts[k], cuts[k+1]
        g(θ) = (η = e1 + (e2 - e1) * (1 - cos(θ)) / 2; z = z0(Y, X, η); z <= 1 ? 0.0 : (4z^3 - 3z) / (2P * sqrt(z^2 - 1)) * R(X, η) * (e2 - e1) * sin(θ) / 2)
        tot += quadgk(g, 0, P; rtol = 1e-10)[1]
    end
    tot
end
for j in (10, 25, 40, 55)
    Y = Ygrid[j]
    Bk = quadgk(X -> inner(Y, X), a, Xm, b; rtol = 1e-9)[1]
    println("Y=", round(Y, digits = 4), "  kernel B = ", Bk, "   chain (-L coeff) = ", -Bchain[j], "   chain L coeff = ", Bchain[j])
end
