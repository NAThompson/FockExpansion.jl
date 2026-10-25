# Z-sector log coefficients A_m(X0, Y) of the continued psi_kp at the far edge X0 (from extend.jl), as Chebyshev series in ty.
include(joinpath(@__DIR__, "extend.jl"))
using Serialization, LinearAlgebra, Printf

function sector_logcoefs(R; ks = [(2, 0), (3, 0), (3, 1), (4, 0), (4, 1)])
    T = BigFloat
    Zs = sort(collect(keys(R.results)))
    nz = length(Zs)
    V = inv([T(z)^m for z in Zs, m = 0:nz-1])
    out = Dict{Any,Vector{T}}()
    for key in ks, m = 1:3
        vals = [begin
            v = R.results[z].ψ1[key]
            m + 1 <= length(v) ? v[m+1][R.e1.iA, :] : zeros(T, R.g.n)
        end for z in Zs]
        for j = 0:nz-1
            s = sum(V[j+1, i] * vals[i] for i = 1:nz)
            maximum(abs.(s)) > T(10)^-40 && (out[(key, m, j)] = s)
        end
    end
    out
end
# Chebyshev coefficients of values on the Gauss grid in ty (x_j = cos(pi (j+1/2)/n))
function chebcoef(h::Vector{T}) where {T}
    n = length(h); P = T(π)
    [2 / T(n) * sum(h[j+1] * cos(P * k * (j + T(1) / 2) / n) for j = 0:n-1) / (k == 0 ? 2 : 1) for k = 0:n-1]
end
chebeval(c, x) = (b1 = zero(x); b2 = zero(x); for k = length(c)-1:-1:1; b1, b2 = 2x * b1 - b2 + c[k+1], b1; end; x * b1 - b2 + c[1])
