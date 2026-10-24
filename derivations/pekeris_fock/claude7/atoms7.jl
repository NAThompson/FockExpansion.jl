# Phase 7 numerics: split psi30 into its structured pieces (classical + four elliptic panels), each a separate atom.
include("../claude3/atoms3.jl")
const ZZ = 2
function pieces_xy(x, y)
    # geometry at s = 1 in the chamber r1 >= r2 (x >= 0), alpha <= pi/2 (no fold needed)
    α, μ = angles(x, y); θ = theta(x, y)
    g = FockExpansion.geometry(α, θ)
    (; r1, r2, v, ξ) = g
    rr = sqrt(oftype(x, 2)); yv = -ξ / rr
    c = cos(α); d = v - 2r2^2
    rt = BigFloat(10)^(-PSI30_DIGITS[])
    out = Vector{typeof(x)}(undef, 5)
    k = 1
    for (b, opb, omb, bmy, bpy) in ((c, 2r1^2, 2r2^2, d, c + ξ^2), (-c, 2r2^2, 2r1^2, -c - ξ^2, -d))
        c0 = -rr * ZZ / 576 * (9ZZ * b^2 - 18ZZ - 8b^2 + 8)
        c2 = -rr * ZZ / 288 * omb * (90ZZ * yv^2 - 45ZZ - 18b - 32yv^2 + 16)
        out[k] = FockExpansion.panel(b, yv, c0, c2; rtol = rt, opb, omb, bmy, bpy)[1]; k += 1
    end
    b = -v; opb, omb = ξ^2, g.opv
    for yr in (-r1, -r2)
        out[k] = FockExpansion.panel(b, yr, ZZ^2 * opb * omb / 18, -2ZZ^2 * omb * (2yr^2 - 1) / 9; rtol = rt, opb, omb)[1]; k += 1
    end
    out[5] = FockExpansion.classical(α, 2atan(ξ, sqrt(g.opv)), ZZ)
    hx3 = sqrt((1 + x^2) / 2)^3
    hx3 .* out
end
piece_xy(i) = (x, y) -> pieces_xy(x, y)[i]
const PIECES = (:pan1, :pan2, :pan3, :pan4, :clas)
function build7(N, variant; ε, Na = N - 2)
    B = build3(N, :min5; ε)
    if variant == :min5split
        for g in PIECES, n = 3:Na, a = 0:2:n-3, b = 0:n-3-a
            push!(B, [Term(n, 0, mono(a, b, g))])
        end
    end
    B
end
