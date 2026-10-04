"""Cached corrected Langner double-series coefficients for a specified Z.
The distributed table supports N≤60, J≤320; loading is excluded from warmed timings.
"""
struct LangnerTable{T}
    Z::T
    unrotated::Matrix{T}
    rotated::Matrix{T}
    bminus::Matrix{T}
end
function LangnerTable(; Z = 2.0, T = Float64)
    arrays=ntuple(_->zeros(T, 61, 321), 3)
    z=T(Z)
    path=joinpath(@__DIR__, "..", "data", "langner_coefficients.tsv")
    for line in eachline(path)
        startswith(line, "#") && continue
        kind, n, k, l, q=split(line)
        j=parse(Int, kind)+1
        arrays[j][parse(Int, n)+1, parse(Int, k)+1]=j==3 ? parse(T, l) :
                                                    z*parse(T, l)+z^2*parse(T, q)
    end
    LangnerTable(z, arrays...)
end
function double_horner(c, x, y, N, J)
    outer=zero(x+y+c[1, 1])
    for n = (N+1):-1:1
        inner=zero(outer)
        for k = (J+1):-1:1
            inner=muladd(inner, y, c[n, k])
        end
        outer=muladd(outer, x, inner)
    end
    outer
end
function double_direct(c, x, y, N, J)
    xp=[x^n for n = 0:N]
    yp=[y^k for k = 0:J]
    sum(c[n+1, k+1]*xp[n+1]*yp[k+1] for n = 0:N for k = 0:J)
end
# Cached-power dot products expose independent terms to SIMD.
function double_powers(c, x, y, N, J)
    xp=Vector{typeof(x+y+c[1, 1])}(undef, N+1)
    yp=similar(xp, J+1)
    xp[1]=one(x)
    yp[1]=one(y)
    for n = 2:(N+1)
        xp[n]=xp[n-1]*x
    end
    for k = 2:(J+1)
        yp[k]=yp[k-1]*y
    end
    total=zero(xp[1])
    for k = 1:(J+1)
        row=zero(total)
        @inbounds @simd for n = 1:(N+1)
            row+=c[n, k]*xp[n]
        end
        total+=row*yp[k]
    end
    total
end
function langner_head(a, d, Z, E, a21, bm)
    π=typedpi(a)
    r1=sqrt((1+a)/2)
    r2=sqrt((1-a)/2)
    s=r1+r2
    z=sqrt(2-d*d)
    w=1-d*d
    c=a21+Z*(-oftype(a, 17)/72+(24clausen2(π/2)-31)/(36π))
    aw=iszero(w) ? one(w) : asin(w)/w
    anglepart=aw/12+w*aw^2/(12π)-d*z*acos(d/sqrt(oftype(a, 2)))/(3π)
    ans=E*(Z*s*(2+r1*r2)/18-(6-d*d)*d/72)-Z^3*s*(1+5r1*r2)/18
    ans-=c*(Z*s*w/2-(6-5d*d)*d/12)
    ans+=(3-2d*d)*d/72
    # Corrected cubic sign and block-to-total-projection convention.
    ans+=Z*(
        -(6-5d*d)*d*log(s+d)/36-a*d*bm/24-5d^3/(108π)-d^3/54-5s*d*d/36+(1+10r1*r2)*d/72+s^3/108+d*d*z*acos(
            d/sqrt(oftype(a, 2)),
        )/(6π)
    )
    ans+=Z^2*(
        s*w*(log(s+d)/6+1/(12π))+(r1-r2)*(1+oftype(a, 2.5)*r1*r2)*bm/9+s*anglepart-29d^3/216+s*d*d/24+(
            17+20r1*r2
        )*d/36-s*(47-2r1*r2)/216
    )
    ans
end
"""Corrected Langner rectangular truncation. No automatic convergence guarantee.
`method=:powers` uses vectorized cached powers; `:horner` avoids power arrays;
`:direct` retains the baseline summation.
"""
function psi30_langner(α, θ, table::LangnerTable; E, a21, N = 60, J = 320, method = :powers)
    0<=N<=60 && 0<=J<=320 || throw(ArgumentError("Table supports 0≤N≤60, 0≤J≤320"))
    method in (:horner, :direct, :powers) || throw(ArgumentError("Unknown series method"))
    g=geometry(α, θ)
    a=cos(g.α)
    d=g.ξ
    y=-d/sqrt(oftype(a, 2))
    x=(a/2)^2
    ev=method==:horner ? double_horner : method==:powers ? double_powers : double_direct
    bm=(a/2)*ev(table.bminus, x, y, N, J)
    langner_head(a, d, table.Z, E, a21, bm)+ev(table.unrotated, x, y, N, J)+ev(
        table.rotated,
        g.v/2,
        -g.r1,
        N,
        J,
    )+ev(table.rotated, g.v/2, -g.r2, N, J)
end
