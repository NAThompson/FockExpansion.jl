# Multipoint QR ODE search for B_A (u variable, arc b), then dump the null vector for the best (r,d) gap.
include(joinpath(@__DIR__, "guess_ode.jl"))
include(joinpath(@__DIR__, "logcoef.jl"))
R1 = deserialize(ENV["F1"]); kp = (4, 0)
v = R1.out[(0, 0)][kp]; c1 = chebcoef(v)
println("tail |c| (last 5): ", [Float64(abs(x)) for x in c1[end-4:end]])
if haskey(ENV, "F2")
    c2 = chebcoef(deserialize(ENV["F2"]).out[(0, 0)][kp])
    println("|F1-F2| = ", Float64(maximum(abs(chebeval(c1, T(x)) - chebeval(c2, T(x))) for x in (-0.95, -0.5, 0.0, 0.5, 0.95))))
end
tlo = T(parse(Float64, get(ENV, "TLO", "-0.95"))); thi = T(parse(Float64, get(ENV, "THI", "-0.05"))); cen = (tlo + thi) / 2; hw = (thi - tlo) / 2
best = nothing
for (r, d) in eval(Meta.parse(ENV["RD"]))
    tot = (r + 1) * (d + 1); np = tot + 60
    rows = Vector{Vector{T}}()
    for i = 1:np
        t0 = tlo + (thi - tlo) * (i - 1) / (np - 1)
        s = series_at(c1, t0, r)
        push!(rows, [((t0 - cen) / hw)^j * s[k+1] * factorial(big(k)) for k = 0:r for j = 0:d])
    end
    A = permutedims(hcat(rows...)); sc = [norm(A[:, k]) for k = 1:size(A, 2)]; A ./= sc'
    F = qr(A, ColumnNorm()); dg = abs.(diag(F.R)) ./ abs(F.R[1, 1])
    ratio = dg[end] / dg[end-1]
    @printf("r=%2d d=%3d unk=%3d  last %.1e  prev %.1e  ratio %.1e\n", r, d, tot, Float64(dg[end]), Float64(dg[end-1]), Float64(ratio))
    flush(stdout)
    if ratio < T(10)^-25
        Rm = F.R; m = size(Rm, 2); x = zeros(T, m); x[m] = 1
        x[1:m-1] = -(UpperTriangular(Rm[1:m-1, 1:m-1]) \ Rm[1:m-1, m])
        w = zeros(T, m); w[F.p] = x; w ./= sc; w ./= w[argmax(abs.(w))]
        open("nullvec_A240_r$(r)_d$(d).txt", "w") do io
            for k = 0:r, j = 0:d
                println(io, k, " ", j, " ", string(w[k*(d+1)+j+1]))
            end
        end
        println("FOUND r=$r d=$d; null vector saved"); break
    end
end
