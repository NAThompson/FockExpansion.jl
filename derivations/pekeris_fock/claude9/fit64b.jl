include("fit64.jl")
n = 6
for sel in (1:1, 1:7)
    rows = map(geo) do (r1, r2, r12)
        A = atoms(r1, r2, r12)[sel]
        m = monos(r1, r2, r12, n)
        vcat([a .* m for a in A]...)
    end
    M = permutedims(hcat(rows...))
    c = qr(M, ColumnNorm()) \ vals
    r = M * c - vals
    println("sel=$sel rel ", norm(r) / norm(vals))
    p = sortperm(abs.(r), rev = true)
    for i in p[1:8]
        println("  ", pts[i], "  r=", r[i], "  val=", vals[i])
    end
    println("  median |r| = ", sort(abs.(r))[end÷2])
end
