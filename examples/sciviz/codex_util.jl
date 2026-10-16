# Small helpers shared by the codex_*.jl scripts (CSV in/out, symmetric log colour).
const latexstring = Makie.LaTeXStrings.latexstring

"""sign(x) log₁₀(1 + |x|): a symmetric log scale for colour maps of signed data."""
symlog(x) = sign(x) * log10(1 + abs(x))
symlog_ticks(n) = (collect(-n:n), [k == 0 ? L"0" : latexstring(k < 0 ? "-10^{$(-k)}" : "10^{$k}") for k = -n:n])

function write_rows_csv(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ","))
        for r in rows
            println(io, join(r, ","))
        end
    end
    println("wrote ", path)
end

"""Write same-size matrices as long-form CSV: i, j, then one column per name."""
function write_matrix_csv(path, mats::Dict)
    names = sort(collect(keys(mats)))
    A = first(values(mats))
    rows = [(i, j, (mats[k][i, j] for k in names)...) for j in axes(A, 2) for i in axes(A, 1)]
    write_rows_csv(path, ["i"; "j"; names], rows)
end

"""Read a CSV with a header into a Dict of column vectors (Float64 where parseable)."""
function read_csv(path)
    lines = readlines(path)
    header = split(lines[1], ",")
    cols = [split(l, ",") for l in lines[2:end]]
    # BigFloat fallback: mpmath writes values like 1.0e-324 that underflow Float64.
    num(s) = something(tryparse(Float64, s), (b = tryparse(BigFloat, s); b === nothing ? String(s) : Float64(b)))
    Dict(String(h) => [num(c[k]) for c in cols] for (k, h) in enumerate(header))
end

function read_matrix_csv(path)
    d = read_csv(path)
    n, m = Int(maximum(d["i"])), Int(maximum(d["j"]))
    Dict(k => reshape(Float64.(v), n, m) for (k, v) in d if k ∉ ("i", "j"))
end
