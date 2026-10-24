# Phase 2: E error vs basis size for polynomial / +log / +log² / +order-2 atoms / +order-3,4 atoms; log-log and semilog.
using CairoMakie
function rows(path, v; Nmax = 99)
    out = Dict{Int,Float64}()
    for l in eachline(path)
        m = match(Regex("^" * v * "\\s+N=\\s*(\\d+)\\s+size=\\s*(\\d+)\\s+E-E0=\\s*(\\S+)"), l)
        m === nothing && continue
        parse(Int, m[1]) <= Nmax && (out[parse(Int, m[2])] = parse(Float64, m[3]))
    end
    k = sort(collect(keys(out)))
    k, [out[i] for i in k]
end
ink, muted, bg = colorant"#0b0b0b", colorant"#52514e", colorant"#fcfcfb"
series = [("none", colorant"#1baf7a", "polynomial"), ("log1", colorant"#eb6834", "+ log σ"),
    ("log2", colorant"#eda100", "+ log² σ"), ("atoms2", colorant"#2a78d6", "+ order-2 atoms (χ, log R/s)"),
    ("atoms4", colorant"#6250d6", "+ order-3/4 atoms (ψ₃₀ Float64; N ≤ 7 only)")]
fig = Figure(size = (1200, 520), backgroundcolor = bg, fontsize = 15)
for (col, (xs, ttl)) in enumerate(((log10, "log–log"), (identity, "semilog")))
    ax = Axis(fig[1, col]; xscale = xs, yscale = log10, xlabel = "number of basis functions",
        ylabel = col == 1 ? "E − E₀  (hartree)" : "", title = ttl, xticklabelcolor = muted, yticklabelcolor = muted,
        xgridcolor = (:black, 0.06), ygridcolor = (:black, 0.06), backgroundcolor = bg,
        topspinevisible = false, rightspinevisible = false)
    for (v, c, lab) in series
        n, e = v == "atoms4" ? rows("run_atoms4.out", v; Nmax = 7) : rows("run_atoms.out", v)
        isempty(n) && continue
        keep = e .> 0
        lines!(ax, n[keep], e[keep]; color = c, linewidth = 2, label = lab)
        scatter!(ax, n[keep], e[keep]; color = c, markersize = 9, strokecolor = bg, strokewidth = 2)
    end
    col == 1 && axislegend(ax; position = :lb, framevisible = false, labelcolor = ink)
end
Label(fig[0, :], "Helium ground state, fixed ε = 1.704: Ritz error with exact Fock angular atoms"; fontsize = 17, color = ink, font = :bold)
save("convergence_atoms.png", fig; px_per_unit = 2)
println("saved")
