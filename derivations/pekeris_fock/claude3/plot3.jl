# Phase 3 decisive plot: E error vs basis size (log–log and semilog), all variants at 320–384 bits.
using CairoMakie
function rows(path, v)
    out = Dict{Int,Float64}()
    isfile(path) || return Int[], Float64[]
    for l in eachline(path)
        m = match(Regex("^" * v * "\\s+N=\\s*(\\d+)\\s+size=\\s*(\\d+)\\s+E-E0=\\s*(\\S+)"), l)
        m === nothing || (out[parse(Int, m[2])] = parse(Float64, m[3]))
    end
    k = sort(collect(keys(out)))
    k, [out[i] for i in k]
end
ink, muted, bg = colorant"#0b0b0b", colorant"#52514e", colorant"#fcfcfb"
series = [("../claude2/run_atoms.out", "none", colorant"#1baf7a", "polynomial", :solid),
    ("../claude2/run_atoms.out", "log1", colorant"#eb6834", "+ log σ", :solid),
    ("../claude2/run_atoms.out", "log2", colorant"#eda100", "+ log² σ", :solid),
    ("../claude2/run_atoms.out", "atoms2", colorant"#2a78d6", "+ χ, log(R/s)  (order 2)", :solid),
    ("run3_a2psi.out", "atoms2psi", colorant"#2a78d6", "+ χ, log(R/s), ψ₃₀  (orders 2–3)", :dash),
    ("run3_min5.out", "min5", colorant"#6250d6", "+ minimal order-2 atoms (L̃, α, V, U)", :solid),
    ("run3_min5psi.out", "min5psi", colorant"#6250d6", "+ minimal order-2 atoms + ψ₃₀", :dash)]
fig = Figure(size = (1200, 540), backgroundcolor = bg, fontsize = 15)
for (col, (xs, ttl)) in enumerate(((log10, "log–log"), (identity, "semilog")))
    ax = Axis(fig[1, col]; xscale = xs, yscale = log10, xlabel = "number of basis functions",
        ylabel = col == 1 ? "E − E₀  (hartree)" : "", title = ttl, xticklabelcolor = muted, yticklabelcolor = muted,
        xgridcolor = (:black, 0.06), ygridcolor = (:black, 0.06), backgroundcolor = bg,
        topspinevisible = false, rightspinevisible = false)
    for (f, v, c, lab, ls) in series
        n, e = rows(f, v)
        isempty(n) && continue
        lines!(ax, n, e; color = c, linewidth = 2, label = lab, linestyle = ls)
        scatter!(ax, n, e; color = c, markersize = 9, strokecolor = bg, strokewidth = 2)
    end
    col == 1 && axislegend(ax; position = :lb, framevisible = false, labelcolor = ink, labelsize = 13)
end
Label(fig[0, :], "Helium ground state, fixed ε = 1.704: Ritz error with exact Fock angular atoms (all E ≥ E₀)"; fontsize = 17, color = ink, font = :bold)
save("convergence3.png", fig; px_per_unit = 2)
println("saved")
