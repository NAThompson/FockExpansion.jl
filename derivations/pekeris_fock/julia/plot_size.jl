# Same data as convergence.png, against basis size (review point 5: fairer comparison).
using CairoMakie
rows(path, rx) = (m = [match(rx, l) for l in eachline(path)]; m = filter(!isnothing, m);
    (parse.(Int, getindex.(m, 1)), parse.(Float64, getindex.(m, 2))))
p33 = rows("pekeris33.out", r"^\s*\d+\s+(\d+)\s+.*\s(\S+)$")
none = rows("galerkin_log.out", r"^none\s+\d+\s+(\d+)\s+\S+\s+(\S+)")
lg1 = rows("galerkin_log.out", r"^fixedlog\s+\d+\s+(\d+)\s+\S+\s+(\S+)")
lg2 = rows("log2_test.out", r"^fixedlog2\s+N=\s*\d+\s+size=\s*(\d+)\s+E-E0=(\S+)")
ink, muted = colorant"#0b0b0b", colorant"#52514e"
fig = Figure(size = (820, 520), backgroundcolor = colorant"#fcfcfb", fontsize = 15)
ax = Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "number of basis functions (singlet)",
    ylabel = "E − E₀  (hartree)", title = "Same runs against basis size",
    xticklabelcolor = muted, yticklabelcolor = muted, xgridcolor = (:black, 0.06), ygridcolor = (:black, 0.06),
    backgroundcolor = colorant"#fcfcfb", topspinevisible = false, rightspinevisible = false)
for ((n, e), c, lab, ls) in [(none, colorant"#1baf7a", "polynomial, ε = 1.704", :solid),
    (p33, colorant"#2a78d6", "33-term recurrence, ε² = −E", :dash),
    (lg1, colorant"#eb6834", "+ log σ sector", :solid), (lg2, colorant"#eda100", "+ log σ and log² σ sectors", :solid)]
    lines!(ax, n, e; color = c, linewidth = 2, label = lab, linestyle = ls)
    scatter!(ax, n, e; color = c, markersize = 9, strokecolor = colorant"#fcfcfb", strokewidth = 2)
end
axislegend(ax; position = :lb, framevisible = false, labelcolor = ink)
save("convergence_size.png", fig; px_per_unit = 2)
println("saved")
