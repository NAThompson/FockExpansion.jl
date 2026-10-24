# Convergence figure for a21: Steklov–Ritz matching versus Codex's hyperspherical channels.
using CairoMakie
using Printf

const REF = 0.476747879
const BLUE, ORANGE, AQUA, YELLOW, MAGENTA = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
const INK, MUTED = "#0b0b0b", "#52514e"

# Codex fixed-energy channel sequence (r = 0.75), d = 4, 8, ..., 24.
codex_d = [4, 8, 12, 16, 20, 24]
codex_a21 = [0.476594728101, 0.476737464886, 0.476746090427, 0.476747395497,
             0.476747707937, 0.476747806703]

# Parse "R*=... L=.. <basis> d=..: ... a21=..." lines.
function parse_log(file)
    rows = NamedTuple[]
    for line in eachline(file)
        m = match(r"R\*=([\d.]+) L=\s*(\d+) (\S+)\s+(?:nb=\s*\d+ )?d=\s*(\d+):.*a21=([\d.]+)", line)
        m === nothing && continue
        push!(rows, (Rs = parse(Float64, m[1]), L = parse(Int, m[2]), basis = m[3],
                     d = parse(Int, m[4]), a21 = parse(Float64, m[5])))
    end
    return rows
end

rows = vcat(parse_log("sweep_d_R1.log"), parse_log("sweep_d_R075.log"),
            parse_log("sweep_d_R15.log"), parse_log("sweep_d_R075_L12.log"),
            parse_log("sweep_d_R15_L12.log"), parse_log("sweep_d_R1_B10P6.log"),
            parse_log("sweep_d_R1_B10P6_L12.log"), parse_log("extbasis_R1.log"),
            parse_log("extbasis_R1_hiL.log"), parse_log("extbasis_R1_L14.log"))
pick(f) = sort(unique(r -> (r.d, r.L, r.Rs, r.basis), filter(f, rows)); by = r -> r.d)

fig = Figure(size = (1300, 980), fontsize = 15, backgroundcolor = "#fcfcfb")

# (a) a21 versus interior degree d, against the reference.
ax1 = Axis(fig[1, 1]; xlabel = "interior channel degree d  (channels ≈ (d+2)²/4)",
           ylabel = "a21 − 0.4767478   (×10⁻⁹)",
           title = "(a) a21 vs interior degree, R* = 1",
           titlealign = :left, xgridcolor = (:gray, 0.15), ygridcolor = (:gray, 0.15))
s(x) = (x - 0.4767478) * 1e9
band!(ax1, [3, 50], fill(s(REF - 1e-9), 2), fill(s(REF + 1e-9), 2); color = (:gray, 0.25))
hlines!(ax1, [s(REF)]; color = MUTED, linestyle = :dash, linewidth = 1.5)
text!(ax1, 49, s(REF) - 30; text = "reference 0.476747879 (±1e-9 band)", align = (:right, :top),
      color = MUTED, fontsize = 13)
lines!(ax1, codex_d, s.(codex_a21); color = ORANGE, linewidth = 2)
scatter!(ax1, codex_d, s.(codex_a21); color = ORANGE, markersize = 10, label = "Codex channels (r=0.75)")
for (L, c) in ((6, AQUA), (8, YELLOW), (12, BLUE))
    pts = pick(r -> r.Rs == 1.0 && r.L == L && r.basis == "B12P6")
    lines!(ax1, [p.d for p in pts], s.([p.a21 for p in pts]); color = c, linewidth = 2)
    scatter!(ax1, [p.d for p in pts], s.([p.a21 for p in pts]); color = c, markersize = 10,
             label = "Steklov L = $L")
end
ylims!(ax1, -450, 250)
xlims!(ax1, 3, 50)
axislegend(ax1; position = :rb, framevisible = false)

# (b) successive increments |a21(d) − a21(d−4)| on log–log axes.
ax2 = Axis(fig[1, 2]; xscale = log10, yscale = log10, xlabel = "interior channel degree d",
           ylabel = "|a21(d) − a21(d−4)|", title = "(b) increments in d (algebraic, ~d^-5.9)",
           titlealign = :left, xgridcolor = (:gray, 0.15), ygridcolor = (:gray, 0.15),
           xticks = [8, 12, 16, 24, 32, 48])
lines!(ax2, codex_d[2:end], abs.(diff(codex_a21)); color = ORANGE, linewidth = 2)
scatter!(ax2, codex_d[2:end], abs.(diff(codex_a21)); color = ORANGE, markersize = 10,
         label = "Codex channels")
pts = pick(r -> r.Rs == 1.0 && r.L == 12 && r.basis == "B12P6")
dd = [p.d for p in pts]
aa = [p.a21 for p in pts]
lines!(ax2, dd[2:end], abs.(diff(aa)); color = BLUE, linewidth = 2)
scatter!(ax2, dd[2:end], abs.(diff(aa)); color = BLUE, markersize = 10, label = "Steklov L = 12")
xg = [16.0, 48.0]
lines!(ax2, xg, 2.5e-7 .* (xg ./ 18) .^ -5.87; color = MUTED, linestyle = :dot, linewidth = 1.5,
       label = "∝ d^-5.87")
axislegend(ax2; position = :lb, framevisible = false)

# (c) trace degree L at fixed d = 24 and R* = 0.75, 1, 1.5.
ax3 = Axis(fig[2, 1]; xlabel = "trace degree L (harmonics up to K = 2L)",
           ylabel = "a21 − 0.4767478   (×10⁻⁹)", title = "(c) a21 vs trace degree, d = 24 and 40",
           titlealign = :left, xgridcolor = (:gray, 0.15), ygridcolor = (:gray, 0.15))
hlines!(ax3, [s(REF)]; color = MUTED, linestyle = :dash, linewidth = 1.5)
for (Rs, bas, d, c, lab) in ((1.0, "B12P6", 24, BLUE, "R* = 1, d = 24"),
                             (1.0, "B10P6", 40, AQUA, "R* = 1, d = 40"),
                             (0.75, "B10P6", 40, MAGENTA, "R* = 0.75, d = 40"),
                             (1.5, "B10P6", 40, YELLOW, "R* = 1.5, d = 40"))
    pts = sort(unique(r -> r.L, filter(r -> r.Rs == Rs && r.basis == bas && r.d == d && r.L >= 6, rows));
               by = r -> r.L)
    isempty(pts) && continue
    lines!(ax3, [p.L for p in pts], s.([p.a21 for p in pts]); color = c, linewidth = 2)
    scatter!(ax3, [p.L for p in pts], s.([p.a21 for p in pts]); color = c, markersize = 10, label = lab)
end
axislegend(ax3; position = :rt, framevisible = false)

# (d) exterior basis: spread at fixed (R* = 1, d = 24) for each L.
ax4 = Axis(fig[2, 2]; xlabel = "trace degree L", ylabel = "a21(basis) − a21(B12P6)   (×10⁻⁹)",
           title = "(d) exterior basis dependence, R* = 1, d = 24 (B8P2 off scale, +0.3 to +1.6e-6)", titlealign = :left,
           xgridcolor = (:gray, 0.15), ygridcolor = (:gray, 0.15))
hlines!(ax4, [0.0]; color = MUTED, linewidth = 1)
for (bas, c) in (("B8P2", ORANGE), ("B10P4", AQUA), ("B10P6", YELLOW), ("B12P4", MAGENTA))
    pts = []
    for L in (4, 6, 8, 10, 12)
        a = filter(r -> r.Rs == 1.0 && r.L == L && r.basis == bas && r.d == 24, rows)
        b = filter(r -> r.Rs == 1.0 && r.L == L && r.basis == "B12P6" && r.d == 24, rows)
        (isempty(a) || isempty(b)) && continue
        push!(pts, (L, 1e9 * (a[1].a21 - b[1].a21)))
    end
    isempty(pts) && continue
    scatterlines!(ax4, first.(pts), last.(pts); color = c, linewidth = 2, markersize = 10, label = bas)
end
ylims!(ax4, -40, 40)
axislegend(ax4; position = :rb, framevisible = false, nbanks = 2)

save("a21_convergence.png", fig; px_per_unit = 2)
println("saved a21_convergence.png")
