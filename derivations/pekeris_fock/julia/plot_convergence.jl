# Stage 4(d): energy error vs shell N, with and without the log(sigma) sector.
using CairoMakie

function read_pekeris33(path)
    N = Int[]; err = Float64[]
    for line in eachline(path)
        f = split(line)
        (isempty(f) || tryparse(Int, f[1]) === nothing) && continue
        push!(N, parse(Int, f[1])); push!(err, parse(Float64, f[end]))
    end
    N, err
end

function read_log2(path)
    N = Int[]; err = Float64[]
    for line in eachline(path)
        mt = match(r"^fixedlog2\s+N=\s*(\d+).*E-E0=(\S+)", line)
        mt === nothing && continue
        push!(N, parse(Int, mt[1])); push!(err, parse(Float64, mt[2]))
    end
    N, err
end

function read_galerkin(path, variant)
    N = Int[]; err = Float64[]
    for line in eachline(path)
        f = split(line)
        (length(f) < 5 || f[1] != variant) && continue
        push!(N, parse(Int, f[2])); push!(err, parse(Float64, f[5]))
    end
    N, err
end

blue, orange, aqua, yellow = colorant"#2a78d6", colorant"#eb6834", colorant"#1baf7a", colorant"#eda100"
ink, muted = colorant"#0b0b0b", colorant"#52514e"

fig = Figure(size = (820, 540), backgroundcolor = colorant"#fcfcfb", fontsize = 15)
ax = Axis(fig[1, 1]; yscale = log10, xlabel = "shell N (polynomial degree l+m+n ≤ N)",
    ylabel = "E(N) − E₀  (hartree)", title = "Helium ground state: Pekeris space with and without log(u+v+w) sectors",
    titlecolor = ink, xlabelcolor = ink, ylabelcolor = ink, xticklabelcolor = muted, yticklabelcolor = muted,
    xgridcolor = (:black, 0.06), ygridcolor = (:black, 0.06), xticks = 2:2:16,
    backgroundcolor = colorant"#fcfcfb", topspinevisible = false, rightspinevisible = false)
series = [
    (read_galerkin("galerkin_log.out", "none"), aqua, "Rayleigh–Ritz, ε = 1.704 (no log)"),
    (read_pekeris33("pekeris33.out"), blue, "33-term recurrence, ε² = −E (no log; overlaps green)"),
    (read_galerkin("galerkin_log.out", "fixedlog"), orange, "Rayleigh–Ritz, ε = 1.704, + log σ sector (Fock R² log R fixed)"),
    (read_log2("log2_test.out"), yellow, "  … + log² σ sector (polynomials of degree 4..N−2)"),
]
for ((N, err), c, lab) in series
    lines!(ax, N, err; color = c, linewidth = 2, label = lab, linestyle = c == blue ? :dash : :solid)
    scatter!(ax, N, err; color = c, markersize = 9, strokecolor = colorant"#fcfcfb", strokewidth = 2)
end
N, err = series[4][1]
text!(ax, N[end], err[end]; text = "  4.7e-14", color = ink, align = (:left, :center), fontsize = 13)
N, err = series[3][1]
text!(ax, N[end], err[end]; text = "  9.3e-13", color = ink, align = (:left, :center), fontsize = 13)
N, err = series[2][1]
text!(ax, N[end], err[end]; text = "  1.4e-8", color = ink, align = (:left, :center), fontsize = 13)
xlims!(ax, 1.5, 17.5)
axislegend(ax; position = :lb, framevisible = false, labelcolor = ink)
save("convergence.png", fig; px_per_unit = 2)
println("saved")
