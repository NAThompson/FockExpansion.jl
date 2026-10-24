# Run with: julia --project=/Users/nathompson7/FockExpansion.jl/examples plucker_plot.jl
using CairoMakie, DelimitedFiles
const DIR = @__DIR__
raw = readdlm(joinpath(DIR, "plucker_conv.csv"), ','; skipstart = 1)
sel(chart, om, R) = (r = [i for i = 1:size(raw, 1) if raw[i, 1] == chart && raw[i, 2] == om && raw[i, 3] == R];
    (Float64.(raw[r, 4]), max.(Float64.(raw[r, 5]), 1e-17)))
fig = Figure(size = (1300, 950), fontsize = 14)
colors = Makie.wong_colors()
for (col, om) in enumerate((0.224, 0.4))
    ax = Axis(fig[1, col]; xscale = log10, yscale = log10, xlabel = "terms K", ylabel = "error",
        title = "d=2, ω=$om: Plücker angle (solid) vs Y-chart ‖ΔY‖ (dashed)")
    for (i, R) in enumerate((0.5, 1.0, 1.5, 2.0, 3.0))
        K, e = sel("plucker_angle", om, R); lines!(ax, K[2:end], e[2:end]; color = colors[i], label = "R=$R")
        K, e = sel("Ychart", om, R); lines!(ax, K[2:end], e[2:end]; color = colors[i], linestyle = :dash)
    end
    ylims!(ax, 1e-16, 1e5)
    axislegend(ax; position = :lb)
end
# fitted abscissa p = R + omega * slope(K in 400..800)
ax = Axis(fig[2, 1]; xlabel = "R", ylabel = "fitted p = R + ω·dlog(err)/dlog(K)",
    title = "Observed abscissa (K = 400–800)")
for (j, om) in enumerate((0.224, 0.3, 0.4)), (chart, ls) in (("plucker_angle", :solid), ("Ychart", :dash))
    Rs = [1.0, 1.5, 2.0, 3.0]; ps = Float64[]
    for R in Rs
        K, e = sel(chart, om, R); k = 401:801
        x = log.(K[k]); y = log.(e[k]); s = sum((x .- sum(x) / length(x)) .* (y .- sum(y) / length(y))) / sum((x .- sum(x) / length(x)) .^ 2)
        push!(ps, R + om * s)
    end
    scatterlines!(ax, Rs, ps; color = colors[j], linestyle = ls, label = "$(chart == "Ychart" ? "Y" : "Plücker"), ω=$om")
end
hlines!(ax, [0.6373]; color = :gray, linestyle = :dot)
text!(ax, 2.2, 0.66; text = "Re(Riccati pole) = 0.637", color = :gray)
axislegend(ax; position = :lt)
P = readdlm(joinpath(DIR, "plucker_pade.csv"), ','; skipstart = 1)
ax = Axis(fig[2, 2]; xlabel = "Re t", ylabel = "Im t", title = "Borel plane: Padé [60/60] poles", aspect = DataAspect())
for (name, mk, c) in (("plucker", :circle, colors[1]), ("Ychart", :xcross, colors[2]))
    r = [i for i = 1:size(P, 1) if P[i, 1] == name && P[i, 2] == 60]
    scatter!(ax, Float64.(P[r, 3]), Float64.(P[r, 4]); marker = mk, color = c, label = name)
end
vlines!(ax, [-3.556]; color = :gray, linestyle = :dot)
xlims!(ax, -8, 3); ylims!(ax, -4, 4)
axislegend(ax; position = :rt)
save(joinpath(DIR, "plucker.png"), fig)
println("saved")
