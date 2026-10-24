using CairoMakie, Serialization
data = deserialize("goal1_plot.jls"); Ns = (10, 14, 18, 22)
bg = colorant"#fcfcfb"; ink = colorant"#0b0b0b"; muted = colorant"#52514e"
cs = [colorant"#1baf7a", colorant"#2a78d6", colorant"#6250d6", colorant"#eb6834"]
fig = Figure(size = (1250, 480), backgroundcolor = bg, fontsize = 14)
axs = [Axis(fig[1, i]; yscale = log10, xlabel = "shell k = l+m+n", backgroundcolor = bg, xgridcolor = (:black, 0.06), ygridcolor = (:black, 0.06),
    topspinevisible = false, rightspinevisible = false, xticklabelcolor = muted, yticklabelcolor = muted, title = t)
       for (i, t) in enumerate(("Miller solution at E₀ (zero tail)", "|response| to unit random tail data (shells N+1..N+3)", "Miller solution at E₀ + 0.5i"))]
axs[1].ylabel = "shell norm"
for (j, N) in enumerate(Ns)
    sn, resp, snc, g = data[N]
    k = 0:N
    lines!(axs[1], k, sn; color = cs[j], linewidth = 2, label = "N_max = $N")
    for r in resp
        lines!(axs[2], 1:N, r[2:end]; color = cs[j], linewidth = 1.5)
    end
    lines!(axs[3], k, snc; color = cs[j], linewidth = 2)
end
axislegend(axs[1]; position = :lb, framevisible = false, labelcolor = ink)
Label(fig[0, :], "Plain Pekeris recurrence (57-term, ε = 1.704): no exponential dichotomy — tail data propagate undamped to low shells"; fontsize = 16, font = :bold, color = ink)
save("goal1_separation.png", fig; px_per_unit = 2)
println("saved")
