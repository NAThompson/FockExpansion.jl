# Run with: julia --project=/Users/nathompson7/FockExpansion.jl/examples plot_figure.jl
using CairoMakie, DelimitedFiles
const DIR = @__DIR__
readcsv(name) = readdlm(joinpath(DIR, name), ',', Float64; skipstart = 1)
floor_err(e) = max(e, 1e-17)

fig = Figure(size = (1200, 900), fontsize = 14)
for (col, d) in enumerate((2, 6))
    ax = Axis(fig[1, col]; yscale = log10, xlabel = "terms K", ylabel = "‖Y_K(R) − Y_panels(R)‖",
        title = "Inverse factorial series of Y, d = $d")
    A = readcsv("factorial_d$d.csv")
    for R in unique(A[:, 1])
        rows = A[A[:, 1].==R, :]
        lines!(ax, rows[:, 2], floor_err.(rows[:, 3]); label = "R = $R")
    end
    ylims!(ax, 1e-17, 1e6)
    axislegend(ax; position = :rt)
end
ax = Axis(fig[2, 1]; yscale = log10, xlabel = "Taylor order", ylabel = "‖Y(R0) − Y_panels(R0)‖",
    title = "Factorial series at R1, then one Taylor expansion about R1")
for (d, R1) in ((2, 4), (4, 6), (6, 9))
    A = readcsv("taylor_d$d.csv")
    for R0 in (1.0, 0.75)
        rows = A[A[:, 1].==R0, :]
        lines!(ax, rows[:, 2], floor_err.(rows[:, 3]); label = "d=$d, R1=$R1, R0=$R0",
            linestyle = R0 == 1.0 ? :solid : :dash)
    end
end
ylims!(ax, 1e-17, 1e2)
axislegend(ax; position = :rt)
ax = Axis(fig[2, 2]; yscale = log10, xlabel = "R", ylabel = "relative change / 1e-12",
    title = "Fock-sum sensitivity to 1e-12 perturbations")
for d in (2, 4)
    A = readcsv("sensitivity_d$d.csv")
    scatterlines!(ax, A[:, 1], A[:, 2]; label = "d=$d, δa")
    scatterlines!(ax, A[:, 1], A[:, 3]; label = "d=$d, δE", linestyle = :dash)
end
axislegend(ax; position = :lt)
save(joinpath(DIR, "exterior_factorial.png"), fig)
println("saved")
