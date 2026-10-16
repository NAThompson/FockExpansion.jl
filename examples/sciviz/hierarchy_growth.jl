# Growth (decay) of the Fock coefficients with k, from the spectral chain on its native
# grid: max |ψ_kp| and the S³ root-mean-square, and two empirical radius-of-convergence
# estimates (ratio and root tests) compared with the rigorous bound of the boundary-matching
# notes (|p_n| ≤ 2^153 (n+1)^8 4096^n, i.e. radius ≥ 1/4096 per order).
# Run: julia -t 4 --project=examples/sciviz examples/sciviz/hierarchy_growth.jl
include(joinpath(@__DIR__, "hierarchy_common.jl"))

const KMAX = 32

g, ψ = he_chain(; n = 80, kmax = KMAX)
g2, ψ2 = he_chain(; n = 64, kmax = KMAX)
ks = keys_upto(KMAX)
mx = Dict(key => maximum(abs, ψ[key]) for key in ks)
rm_ = Dict(key => rms(g, ψ[key]) for key in ks)
# grid-independence: relative change of the RMS between n = 64 and n = 80
δ = Dict(key => abs(rms(g2, ψ2[key]) - rm_[key]) / rm_[key] for key in ks)
for k = 0:KMAX
    @printf(
        "k=%2d  max|ψ_k0|=%.3e  rms=%.3e  rel. n-change (worst p)=%.1e\n",
        k,
        mx[(k, 0)],
        rm_[(k, 0)],
        maximum(δ[(k, p)] for p = 0:(k÷2))
    )
end

# Coefficient size per order: M_k = Σ_p max|ψ_kp| (bounds |Σ_p log^p R ψ_kp| at R = 1/e … e).
Mk = [sum(mx[(k, p)] for p = 0:(k÷2)) for k = 0:KMAX]
Rk = [sum(rm_[(k, p)] for p = 0:(k÷2)) for k = 0:KMAX]

fig = Figure(size = (1200, 520), fontsize = 16)
ax1 = Axis(
    fig[1, 1];
    yscale = log10,
    xlabel = L"k",
    ylabel = L"\max_{\alpha,\theta}|\psi_{k,p}|",
    title = L"\mathrm{(a)\ size\ of\ each\ coefficient}",
)
for p = 0:6
    kk = [k for k = 2p:KMAX]
    lines!(
        ax1,
        kk,
        [mx[(k, p)] for k in kk];
        color = SERIES[p+1],
        linewidth = 2,
        label = latexstring("p=$p"),
    )
    scatter!(ax1, kk, [mx[(k, p)] for k in kk]; color = SERIES[p+1], markersize = 7)
end
kk = 0:KMAX
lines!(
    ax1,
    kk,
    [rm_[(k, 0)] for k in kk];
    color = SERIES[1],
    linestyle = :dash,
    linewidth = 1.5,
    label = L"p=0,\ \mathrm{RMS\ over\ } S^3",
)
axislegend(ax1; position = :lb, framevisible = false, labelsize = 13, nbanks = 2)

ax2 = Axis(
    fig[1, 2];
    yscale = log10,
    xlabel = L"k",
    ylabel = L"\mathrm{radius\ estimate}",
    title = L"\mathrm{(b)\ empirical\ radius\ of\ convergence}",
)
k2 = 2:KMAX
lines!(
    ax2,
    k2,
    [Mk[k] / Mk[k+1] for k in k2];
    color = SERIES[1],
    linewidth = 2,
    label = L"\mathrm{ratio}\ M_{k-1}/M_k",
)
scatter!(ax2, k2, [Mk[k] / Mk[k+1] for k in k2]; color = SERIES[1], markersize = 7)
lines!(
    ax2,
    k2,
    [Mk[k+1]^(-1 / k) for k in k2];
    color = SERIES[2],
    linewidth = 2,
    label = L"\mathrm{root}\ M_k^{-1/k}",
)
scatter!(ax2, k2, [Mk[k+1]^(-1 / k) for k in k2]; color = SERIES[2], markersize = 7)
lines!(
    ax2,
    k2,
    [Rk[k] / Rk[k+1] for k in k2];
    color = SERIES[3],
    linewidth = 2,
    linestyle = :dash,
    label = L"\mathrm{ratio\ of\ RMS\ sums}",
)
# Fit ratio ≈ c k^s over k ≥ 10 (s > 0: the radius estimate keeps growing, as it must for
# a series converging at every R, Morgan 1986).
kf = 10:KMAX
ratio = [Mk[k] / Mk[k+1] for k in kf]
A = [ones(length(kf)) log.(kf)]
lc, s = A \ log.(ratio)
@printf("ratio fit over k = 10..%d: %.3f k^%.3f\n", KMAX, exp(lc), s)
hlines!(ax2, [1 / 4096]; color = :black, linestyle = :dot, linewidth = 1.5)
text!(
    ax2,
    2,
    1 / 4096;
    text = L"\mathrm{rigorous\ bound\ } 1/4096",
    align = (:left, :bottom),
    fontsize = 13,
)
lines!(
    ax2,
    k2,
    exp(lc) .* k2 .^ s;
    color = :gray30,
    linestyle = :dashdot,
    linewidth = 1.5,
    label = latexstring(
        @sprintf("\\mathrm{fit}\\ %.2f\\,k^{%.2f}\\ (k\\geq 10)", exp(lc), s)
    ),
)
ylims!(ax2, 1e-4, 30)
axislegend(ax2; position = :rc, framevisible = false, labelsize = 13)
Label(
    fig[2, 1:2],
    L"M_k=\sum_p \max|\psi_{k,p}|.\ \ Z=2,\ E_{\mathrm{He}},\ a_{21}=0.4767;\ \ \mathrm{free\ harmonic\ parts\ of\ }\psi_{k,0}\ (k\geq4)\ \mathrm{set\ to\ zero}.";
    fontsize = 13,
)
savefig("hierarchy_growth", fig)

isfile(joinpath(OUT, "hierarchy_growth.nc")) && rm(joinpath(OUT, "hierarchy_growth.nc"))
NCDataset(
    joinpath(OUT, "hierarchy_growth.nc"),
    "c";
    attrib = ["title" => "Size of Fock coefficients vs k (spectral chain, n = 80)"],
) do ds
    defVar(ds, "k", collect(0:KMAX), ("k",); attrib = ["long_name" => "power of R"])
    defVar(
        ds,
        "p",
        collect(0:(KMAX÷2)),
        ("p",);
        attrib = ["long_name" => "power of log R"],
    )
    # Julia (k, p) is HDF5 (p, k): k is the line view's x axis, p a slider.
    M = fill(NaN, KMAX + 1, KMAX ÷ 2 + 1)
    R = fill(NaN, size(M))
    for (k, p) in ks
        M[k+1, p+1] = mx[(k, p)]
        R[k+1, p+1] = rm_[(k, p)]
    end
    hint = view_hint(vis = "line", scale = "log")
    # All p on one plot, k on the x axis.
    multi = view_hint(x_dim = "k", multi_dim = "p")
    defVar(
        ds,
        "max_abs_psi",
        M,
        ("k", "p");
        attrib = ["long_name" => "max |ψ_{k,p}| (NaN where p > k/2)", hint..., multi...],
    )
    defVar(
        ds,
        "rms_psi",
        R,
        ("k", "p");
        attrib = ["long_name" => "RMS of ψ_{k,p} over S^{3}", hint..., multi...],
    )
    defVar(
        ds,
        "M",
        Mk,
        ("k",);
        attrib = ["long_name" => "Σ_{p} max|ψ_{k,p}|", hint...],
    )
end
