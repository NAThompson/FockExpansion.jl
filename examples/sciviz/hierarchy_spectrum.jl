# Hyperspherical-harmonic content of each ψ_kp: projections onto the singlet S-state
# harmonics Y_Kl (Λ² Y_Kl = K(K+4) Y_Kl, K = 2n, n = S³ degree), computed with the chain's
# own quadrature on its Chebyshev grid. The algebraic decay with K comes from the
# coalescence (Kato cusp-type) singularities of the coefficients.
# Run: julia -t 4 --project=examples/sciviz examples/sciviz/hierarchy_spectrum.jl
include(joinpath(@__DIR__, "hierarchy_common.jl"))

const KMAX = 12
const KDEG = 80   # largest harmonic degree K (n = 40)

g, ψ = he_chain(; n = 96, kmax = KMAX)
ks = keys_upto(KMAX)
Ks = 0:2:KDEG
# normalised harmonics (unit mean square over S³), singlet only
Hs = Dict{Int,Vector{Matrix{Float64}}}()
for K in Ks
    Hs[K] = [Y ./ rms(g, Y) for (l, Y) in singlet_harmonics(g, K)]
end
# orthonormality check of the quadrature
W = g.W ./ sum(g.W)
allH = reduce(vcat, [Hs[K] for K in Ks])
G = [sum(W .* a .* b) for a in allH, b in allH]
orth = maximum(abs, G - I)
@printf(
    "%d singlet harmonics up to K = %d; max |Gram - I| = %.1e\n",
    length(allH),
    KDEG,
    orth
)

# amplitude per degree: A_K = sqrt(Σ_l c_Kl²) relative to the RMS of ψ_kp
amp = Dict{Tuple{Int,Int},Vector{Float64}}()
for key in ks
    f = ψ[key]
    nf = rms(g, f)
    amp[key] = [sqrt(sum(sum(W .* f .* Y)^2 for Y in Hs[K]; init = 0.0)) / nf for K in Ks]
    captured = sum(abs2, amp[key])
    key[1] <= 4 && @printf(
        "ψ%d%d: fraction of ‖ψ‖² captured up to K = %d: 1 - %.1e\n",
        key...,
        KDEG,
        1 - captured
    )
end

# fitted power law A_K ~ K^-s over 20 ≤ K ≤ 60
function slope(a)
    sel = [i for (i, K) in enumerate(Ks) if 20 <= K <= 60 && a[i] > 1e-15]
    x = log.(Ks[sel])
    y = log.(a[sel])
    A = [ones(length(x)) x]
    -(A\y)[2]
end

fig = Figure(size = (1250, 640), fontsize = 15)
ax1 = Axis(
    fig[1, 1];
    xlabel = L"\mathrm{harmonic\ degree\ } K",
    ylabel = "",
    yticks = (1:length(ks), [latexstring("\\psi_{$k,$p}") for (k, p) in ks]),
    yticklabelsize = 10,
    yreversed = true,
    title = L"\mathrm{(a)}\ \log_{10}\left(A_K/\|\psi_{k,p}\|\right)",
)
M = [log10(max(amp[key][i], 1e-17)) for i in eachindex(Ks), key in ks]
hm = heatmap!(ax1, collect(Ks), 1:length(ks), M; colormap = :viridis, colorrange = (-12, 0))
Colorbar(fig[1, 2], hm; label = L"\log_{10}(A_K/\|\psi\|)")
for (j, (k, p)) in enumerate(ks)
    p == 0 && k > 0 && hlines!(ax1, [j - 0.5]; color = :white, linewidth = 0.6)
end

# Tail norm T_K = sqrt(Σ_{K' > K} A_K'²): the relative L² error of truncating the harmonic
# expansion after degree K. Smooth, unlike A_K itself (A_K alternates with K mod 4).
tail(a) = [sqrt(sum(abs2, a[(i+1):end]; init = 0.0)) for i in eachindex(a)]
function tail_slope(t)
    sel = [i for (i, K) in enumerate(Ks) if 16 <= K <= 56 && t[i] > 1e-14]
    A = [ones(length(sel)) log.(Ks[sel])]
    -(A\log.(t[sel]))[2]
end
ax2 = Axis(
    fig[1, 3];
    xscale = log10,
    yscale = log10,
    xlabel = L"K",
    xticks = [2, 4, 8, 16, 32, 56],
    ylabel = L"\|\psi-\mathcal{P}_{\leq K}\psi\|/\|\psi\|",
    title = L"\mathrm{(b)\ truncation\ error\ of\ the\ harmonic\ series}",
)
for (i, key) in enumerate([(1, 0), (2, 0), (3, 0), (4, 0), (6, 0), (4, 1), (3, 1)])
    t = tail(amp[key])
    sel = [j for j in eachindex(Ks) if Ks[j] > 0 && t[j] > 1e-14 && Ks[j] <= KDEG - 24]
    lines!(
        ax2,
        Ks[sel],
        t[sel];
        color = SERIES[i],
        linewidth = 2,
        linestyle = key[2] > 0 ? :dash : :solid,
        label = latexstring(@sprintf("\\psi_{%d,%d}:\\ K^{-%.2f}", key..., tail_slope(t))),
    )
end
axislegend(ax2; position = :lb, framevisible = false, labelsize = 13)
colsize!(fig.layout, 3, Relative(0.42))
savefig("hierarchy_spectrum", fig)
for key in ks
    @printf(
        "ψ%d%d amplitude slope %.2f, tail slope %.2f\n",
        key...,
        slope(amp[key]),
        tail_slope(tail(amp[key]))
    )
end

path = joinpath(OUT, "hierarchy_spectrum.nc")
isfile(path) && rm(path)
NCDataset(
    path,
    "c";
    attrib = [
        "title" => "Hyperspherical harmonic amplitudes A_K of psi_kp (relative to RMS)",
        "gram_error" => orth,
    ],
) do ds
    defVar(
        ds,
        "K",
        collect(Ks),
        ("K",);
        attrib = ["long_name" => "harmonic degree K (Λ^{2} = K(K+4))"],
    )
    defVar(ds, "k", collect(0:KMAX), ("k",); attrib = ["long_name" => "power of R"])
    defVar(
        ds,
        "p",
        collect(0:(KMAX÷2)),
        ("p",);
        attrib = ["long_name" => "power of log R"],
    )
    A = fill(NaN, length(Ks), KMAX ÷ 2 + 1, KMAX + 1)
    for (k, p) in ks
        A[:, p+1, k+1] = amp[(k, p)]
    end
    defVar(
        ds,
        "amplitude",
        A,
        ("K", "p", "k");
        attrib = field_attrib(
            A,
            [
                "long_name" => "A_{K}/‖ψ_{k,p}‖",
                view_hint(
                    vis = "line",
                    scale = "log",
                    x_dim = "K",
                    multi_dim = "p",
                    envelope = "cummax",
                    envelope_direction = "reverse",
                )...,
            ],
        ),
    )
end
