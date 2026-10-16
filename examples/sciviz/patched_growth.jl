# Proven coefficient growth bounds of the patched construction versus the actual
# sizes of the Fock coefficients computed by FockExpansion.jl (k ≤ 4).
#
#     julia -t 4 --project=examples/sciviz examples/sciviz/patched_growth.jl
#
# The coefficients are evaluated here on this figure's own grids (41², ψ40: 15²) and only
# their maxima are written to output/patched.nc; the fields themselves are in
# output/hierarchy_package.nc (from hierarchy_compute.jl).
include(joinpath(@__DIR__, "patched_common.jl"))

# The bounds are for the L2-normalized ground state, the package uses ψ00 = 1.
# Rescale by Ψ(0) of the normalized Codex witness (Ψ(0) = 1/‖v‖₂ with v(0) = 1).
Wit = load_witness()
psi0 = 1 / sqrt(witness_norm2(Wit))
@printf("normalized Psi(0) = %.5f\n", psi0)

coeffs = [
    (0, 0, (a, t) -> psi00(a, t)),
    (1, 0, (a, t) -> psi10(a, t; Z = HE.Z)),
    (2, 0, (a, t) -> psi20(a, t; HE...)),
    (2, 1, (a, t) -> psi21(a, t; Z = HE.Z)),
    (3, 0, (a, t) -> psi30(a, t; HE..., rtol = 1e-12)),
    (3, 1, (a, t) -> psi31(a, t; Z = HE.Z)),
    (4, 1, (a, t) -> psi41(a, t; HE...)),
    (4, 2, (a, t) -> psi42(a, t; Z = HE.Z)),
]
αs = interior(41; ε = 1e-2)
θs = interior(41; ε = 1e-2)
fields = Dict{Tuple{Int,Int},Matrix{Float64}}()
for (k, p, f) in coeffs
    t = @elapsed fields[(k, p)] = grid(f, αs, θs)
    @printf("psi%d%d: max |.| = %.4f  (%.1f s)\n", k, p, maximum(abs, fields[(k, p)]), t)
end
# ψ40 costs ~0.3 s per point at step = 0.125; use a coarse grid. Its free harmonic
# part (a40, a42) is not fixed by the recurrence; it is set to zero here.
αc = interior(15; ε = 0.05)
θc = interior(15; ε = 0.05)
t40 = @elapsed psi40grid = grid((a, t) -> psi40(a, t; HE..., step = 0.125), αc, θc)
@printf("psi40: max |.| = %.4f  (%.1f s, 15x15 grid)\n", maximum(abs, psi40grid), t40)
maxabs = Dict(kp => psi0 * maximum(abs, A) for (kp, A) in fields)
maxabs[(4, 0)] = psi0 * maximum(abs, psi40grid)

ns = 0:8
mean(x) = sum(x) / length(x)
log2_h_bound(n) = 54 + 12n                               # ‖h_n‖_{H1} ≤ 2^54 4096^n
log2_p_bound(n) = 155 + 8log2(n + 1) + 13n               # ‖p_n‖_{C0} ≤ 2^155 (n+1)^8 8192^n
log2_p_sharp(n) = 153 + 8log2(n + 1) + 12n               # 2026-10-14: 2^153 (n+1)^8 4096^n
emp(n) = maximum(v for ((k, p), v) in maxabs if k == n)

# Empirical geometric rate from k = 1..4 (least squares on log2 max|ψ_k|).
kk = 1:4
y = [log2(emp(k)) for k in kk]
slope = sum((kk .- mean(kk)) .* (y .- mean(y))) / sum((kk .- mean(kk)) .^ 2)
@printf(
    "empirical growth rate 2^%.3f per order, i.e. radius scale %.3f bohr\n",
    slope,
    2.0^-slope
)

# ---------------------------------------------------------------- figure
fig = Figure(size = (1350, 600))
ax = Axis(
    fig[1, 1];
    xlabel = L"order $n$ (power of $\rho$)",
    ylabel = L"\log_2\,(\mathrm{norm})",
    title = "(a) proven growth bounds vs actual coefficient sizes",
    titlefont = :regular,
    xticks = ns,
)
lines!(
    ax,
    ns,
    log2_p_bound.(ns);
    color = C.blue,
    linewidth = 2,
    label = L"$\|p_n\|_{C^0} \le 2^{155}(n+1)^8\,8192^n$ (radius $1/16384$)",
)
lines!(
    ax,
    ns,
    log2_p_sharp.(ns);
    color = C.blue,
    linewidth = 2,
    linestyle = :dash,
    label = L"$\|p_n\|_{C^0} \le 2^{153}(n+1)^8\,4096^n$ (sharpened, radius $1/8192$)",
)
lines!(
    ax,
    ns,
    log2_h_bound.(ns);
    color = C.orange,
    linewidth = 2,
    label = L"$\|h_n\|_{H^1} \le 2^{54}\,4096^n$ (free harmonic data)",
)
markers = (:circle, :rect, :diamond)
for p = 0:2
    pts = [(k, log2(v)) for ((k, q), v) in maxabs if q == p]
    sort!(pts)
    scatter!(
        ax,
        first.(pts),
        last.(pts);
        color = C.aqua,
        marker = markers[p+1],
        markersize = 11,
        strokecolor = :white,
        strokewidth = 1,
        label = L"FockExpansion.jl $\max_\Omega |\psi_{n,%$p}|$",
    )
end
text!(
    ax,
    4,
    log2(emp(4));
    text = L"\approx 2^{%$(round(Int, log2_p_bound(4) - log2(emp(4))))}\ \mathrm{gap\ at}\ n=4",
    align = (:center, :bottom),
    offset = (0, 14),
    color = C.ink2,
    fontsize = 14,
)
ylims!(ax, -12, 280)
axislegend(ax; position = :lt, framevisible = false, labelsize = 13, rowgap = 0)

# (b) term sizes at the patch radius
rstar = 2.0^-14
ax2 = Axis(
    fig[1, 2];
    xlabel = L"order $n$",
    ylabel = L"\log_2\,(\rho^n \times \mathrm{norm})",
    title = L"(b) size of the $n$-th Fock term at $\rho = r_* = 2^{-14}$ and at $\rho = 1/2$",
    xticks = ns,
)
lines!(
    ax2,
    ns,
    log2_p_bound.(ns) .- 14ns;
    color = C.blue,
    linewidth = 2,
    label = L"proven bound at $\rho=r_*$",
)
ks = 0:4
scatterlines!(
    ax2,
    ks,
    [log2(emp(k)) - 14k for k in ks];
    color = C.aqua,
    linewidth = 2,
    markersize = 10,
    label = L"actual $\max_{p,\Omega}|\psi_{n,p}|\,r_*^n$",
)
scatterlines!(
    ax2,
    ks,
    [log2(emp(k)) - k for k in ks];
    color = C.violet,
    linewidth = 2,
    markersize = 10,
    label = L"actual at $\rho = 1/2$ (log factors omitted)",
)
hlines!(ax2, [-52]; color = C.gray, linestyle = :dot)
text!(
    ax2,
    8,
    -52;
    text = L"2^{-52}",
    align = (:right, :bottom),
    color = C.ink2,
    fontsize = 13,
)
ylims!(ax2, -70, 175)
axislegend(ax2; position = :rc, framevisible = false, labelsize = 13)
Label(
    fig[2, 1:2],
    @sprintf(
        "Coefficients scaled to the L²-normalized state by Ψ(0) = %.4f (Codex witness). Norms differ:\n the bounds use C⁰,* / H¹, the points are maxima over a 41×41 interior grid (ψ₄₀: 15×15, harmonic part a₄₀ = a₄₂ = 0). Empirical growth ≈ 2^(%.2f) per order (proven: 2^13).",
        psi0,
        slope
    );
    fontsize = 13,
    color = C.ink2,
    tellwidth = false,
)
savefig("patched_growth", fig)

# ---------------------------------------------------------------- NetCDF
ncgroup(
    "growth";
    attrib = [
        "comment" => "Proven growth bounds (EFFECTIVE_UNIFORM_FOCK, EFFECTIVE_FOCK_AVERAGING) vs FockExpansion.jl coefficients, scaled to the L2-normalized state",
    ],
) do g
    putvar!(g, "n", collect(ns), ("n",); attrib = ["long_name" => "coefficient index n"])
    putvar!(
        g,
        "log2_p_bound",
        log2_p_bound.(ns),
        ("n",);
        attrib = [
            "long_name" => "2^{155} (n+1)^{8} 8192^{n}",
            "H5WEB_label" => "proven bound",
            "H5WEB_color" => "#eb6834",
            view_hint(
                vis = "line",
                value_transform = "exp2",
                series = "/growth/log2_p_bound_sharp /growth/log2_h_bound",
            )...,
        ],
    )
    putvar!(
        g,
        "log2_p_bound_sharp",
        log2_p_sharp.(ns),
        ("n",);
        attrib = [
            "long_name" => "2^{153} (n+1)^{8} 4096^{n}",
            "H5WEB_label" => "sharpened bound",
            "H5WEB_color" => "#2a78d6",
            view_hint(value_transform = "exp2")...,
        ],
    )
    putvar!(
        g,
        "log2_h_bound",
        Float64.(log2_h_bound.(ns)),
        ("n",);
        attrib = [
            "long_name" => "2^{54} 4096^{n}",
            "H5WEB_label" => "H^{1} bound",
            "H5WEB_color" => "#1baf7a",
            view_hint(value_transform = "exp2")...,
        ],
    )
    M = fill(NaN, 5, 3)
    for ((k, p), v) in maxabs
        M[k+1, p+1] = v
    end
    putvar!(g, "k", collect(0:4), ("k",); attrib = ["long_name" => "power of R"])
    putvar!(g, "p", collect(0:2), ("p",); attrib = ["long_name" => "power of log R"])
    putvar!(
        g,
        "max_abs_psi",
        M,
        ("k", "p");
        attrib = ["long_name" => "Ψ(0) max over grid |ψ_{k,p}| (NaN: absent)"],
    )
end
