# Error budgets of the patched construction, re-derived from the formulas in
# EFFECTIVE_FOCK_AVERAGING / EFFECTIVE_UNIFORM_FOCK / EFFECTIVE_HEAT_SELECTION and
# checked against the JSON files written by Codex's budget scripts (read only).
#
#     julia -t 4 --project=examples/sciviz examples/sciviz/patched_budget.jl
include(joinpath(@__DIR__, "patched_common.jl"))

# ---------------------------------------------------------------- (a) residual
# Coefficient n to H1 accuracy ε needs full residual δ ≤ ε 2^-(56+12n).
log2_delta(n, k) = -(k + 56 + 12n)              # ε = 2^-k
fock_json = read(joinpath(BM, "EFFECTIVE_FOCK_BUDGET.json"), String)
for (n, r) in zip(
    parse.(Int, jsonall(fock_json, "n")),
    jsonall(fock_json, "sufficient_full_residual"),
)
    @assert r == "2^-$(-log2_delta(n, 20))" (n, r)
end

seed_json = read(joinpath(BM, "EFFECTIVE_SEED_CERTIFICATE.json"), String)
seed_res = parse(Float64, jsonstring(seed_json, "full_residual_norm_upper"))
E_lo = parse(BigFloat, jsonstring(seed_json, "energy_lower"))
E_hi = parse(BigFloat, jsonstring(seed_json, "energy_upper"))
E0 = parse(BigFloat, "-2.903724377034119598311159245194404446696925310")
@assert E_lo < E0 < E_hi

# ---------------------------------------------------------------- (b) Fock tail
log2_tail(N) = 160 + 8log2(N + 2) - N           # 2^160 (N+2)^8 2^-N, N ≥ 31
function select_N(k)
    N = 31
    while N < 160 + k || big(N + 2)^8 > big(2)^(N - 160 - k)
        N += 1
    end
    N
end
unif_json = read(joinpath(BM, "EFFECTIVE_UNIFORM_FOCK_BUDGET.json"), String)
targets = [parse(Int, split(s, "^-")[2]) for s in jsonall(unif_json, "target")]
Njson = parse.(Int, jsonall(unif_json, "N"))
Nsel = select_N.(targets)
@assert Nsel == Njson (Nsel, Njson)

# ---------------------------------------------------------------- (c) heat schedule
cdiv(x) = cld(numerator(x), denominator(x))
function log2ceil(x)
    x <= 1 && return 0
    k = 0
    while big(2)^k < x
        k += 1
    end
    k
end
function sqrtceil(x)
    k = isqrt(fld(numerator(x), denominator(x)))
    while k^2 < x
        k += 1
    end
    k
end
"""EFFECTIVE_HEAT_BUDGET.schedule: V bounds ‖v‖₂ and ‖v‖∞, G both gradient norms."""
function schedule(ε; V = big(10) // 1, G = big(100) // 1)
    t = big(1) // 1024
    M = cdiv(36V / ε)
    a = 5M * t
    J = max(cdiv(6a), log2ceil(4V / ε))
    K = log2ceil(48(J + 1) * V / ε)
    L = sqrtceil(4a + 2K)
    Cc = (J + 1) * (4L + 3) * (t * 6M^2 * V + G)
    m = cdiv(2a) + log2ceil(4Cc / ε)
    (M = M, J = J, L = L, log2_inv_h = 2m, dim = 7J + 6)
end
heat_json = read(joinpath(BM, "EFFECTIVE_HEAT_BUDGET.json"), String)
heat_k =
    [round(Int, log2(parse(Int, split(s, "/")[2]))) for s in jsonall(heat_json, "epsilon")]
for (i, k) in enumerate(heat_k)
    s = schedule(big(1) // big(2)^k)
    @assert s.M == parse(Int, jsonall(heat_json, "clipping_M")[i])
    @assert s.J == parse(Int, jsonall(heat_json, "Dyson_order_J")[i])
    @assert s.L == parse(Int, jsonall(heat_json, "Gaussian_halfwidth_L")[i])
    @assert "2^(-$(s.log2_inv_h))" == jsonall(heat_json, "mesh_h")[i]
    @assert s.dim == parse(Int, jsonall(heat_json, "max_integration_dimension")[i])
end
ks = 1:14
sched = [schedule(big(1) // big(2)^k) for k in ks]
println("all budget numbers reproduce the Codex JSON files")

# ---------------------------------------------------------------- figure
fig = Figure(size = (1350, 1000))

# (a) required residual
ax = Axis(
    fig[1, 1];
    xlabel = L"coefficient index $n$",
    ylabel = L"\log_2 \delta",
    title = L"(a) full PDE residual $\delta$ sufficient for $\|h_n - h_n^{\mathrm{exact}}\|_{H^1} \le \varepsilon$",
    xticks = 0:2:16,
)
ns = 0:16
for (k, col) in ((10, C.aqua), (20, C.blue), (40, C.violet))
    lines!(ax, ns, log2_delta.(ns, k); color = col, linewidth = 2)
    text!(
        ax,
        16,
        log2_delta(16, k);
        text = L"\varepsilon = 2^{-%$k}",
        color = C.ink2,
        align = (:left, :center),
        offset = (6, 0),
        fontsize = 14,
    )
end
jn = parse.(Int, jsonall(fock_json, "n"))
scatter!(
    ax,
    jn,
    log2_delta.(jn, 20);
    color = C.blue,
    markersize = 11,
    strokecolor = :white,
    strokewidth = 1.5,
    label = L"EFFECTIVE\_FOCK\_BUDGET.json ($\varepsilon=2^{-20}$)",
)
hlines!(ax, [log2(seed_res)]; color = C.orange, linewidth = 2, linestyle = :dash)
text!(
    ax,
    0.2,
    log2(seed_res);
    text = L"certified seed residual $2.42\times10^{-4} \approx 2^{-12}$",
    color = C.ink2,
    align = (:left, :bottom),
    offset = (0, 3),
    fontsize = 14,
)
hlines!(ax, [-52]; color = C.gray, linewidth = 1.5, linestyle = :dot)
text!(
    ax,
    18.4,
    -52;
    text = L"Float64 $\epsilon = 2^{-52}$",
    color = C.ink2,
    align = (:right, :bottom),
    offset = (0, 3),
    fontsize = 14,
)
xlims!(ax, -0.6, 18.6)
ylims!(ax, -300, 5)
arrows2d!(
    ax,
    [1.0],
    [log2(seed_res) - 2],
    [0.0],
    [log2_delta(0, 20) - log2(seed_res) + 6];
    color = C.ink2,
    shaftwidth = 1.5,
    tipwidth = 9,
    tiplength = 9,
)
text!(
    ax,
    1.1,
    -45;
    text = L"$2^{64}$ smaller than the seed residual (already at $n=0$)",
    color = C.ink2,
    align = (:left, :center),
    fontsize = 14,
)
axislegend(ax; position = :lb, framevisible = false, labelsize = 14)

# (b) uniform tail
ax = Axis(
    fig[1, 2];
    xlabel = L"truncation order $N$",
    ylabel = L"\log_2 \,(\mathrm{tail\ bound})",
    title = L"(b) all-angle Fock tail on $\rho\le r_*$: $2^{160}(N+2)^8\,2^{-N}$",
)
Ns = 31:420
lines!(ax, Ns, log2_tail.(Ns); color = C.blue, linewidth = 2)
for (k, N) in zip(targets, Nsel)
    hlines!(ax, [-k]; color = C.gray, linewidth = 1, linestyle = :dot)
    scatter!(
        ax,
        [N],
        [log2_tail(N)];
        color = C.orange,
        markersize = 11,
        strokecolor = :white,
        strokewidth = 1.5,
    )
    right = N < 380
    text!(
        ax,
        N,
        log2_tail(N);
        text = L"N=%$N",
        color = C.ink2,
        align = (right ? :left : :right, :bottom),
        offset = (right ? 8 : -8, 2),
        fontsize = 14,
    )
end
text!(
    ax,
    40,
    -160;
    text = "dotted: targets 2⁰, 2⁻²⁰, 2⁻⁴⁰, 2⁻⁸⁰, 2⁻¹⁶⁰",
    color = C.ink2,
    align = (:left, :bottom),
    offset = (0, 3),
    fontsize = 13,
)
text!(
    ax,
    40,
    2;
    text = "bound exceeds 1 for N < 223",
    color = C.ink2,
    align = (:left, :bottom),
    fontsize = 13,
)
ylims!(ax, -185, 205)

# (c) heat schedule
ax = Axis(
    fig[2, 1];
    xscale = log2,
    yscale = log10,
    xreversed = true,
    xlabel = L"global sup/$L^2$ target $\varepsilon$",
    ylabel = "parameter value",
    title = L"(c) heat-smoothing budget ($t_*=1/1024$, $\|v\|\le 10$, $\|\nabla v\|\le 100$)",
    xticks = (2.0 .^ -(2:2:14), [L"2^{-%$k}" for k = 2:2:14]),
)
εs = 2.0 .^ -ks
series = [
    (:M, L"clipping level $M$", C.blue),
    (:J, L"Dyson order $J$", C.orange),
    (:L, L"Gaussian half-width $L$", C.aqua),
    (:log2_inv_h, L"$\log_2(1/h)$ (cubature mesh $h$)", C.violet),
    (:dim, L"integration dimension $7J+6$", C.magenta),
]
for (key, lab, col) in series
    y = [Float64(getproperty(s, key)) for s in sched]
    scatterlines!(ax, εs, y; color = col, linewidth = 2, markersize = 7, label = lab)
end
for (i, k) in enumerate(heat_k)
    scatter!(
        ax,
        fill(2.0^-k, 5),
        [Float64(getproperty(schedule(big(1) // big(2)^k), s[1])) for s in series];
        color = :transparent,
        strokecolor = C.ink,
        strokewidth = 1.5,
        markersize = 13,
    )
end
text!(
    ax,
    2.0^-14,
    3;
    text = "open circles: rows of EFFECTIVE_HEAT_BUDGET.json",
    color = C.ink2,
    align = (:right, :bottom),
    fontsize = 13,
)
axislegend(ax; position = :lt, framevisible = false, labelsize = 14)

# (d) energy enclosure
g = fig[2, 2] = GridLayout()
Label(
    g[0, 1],
    L"(d) Temple enclosure of $E_0$ from the certified seed ($\beta = -5/2$)";
    fontsize = 16,
    font = :bold,
)
unit = 1e-7
lo, hi = Float64(E_lo - E0), Float64(E_hi - E0)
axd = Axis(
    g[1, 1];
    xlabel = L"E - E_0 \;\; (10^{-7}\ \mathrm{hartree})",
    yticksvisible = false,
    yticklabelsvisible = false,
    height = 120,
)
band!(axd, [lo, hi] ./ unit, [-1, -1], [1, 1]; color = (C.blue, 0.25))
lines!(axd, [lo, hi] ./ unit, [0, 0]; color = C.blue, linewidth = 4)
vlines!(axd, [0]; color = C.orange, linewidth = 2)
text!(
    axd,
    lo / unit,
    1;
    text = @sprintf("lower %.4e", lo),
    align = (:left, :top),
    offset = (4, -4),
    color = C.ink2,
    fontsize = 13,
)
text!(
    axd,
    0,
    -1;
    text = L"$E_0$ (high-precision)",
    align = (:right, :bottom),
    offset = (-6, 4),
    color = C.ink2,
    fontsize = 13,
)
xlims!(axd, -1.6, 0.15)
ylims!(axd, -1.2, 1.2)
axz = Axis(
    g[2, 1];
    xlabel = L"E - E_0 \;\; (10^{-10}\ \mathrm{hartree})",
    yticksvisible = false,
    yticklabelsvisible = false,
    height = 120,
    title = "zoom on the upper end",
    titlesize = 14,
    titlefont = :regular,
)
band!(axz, [-10, hi / 1e-10], [-1, -1], [1, 1]; color = (C.blue, 0.25))
lines!(axz, [-10, hi / 1e-10], [0, 0]; color = C.blue, linewidth = 4)
vlines!(axz, [0]; color = C.orange, linewidth = 2)
text!(
    axz,
    hi / 1e-10,
    1;
    text = @sprintf("upper %.3e\n(Rayleigh quotient)", hi),
    align = (:right, :top),
    offset = (-4, -4),
    color = C.ink2,
    fontsize = 13,
)
xlims!(axz, -1, 4)
ylims!(axz, -1.2, 1.2)
Label(
    g[3, 1],
    @sprintf(
        "width %.3e;  E₀ sits %.1f%% from the upper end",
        hi - lo,
        100 * hi / (hi - lo)
    );
    fontsize = 14,
    color = C.ink2,
)
rowgap!(g, 8)

savefig("patched_budget", fig)

# ---------------------------------------------------------------- NetCDF
ncgroup(
    "budget";
    attrib = [
        "comment" => "Error budgets of the EFFECTIVE_* construction, recomputed in Julia; budgets only, no cubature executed",
    ],
) do grp
    putvar!(grp, "n", collect(ns), ("n",); attrib = ["long_name" => "coefficient index n"])
    putvar!(
        grp,
        "log2_delta_eps20",
        Float64.(log2_delta.(ns, 20)),
        ("n",);
        attrib = [
            "long_name" => "residual for H^{1} error 2^{-20}",
            view_hint(vis = "line", value_transform = "exp2")...,
        ],
    )
    putvar!(
        grp,
        "N",
        collect(Ns),
        ("N",);
        attrib = ["long_name" => "Fock truncation order N"],
    )
    putvar!(
        grp,
        "log2_tail",
        log2_tail.(Ns),
        ("N",);
        attrib = [
            "long_name" => "2^{160} (N+2)^{8} 2^{-N}",
            view_hint(vis = "line", value_transform = "exp2")...,
        ],
    )
    putvar!(
        grp,
        "log2_inv_eps",
        collect(ks),
        ("log2_inv_eps",);
        attrib = ["long_name" => "-log2 epsilon"],
    )
    heat_names = Dict(
        :M => "clipping level M",
        :J => "Dyson order J",
        :L => "Gaussian half-width L",
        :log2_inv_h => "log2(1/h), cubature mesh h",
        :dim => "integration dimension 7J+6",
    )
    for (key, _, _) in series
        putvar!(
            grp,
            "heat_" * string(key),
            [Float64(getproperty(s, key)) for s in sched],
            ("log2_inv_eps",);
            attrib = [
                "long_name" => heat_names[key],
                view_hint(vis = "line", scale = "log")...,
            ],
        )
    end
    putvar!(
        grp,
        "energy_interval_minus_E0",
        [lo, hi],
        ("end",);
        attrib = [
            "units" => "hartree",
            "long_name" => "Temple interval minus E_{0}",
            "E0" => Float64(E0),
        ],
    )
end
