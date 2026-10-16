# Codex analytic series for ψ₃₀ (outputs/psi30_analytic in the Codex workspace):
# partial-wave and Gauss-series convergence, error maps, the certified R1 recurrence,
# the assembly diagram, the three-evaluator comparison, and output/psi30_codex.nc.
# Needs the CSVs from codex_series_data.py and codex_green.jl:
#     julia -t 4 --project=examples/sciviz examples/sciviz/codex_series.jl
include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "codex_util.jl"))
using Printf

const SDIR = joinpath(OUT, "codex_series")
const GDIR = joinpath(OUT, "codex_green")
const PAR = (Z = HE.Z, E = HE.E, a21 = HE.a21)
reference(α, θ) = psi30(α, θ; PAR..., rtol = 1e-13)
pretty(name) = replace(name, "_" => " ", "r12" => "r₁₂", "r1 " => "r₁ ")
const COMPONENTS = [("known", L"Z^2(F_{2c}+F_{2e})"), ("d2", L"Z^2 d_2"), ("d1", L"Z d_1")]

# Convergence factors: the nuclear-frame waves decay like tan(x/2)^l with
# x = min(α, π-α); the e-e-frame waves (d₁) like tan(y/2)^l, cos y = -sin α cos θ.
t_nuc(α, θ) = tan(min(α, π - α) / 2)
function t_ee(α, θ)
    y = acos(clamp(-sin(α) * cos(θ), -1, 1))
    tan(min(y, π - y) / 2)
end

"""Per-l terms of each component at each probe: Dict point => (α, θ, base, Dict comp => terms[l+1])."""
function probe_terms(file = "probe_terms.csv"; key = "point")
    d = read_csv(joinpath(SDIR, file))
    out = Dict{String,Any}()
    for k in eachindex(d[key])
        p = string(d[key][k])
        α = haskey(d, "alpha") ? d["alpha"][k] : π / 2
        θ = haskey(d, "alpha") ? d["theta"][k] : parse(Float64, p)
        e = get!(out, p, (α = α, θ = θ, base = Ref(0.0), terms = Dict(c => zeros(0) for (c, _) in COMPONENTS)))
        c, l = d["component"][k], Int(d["l"][k])
        if c == "base"
            e.base[] = d["term"][k]
        else
            v = e.terms[c]
            length(v) < l + 1 && resize!(v, l + 1)
            v[l+1] = d["term"][k]
        end
    end
    out
end
partial_sums(e) = e.base[] .+ cumsum(sum(e.terms[c] for (c, _) in COMPONENTS))

function fig_convergence(P, line)
    fig = Figure(size = (1150, 480))
    ax1 = Axis(fig[1, 1]; yscale = log10, xlabel = L"truncation $L$", ylabel = L"|\psi_{30}^{(L)} - \psi_{30}|",
        title = L"\mathrm{probe\ points}", yticks = LogTicks(-16:4:0))
    order = ["generic", "collinear", "near_r1_0.05", "equator", "near_r12_0.30", "near_r12_0.10", "near_r12_0.03"]
    for (k, p) in enumerate(order)
        e = P[p]
        err = abs.(partial_sums(e) .- reference(e.α, e.θ))
        lines!(ax1, 0:length(err)-1, max.(err, 1e-17); color = Makie.wong_colors()[k], label = pretty(p))
        ρ = max(t_nuc(e.α, e.θ), t_ee(e.α, e.θ))
        ρ < 0.9 && lines!(ax1, 10:80, 0.3err[11] .* ρ .^ (0:70); color = Makie.wong_colors()[k], linestyle = :dot)
    end
    ylims!(ax1, 1e-17, 1)
    Legend(fig[1, 2], ax1; framevisible = false, labelsize = 14)
    ax2 = Axis(fig[1, 3]; xscale = log10, yscale = log10, xlabel = L"truncation $L$",
        ylabel = L"\max_{L' \geq L} |\psi_{30}^{(L')} - \psi_{30}|", title = L"the line $\alpha = \pi/2$ ($\tan(x/2)=1$)", xticks = [1, 3, 10, 30, 100, 200])
    θs = sort(collect(keys(line)); by = s -> parse(Float64, s))
    cm = cgrad(:plasma, length(θs); categorical = true)
    for (k, s) in enumerate(θs)
        e = line[s]
        err = abs.(partial_sums(e) .- reference(e.α, e.θ))
        env = reverse(accumulate(max, reverse(err)))  # monotone envelope of the oscillating error
        lines!(ax2, 1:length(err)-1, env[2:end]; color = cm[k], label = latexstring("\\theta = $s"))
    end
    Ls = 20:200
    lines!(ax2, Ls, 2.0 .* (Ls ./ 20.0) .^ -2 .* 1e-3; color = :black, linestyle = :dash, label = L"\propto L^{-2}")
    Legend(fig[1, 4], ax2; framevisible = false, labelsize = 14)
    Label(fig[2, 1:4], "Dotted: predicted geometric rate ρᴸ, ρ = max(tan(x/2), tan(y/2)) from the nuclear and e-e frames. Reference: elliptic psi30, rtol = 1e-13.";
        fontsize = 13, color = :gray30)
    savefig("codex_series_convergence", fig)
end

function fig_terms(P)
    pts = ["generic", "equator", "near_r12_0.10", "near_r1_0.05"]
    fig = Figure(size = (1200, 430))
    for (j, (c, lab)) in enumerate(COMPONENTS)
        ax = Axis(fig[1, j]; yscale = log10, xlabel = L"partial wave $l$", ylabel = j == 1 ? L"|\mathrm{term}_l|" : "",
            title = lab, yticks = LogTicks(-40:10:0))
        for (k, p) in enumerate(pts)
            e = P[p]
            v = abs.(e.terms[c])
            ls = c == "d1" ? (0:2:length(v)-1) : (0:length(v)-1)
            scatterlines!(ax, ls, max.(v[ls.+1], 1e-60); color = Makie.wong_colors()[k], markersize = 4,
                label = pretty(p))
            ρ = c == "d1" ? t_ee(e.α, e.θ) : t_nuc(e.α, e.θ)
            ρ < 0.95 && lines!(ax, 0:80, v[3] .* ρ .^ ((0:80) .- 2); color = Makie.wong_colors()[k], linestyle = :dot)
        end
        ylims!(ax, 1e-45, 1)
        j == 3 && axislegend(ax; position = :lb, framevisible = false, labelsize = 13)
    end
    Label(fig[2, 1:3], "Terms of the three Legendre series (Pₗ factor included). Dotted: tan(x/2)ˡ, x the polar angle of the frame (nuclear for F₂c, F₂e, d₂; e-e for d₁).";
        fontsize = 13, color = :gray30)
    savefig("codex_series_terms", fig)
end

function fig_inner_N(P)
    d = read_csv(joinpath(SDIR, "probe_N.csv"))
    fig = Figure(size = (720, 480))
    ax = Axis(fig[1, 1]; yscale = log10, xlabel = L"Gauss-series length $N$ (at $L = 40$)",
        ylabel = L"|\psi_{30}^{(N)} - \psi_{30}^{(300)}|", title = L"inner series in $w = \sin(x/2) \leq 1/\sqrt{2}$")
    for (k, p) in enumerate(unique(d["point"]))
        s = findall(==(p), d["point"])
        Ns, v = d["N"][s], d["psi30"][s]
        sel = 1:length(s)-1
        scatterlines!(ax, Ns[sel], max.(abs.(v[sel] .- v[end]), 1e-17); color = Makie.wong_colors()[k], label = pretty(p))
    end
    lines!(ax, 10:200, 0.1 .* 2.0 .^ (-(10:200) ./ 2); color = :black, linestyle = :dash, label = L"2^{-N/2}")
    ylims!(ax, 1e-17, 1e-1)
    axislegend(ax; position = :rt, framevisible = false, labelsize = 13)
    savefig("codex_series_inner_N", fig)
end

function grid_maps()
    d = read_csv(joinpath(SDIR, "grid_partial_sums.csv"))
    n = Int(maximum(d["ialpha"])) + 1
    LG = Int(maximum(d["L"]))
    S = zeros(n, n, LG + 1)
    for k in eachindex(d["L"])
        S[Int(d["ialpha"][k])+1, Int(d["itheta"][k])+1, Int(d["L"][k])+1] = d["psi30_L"][k]
    end
    αs = interior(n)
    ref = grid(reference, αs, αs)
    αs, S, ref
end

function fig_error_maps(αs, S, ref)
    fig = Figure(size = (1200, 800))
    hm = nothing
    for (k, L) in enumerate((2, 5, 10, 20, 40))
        ax = Axis(fig[(k-1)÷3+1, (k-1)%3+1]; title = latexstring("L = $L"), aspect = DataAspect(),
            xlabel = L"\alpha", ylabel = L"\theta")
        hm = heatmap!(ax, αs, αs, log10.(abs.(S[:, :, L+1] .- ref) .+ 1e-17); colormap = :inferno, colorrange = (-15, -1))
    end
    Colorbar(fig[1:2, 4], hm; label = L"\log_{10}|\psi_{30}^{(L)} - \psi_{30}|")
    ax = Axis(fig[2, 3]; title = L"predicted $\log_{10}\rho^{41}$, $\rho=\max(\tan\frac{x}{2},\tan\frac{y}{2})$",
        aspect = DataAspect(), xlabel = L"\alpha", ylabel = L"\theta")
    ρ = [max(t_nuc(a, t), t_ee(a, t)) for a in αs, t in αs]
    heatmap!(ax, αs, αs, 41 .* log10.(ρ); colormap = :inferno, colorrange = (-15, -1))
    contour!(ax, αs, αs, ρ; levels = [0.5, 0.8], color = :white, linewidth = 0.8)
    savefig("codex_series_error_maps", fig)
    ρ
end

function fig_r1()
    d = read_csv(joinpath(SDIR, "r1.csv"))
    l = Int.(d["l"])
    fig = Figure(size = (1200, 430))
    ax1 = Axis(fig[1, 1]; xlabel = L"l", ylabel = L"R_1(l)", title = L"$R_1(l)$ from the step-two recurrence")
    for (par, c, lab) in ((0, Makie.wong_colors()[1], L"l\ \mathrm{even}"), (1, Makie.wong_colors()[2], L"l\ \mathrm{odd}"))
        s = findall(x -> x % 2 == par && x <= 60, l)
        scatterlines!(ax1, l[s], d["R1"][s]; color = c, markersize = 6, label = lab)
    end
    text!(ax1, [Point2f(18, 0.12 - 0.07k) for k = 1:4]; text = [latexstring("R_1($(l[k])) = $(Int(d["num"][k]) < 0 ? "-" : "")\\frac{$(abs(Int(d["num"][k])))}{$(Int(d["den"][k]))}") for k = 1:4],
        fontsize = 15, align = (:left, :center))
    axislegend(ax1; position = :rt, framevisible = false)
    ax2 = Axis(fig[1, 2]; xscale = log10, yscale = log10, xlabel = L"l", ylabel = L"|R_1(l+2) - R_1(l)|",
        title = L"first difference $d_l$", xticks = [3, 10, 30, 100])
    scatter!(ax2, l, abs.(d["d_l"]); markersize = 5, color = Makie.wong_colors()[3])
    lines!(ax2, 10:200, 1.8 .* (10:200) .^ -1.5; color = :black, linestyle = :dash)
    text!(ax2, 30, 0.045; text = L"\propto l^{-3/2}", fontsize = 15)
    # d_l ~ l^(-3/2), so R1(l) approaches its limit like l^(-1/2); the limit is -log 2 numerically.
    ax3 = Axis(fig[1, 3]; xscale = log10, yscale = log10, xlabel = L"l", ylabel = L"R_1(l) + \log 2",
        title = L"limit $R_1(\infty) = -\log 2$ (numerical)", xticks = [3, 10, 30, 100])
    scatter!(ax3, l, d["R1"] .+ log(2); markersize = 5, color = Makie.wong_colors()[4])
    lines!(ax3, 10:200, 1.1 .* (10:200) .^ -0.5; color = :black, linestyle = :dash)
    text!(ax3, 20, 0.14; text = L"\propto l^{-1/2}", fontsize = 15)
    Label(fig[2, 1:3], "Exact rationals from the certified recurrence (R₁(l)+log 2 > 0 for all l ≤ 200); a 25-digit quadrature of C₁(c) agrees with every printed digit (20) for l = 2…60.";
        fontsize = 13, color = :gray30)
    savefig("codex_r1_recurrence", fig)
end

function fig_assembly()
    fig = Figure(size = (1300, 700))
    ax = Axis(fig[1, 1]; limits = (0, 12.6, 0, 6.6), aspect = DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    col = Dict(:elem => "#cfe3f3", :beta => "#fbe3c3", :gauss => "#d6efd1", :cert => "#ead7f0", :out => "#f3d0d0")
    # A box with one or more centred lines of text.
    function box!(x, y, w, h, kind, lines; fs = 14)
        poly!(ax, Rect(x, y, w, h); color = col[kind], strokecolor = :gray30, strokewidth = 1)
        n = length(lines)
        for (k, txt) in enumerate(lines)
            text!(ax, x + w / 2, y + h / 2 + (n + 1 - 2k) * 0.2; text = txt, align = (:center, :center), fontsize = fs)
        end
    end
    arrow!(p, q) = arrows2d!(ax, [Point2f(p)], [Vec2f(q .- p)]; color = :gray30, shaftwidth = 1.5, tipwidth = 9, tiplength = 9)
    box!(0.2, 5.6, 12.2, 0.8, :out,
        [L"\psi_{30} = F_0 + Z\,(F_{1a} + F_{1b} + F_{1c} + d_1) + Z^2(F_{2a} + F_{2b} + F_{2c} + F_{2e} + d_2) + Z^3 F_3 + a\,r(\alpha,\theta)"]; fs = 17)
    box!(0.2, 3.6, 2.9, 1.3, :elem, [L"F_0, F_{1a}, F_{1b}, F_{1c}, F_{2a}, F_3", L"\mathrm{elementary\ in}\ \xi, \sin\alpha, \varsigma"])
    box!(0.2, 2.0, 2.9, 1.2, :elem, [L"r = \xi(6-5\xi^2)/12 - Z\varsigma v/2", L"a = a_{21} - Z(\pi+4)/(9\pi)"])
    box!(0.2, 0.3, 2.9, 1.2, :elem, [L"\mathrm{exact,\ no\ truncation}"])
    box!(3.5, 3.6, 2.9, 1.3, :beta, [L"F_{2c}, F_{2e} = \sum_l P_l(\cos\theta)\, f_l(t)", L"F_{2b} = \cos\theta\, f_1(t)"])
    box!(3.5, 2.0, 2.9, 1.2, :beta, [L"f_l:\ \mathrm{finite\ sums\ of}\ B_{t^2/(1+t^2)}(a,b)", L"t = \tan(x/2),\ x = \min(\alpha, \pi-\alpha)"]; fs = 13)
    box!(3.5, 0.3, 2.9, 1.2, :beta, [L"\mathrm{terms} \sim \tan(x/2)^l", L"\mathrm{algebraic\ on}\ \alpha = \pi/2"])
    box!(6.8, 3.6, 2.6, 1.3, :gauss, [L"d_2 = \sum_l P_l(\cos\theta)\, f_l(x)", L"\mathrm{nuclear\ frame}"])
    box!(9.7, 3.6, 2.7, 1.3, :gauss, [L"d_1 = 2^{-1/2} \sum_{l\ \mathrm{even}} P_l(\mu_{ee})\, f_l(y)", L"\mathrm{e{-}e\ frame},\ \cos y = -\sin\alpha\cos\theta"]; fs = 13)
    box!(6.8, 2.0, 5.6, 1.2, :gauss, [L"f_l = P_l(z)\,[\,C_l\, Q_l(z) + \sum_{e \leq N/2} g_e z^e\,],\quad z = \sin^2(x/2)",
        L"g_e:\ \mathrm{two{-}term\ recurrence\ (Gauss\ series,\ ratio}\ \leq 1/2)"])
    box!(6.8, 0.3, 2.6, 1.2, :cert, [L"C_0, C_1:\ \mathrm{closed\ constants}", L"(\mathrm{Catalan}\ G,\ \pi,\ \log 2)"])
    box!(9.7, 0.3, 2.7, 1.2, :cert, [L"C_l,\ l \geq 2:\ J{-}\mathrm{multipoles}", L"\mathrm{certified\ step{-}two\ recurrences}"])
    for x in (1.65, 4.95, 8.1, 11.05)
        arrow!((x, 4.9), (x, 5.6))
        arrow!((x, 3.2), (x, 3.6))
    end
    for x in (4.95, 8.1, 11.05)
        x == 4.95 || arrow!((x, 1.5), (x, 2.0))
    end
    Legend(fig[2, 1], [PolyElement(color = col[k], strokecolor = :gray30) for k in (:elem, :beta, :gauss, :cert)],
        ["elementary closed form", "incomplete-beta partial waves (truncated at L)", "Gauss series + Legendre sum (truncated at L, N)",
            "point-independent constants"]; orientation = :horizontal, framevisible = false, labelsize = 14)
    savefig("codex_series_assembly", fig)
end

function fig_comparison()
    probes = ["generic", "equator", "near_r12_0.10", "near_r12_0.03", "near_r1_0.05", "collinear"]
    fig = Figure(size = (900, 680))
    ax = Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "time per point (s)",
        ylabel = "error vs elliptic reference (median over 6 probes; bars: min to max)",
        title = L"three evaluators of $\psi_{30}$ (He ground state)", yticks = LogTicks(-16:2:0))
    summarize(errs) = (median(errs), minimum(errs), maximum(errs))
    median(v) = (s = sort(v); n = length(s); isodd(n) ? s[(n+1)÷2] : (s[n÷2] + s[n÷2+1]) / 2)
    draw!(t, errs, c, mk, lab) = begin
        m, lo, hi = summarize(max.(errs, 1e-16))
        rangebars!(ax, [t], [lo], [hi]; color = (c, 0.5), whiskerwidth = 6)
        scatter!(ax, [t], [m]; color = c, marker = mk, markersize = 13, label = lab)
    end
    e = read_csv(joinpath(GDIR, "elliptic_timing.csv"))
    for (k, r) in enumerate(unique(e["rtol"]))
        s = findall(==(r), e["rtol"])
        draw!(median(e["seconds"][s]), e["error"][s], Makie.wong_colors()[1], :circle, k == 1 ? "elliptic psi30 (rtol 1e-6, 1e-8, 1e-11)" : nothing)
    end
    g = read_csv(joinpath(GDIR, "rtol_sweep.csv"))
    for (m, c, mk) in (("rectangular", Makie.wong_colors()[2], :rect), ("duffy", Makie.wong_colors()[3], :utriangle))
        for (k, r) in enumerate((1e-4, 1e-6, 1e-8, 1e-10, 1e-12))
            s = findall(i -> g["method"][i] == m && g["rtol"][i] ≈ r && g["point"][i] in probes, eachindex(g["rtol"]))
            draw!(median(g["seconds"][s]), g["error"][s], c, mk, k == 1 ? "psi30_green :$m (rtol 1e-4…1e-12)" : nothing)
        end
    end
    P = probe_terms()
    tm = read_csv(joinpath(SDIR, "timing.csv"))
    warm = Dict(Int(tm["L"][k]) => tm["seconds"][k] for k in eachindex(tm["L"]) if tm["what"][k] == "warm")
    cold = Dict(Int(tm["L"][k]) => tm["seconds"][k] for k in eachindex(tm["L"]) if tm["what"][k] == "full_py_uncached")
    for (k, L) in enumerate(sort(collect(keys(warm))))
        errs = [abs(partial_sums(P[p])[L+1] - reference(P[p].α, P[p].θ)) for p in probes]
        draw!(warm[L], errs, Makie.wong_colors()[4], :diamond, k == 1 ? "Codex series, mpmath, cached constants (L = 5…80)" : nothing)
        text!(ax, warm[L], median(max.(errs, 1e-16)); text = latexstring("L=$L"), offset = (-8, 0), fontsize = 12, align = (:right, :center))
        haskey(cold, L) && draw!(cold[L], errs, Makie.wong_colors()[5], :star5, L == 10 ? "Codex series, single cold call full.psi30" : nothing)
    end
    Legend(fig[2, 1], ax; orientation = :horizontal, nbanks = 3, framevisible = false, labelsize = 13)
    savefig("codex_psi30_comparison", fig)
end

function write_netcdf(αs, S, ref, ρ)
    G = read_matrix_csv(joinpath(GDIR, "maps.csv"))
    # Exact agreement is written as log10(0) = -Inf (drawn in the NaN colour).
    logerr(A) = log10.(abs.(A .- ref))
    hint(cm) = view_hint(vis = "heatmap", scale = "linear", colormap = cm)
    fields = [
        "psi30" => ref,
        "log10_err_green_rectangular" => logerr(G["green_rectangular"]),
        "log10_err_green_duffy" => logerr(G["green_duffy"]),
        "green_ms_rectangular" => 1e3 .* G["time_rectangular"],
        "green_ms_duffy" => 1e3 .* G["time_duffy"],
        "log10_err_series_L10" => logerr(S[:, :, 11]),
        "log10_err_series_L40" => logerr(S[:, :, 41]),
        "series_rate_rho" => ρ,
    ]
    attribs = Dict{String,Any}(
        "psi30" => ["long_name" => "ψ_{3,0} (elliptic, rtol 1e-13)", view_hint(vis = "heatmap", colormap = "Viridis")...],
        "series_rate_rho" => ["long_name" => "predicted geometric rate max(tan(x/2), tan(y/2))"; hint("Inferno")],
    )
    names = Dict(
        "log10_err_green_rectangular" => "|ψ_{3,0}^{Green} − ψ_{3,0}|, rectangular rule",
        "log10_err_green_duffy" => "|ψ_{3,0}^{Green} − ψ_{3,0}|, Duffy rule",
        "green_ms_rectangular" => "Green integral time per point, rectangular rule",
        "green_ms_duffy" => "Green integral time per point, Duffy rule",
        "log10_err_series_L10" => "|ψ_{3,0}^{(10)} − ψ_{3,0}|, Codex series",
        "log10_err_series_L40" => "|ψ_{3,0}^{(40)} − ψ_{3,0}|, Codex series",
    )
    for (k, _) in fields[2:end-1]
        cm = startswith(k, "green_ms") ? "Cividis" : "Inferno"
        units = startswith(k, "green_ms") ? ["units" => "ms"] : Pair{String,String}[]
        logs = startswith(k, "log10_") ? view_hint(value_transform = "exp10") : Pair{String,Any}[]
        attribs[k] = ["long_name" => names[k]; units; hint(cm); logs]
    end
    path = joinpath(OUT, "psi30_codex.nc")
    write_angular(path, αs, αs, fields; attribs, global_attribs = ["title" => "ψ₃₀: Green integral and Codex analytic series vs the elliptic formula"])
    NCDataset(path, "a") do ds
        defVar(ds, "L", collect(0:size(S, 3)-1), ("L",); attrib = ["long_name" => "partial-wave truncation"])
        # Dims (alpha, theta, L) in Julia appear as (L, theta, alpha) in HDF5: L becomes the slider.
        E = log10.(abs.(S .- ref))
        defVar(ds, "log10_err_series", E, ("alpha", "theta", "L");
            attrib = field_attrib(E, ["long_name" => "|ψ_{3,0}^{(L)} − ψ_{3,0}|, Codex series",
                view_hint(vis = "heatmap", colormap = "Inferno", domain = [-15.0, 0.0], probe_dim = "L",
                    playback_dim = "L", value_transform = "exp10")...],
                write_coalescence!(ds); kind = :error))
        # Convergence line: worst error over the grid against L, with its running tail
        # envelope max_{L' ≥ L} (as in codex_series_convergence).
        M = [maximum(abs, S[:, :, l] .- ref) for l in axes(S, 3)]
        defVar(ds, "max_abs_err_series", M, ("L",);
            attrib = ["long_name" => "max over the grid |ψ_{3,0}^{(L)} − ψ_{3,0}|",
                view_hint(vis = "line", scale = "log", envelope = "cummax", envelope_direction = "reverse")...])
    end
    println("appended series error stack to ", path)
end

P = probe_terms()
line = probe_terms("line_terms.csv"; key = "theta")
fig_convergence(P, line)
fig_terms(P)
fig_inner_N(P)
αs, S, ref = grid_maps()
ρ = fig_error_maps(αs, S, ref)
fig_r1()
fig_assembly()
fig_comparison()
write_netcdf(αs, S, ref, ρ)
