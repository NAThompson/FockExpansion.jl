# Compute and cache the angular fields used by psi30_figures.jl.
#
#     julia -t 4 --project=examples/sciviz examples/sciviz/psi30_compute.jl
#
# Writes
#   output/psi30.nc            513² interior grid: ψ30, its parts, Z sectors, E/a21
#                              sensitivities, gradient, Λ²ψ30 and the recurrence residual
#   output/psi30_agreement.nc  65² grid: ψ30 against the Green and Langner representations
# Existing files are kept; delete them to recompute.
include("common.jl")
using ForwardDiff, Printf

const N_MAIN = 513
const N_AGREE = 65

psi30_he(α, θ) = psi30(α, θ; HE...)

"""Coulomb potential ρV in hyperspherical angles."""
coulomb(α, θ, Z) = inv(sqrt(1 - sin(α) * cos(θ))) - Z * (sec(α / 2) + csc(α / 2))

"""Value, gradient and Hessian of ψ30 by ForwardDiff, and the pieces of the recurrence
(Λ² − 21)ψ30 = 10ψ31 − 2Vψ20 + 2Eψ10."""
function recurrence_point(α, θ; Z, E, a21)
    f(x) = psi30(x[1], x[2]; Z, E, a21)
    x = [α, θ]
    g = ForwardDiff.gradient(f, x)
    H = ForwardDiff.hessian(f, x)
    ψ = f(x)
    Λψ = -4 * (H[1, 1] + 2cot(α) * g[1] + (H[2, 2] + cot(θ) * g[2]) / sin(α)^2)
    rhs =
        10psi31(α, θ; Z) - 2coulomb(α, θ, Z) * psi20(α, θ; Z, E, a21) +
        2E * psi10(α, θ; Z)
    # Gradient norm in the metric of Λ² = −4Δ: |∇ψ|² = 4(ψ_α² + ψ_θ²/sin²α).
    gradnorm = 2sqrt(g[1]^2 + (g[2] / sin(α))^2)
    (; ψ, Λψ, rhs, gradnorm, lhs = Λψ - 21ψ)
end

function compute_main(path)
    αs = interior(N_MAIN)
    θs = interior(N_MAIN)
    @info "main grid" N_MAIN
    t0 = time()
    P = [psi30_parts(a, t; HE...) for a in αs, t in θs]  # fast; keep serial for clarity
    ψ = getfield.(P, :value)
    classical = getfield.(P, :classical)
    state = getfield.(P, :state)
    elliptic = getfield.(P, :elliptic)
    qerr = getfield.(P, :quadrature_error)
    cancellation =
        log10.((abs.(classical) .+ abs.(state) .+ abs.(elliptic)) ./ abs.(ψ))
    @info "parts" time() - t0

    # ψ30 is affine in E and a21; differentiate anyway (ForwardDiff) as a check.
    dE = grid(αs, θs) do a, t
        ForwardDiff.derivative(e -> psi30(a, t; Z = HE.Z, E = e, a21 = HE.a21), HE.E)
    end
    da21 = grid(αs, θs) do a, t
        ForwardDiff.derivative(c -> psi30(a, t; Z = HE.Z, E = HE.E, a21 = c), HE.a21)
    end
    dZ = grid(αs, θs) do a, t
        ForwardDiff.derivative(z -> psi30(a, t; Z = z, E = HE.E, a21 = HE.a21), HE.Z)
    end

    # At E = a21 = 0, ψ30 is a cubic polynomial in Z: solve for its coefficients.
    Zs = (0.0, 1.0, 2.0, 3.0)
    vals = [grid((a, t) -> psi30(a, t; Z = z, E = 0.0, a21 = 0.0), αs, θs) for z in Zs]
    V = [z^n for z in Zs, n = 0:3]
    Vinv = inv(V)
    sectors = [sum(Vinv[n+1, k] * vals[k] for k = 1:4) for n = 0:3]
    check = grid((a, t) -> psi30(a, t; Z = 4.0, E = 0.0, a21 = 0.0), αs, θs)
    pred = sum(4.0^n * sectors[n+1] for n = 0:3)
    @printf("cubic-in-Z check at Z=4: max |error| = %.2e (max |ψ| = %.2e)\n",
        maximum(abs, pred - check), maximum(abs, check))
    @info "sectors" time() - t0

    R = Matrix{Any}(undef, N_MAIN, N_MAIN)
    Threads.@threads for j in eachindex(θs)
        for i in eachindex(αs)
            R[i, j] = recurrence_point(αs[i], θs[j]; HE...)
        end
    end
    Λψ = [r.Λψ for r in R]
    rhs = [r.rhs for r in R]
    gradnorm = [r.gradnorm for r in R]
    resid = [r.lhs - r.rhs for r in R]
    relresid = abs.(resid) ./ max.(abs.(rhs), abs.(Λψ))
    @info "derivatives" time() - t0

    symlog(cm) = view_hint(vis = "heatmap", scale = "symlog", colormap = cm)
    lin(cm) = view_hint(vis = "heatmap", scale = "linear", colormap = cm)
    logv(cm) = view_hint(vis = "heatmap", scale = "log", colormap = cm)
    longname(s) = ["long_name" => s]
    fields = [
        "psi30" => ψ,
        "psi30_classical" => classical,
        "psi30_state" => state,
        "psi30_elliptic" => elliptic,
        "psi30_quadrature_error" => qerr,
        "log10_cancellation" => cancellation,
        "dpsi30_dE" => dE,
        "dpsi30_da21" => da21,
        "dpsi30_dZ" => dZ,
        "sector_Z0" => sectors[1],
        "sector_Z1" => sectors[2],
        "sector_Z2" => sectors[3],
        "sector_Z3" => sectors[4],
        "grad_norm" => gradnorm,
        "lambda2_psi30" => Λψ,
        "recurrence_rhs" => rhs,
        "recurrence_residual" => resid,
        "recurrence_residual_abs" => abs.(resid),
        "recurrence_residual_rel" => relresid,
    ]
    attribs = Dict{String,Any}(
        "psi30" => [longname("ψ_{3,0}(α, θ), helium ground state"); symlog("auto")],
        "psi30_classical" => [longname("dilogarithmic (classical) part"); symlog("auto")],
        "psi30_state" => [longname("E, a_{21} dependent part"); symlog("auto")],
        "psi30_elliptic" => [longname("four elliptic quadratures"); symlog("auto")],
        "psi30_quadrature_error" => [longname("quadrature error estimate"); logv("Viridis")],
        "log10_cancellation" =>
            [
                longname("sum |parts| / |ψ_{3,0}|")
                lin("Inferno")
                view_hint(value_transform = "exp10")
            ],
        "dpsi30_dE" => [longname("∂ψ_{3,0}/∂E"); symlog("auto")],
        "dpsi30_da21" => [longname("∂ψ_{3,0}/∂a_{21}"); symlog("auto")],
        "dpsi30_dZ" => [longname("∂ψ_{3,0}/∂Z at fixed E, a_{21}"); symlog("auto")],
        "sector_Z0" => [longname("Z^{0} coefficient at E = a_{21} = 0"); symlog("auto")],
        "sector_Z1" => [longname("Z^{1} coefficient at E = a_{21} = 0"); symlog("auto")],
        "sector_Z2" => [longname("Z^{2} coefficient at E = a_{21} = 0"); symlog("auto")],
        "sector_Z3" => [longname("Z^{3} coefficient at E = a_{21} = 0"); symlog("auto")],
        "grad_norm" => [longname("|∇ψ_{3,0}| on the hypersphere"); logv("Viridis")],
        "lambda2_psi30" => [longname("Λ^{2}ψ_{3,0} (ForwardDiff)"); symlog("auto")],
        "recurrence_rhs" => [longname("10ψ_{3,1} − 2Vψ_{2,0} + 2Eψ_{1,0}"); symlog("auto")],
        "recurrence_residual" => [longname("(Λ^{2}−21)ψ_{3,0} − RHS"); symlog("auto")],
        "recurrence_residual_abs" => [longname("|(Λ^{2}−21)ψ_{3,0} − RHS|"); logv("Viridis")],
        "recurrence_residual_rel" => [longname("|(Λ^{2}−21)ψ_{3,0} − RHS| / max(|RHS|, |Λ^{2}ψ_{3,0}|)"); logv("Viridis")],
    )
    write_angular(path, αs, θs, fields; attribs,
        global_attribs = Pair{String,Any}[
            "title" => "Third-order Fock coefficient psi30, helium ground state",
            "Z" => HE.Z, "E" => HE.E, "a21" => HE.a21,
        ])
end

function compute_agreement(path)
    αs = interior(N_AGREE; ε = 0.02)
    θs = interior(N_AGREE; ε = 0.02)
    table = LangnerTable(; Z = HE.Z)
    t0 = time()
    ψ = grid(psi30_he, αs, θs)
    lang = grid((a, t) -> psi30_langner(a, t, table; E = HE.E, a21 = HE.a21), αs, θs)
    @info "langner" time() - t0
    green = grid((a, t) -> psi30_green(a, t; HE..., rtol = 1e-10), αs, θs)
    @info "green" time() - t0
    # Error of the Float64 elliptic formula at its default rtol, against BigFloat.
    exact = setprecision(BigFloat, 256) do
        bigpar = (Z = big(HE.Z), E = big(HE.E), a21 = big(HE.a21))
        grid((a, t) -> Float64(psi30(big(a), big(t); bigpar...)), αs, θs)
    end
    @info "bigfloat" time() - t0
    lv(cm) = view_hint(vis = "heatmap", scale = "log", colormap = cm)
    heatmap = view_hint(vis = "heatmap")
    write_angular(path, αs, θs, [
            "psi30" => ψ,
            "psi30_green" => green,
            "psi30_langner" => lang,
            "abs_diff_green" => abs.(ψ - green),
            "abs_diff_langner" => abs.(ψ - lang),
            "elliptic_error" => abs.(ψ - exact),
        ];
        attribs = Dict{String,Any}(
            "psi30" => ["long_name" => "ψ_{3,0}, elliptic formula"; heatmap],
            "psi30_green" => ["long_name" => "ψ_{3,0}, Green integral (rtol 1e-10)"; heatmap],
            "psi30_langner" => ["long_name" => "ψ_{3,0}, Langner coefficient table"; heatmap],
            "abs_diff_green" => ["long_name" => "|ψ_{3,0} − ψ_{3,0}^{Green}|"; lv("Viridis")],
            "abs_diff_langner" =>
                ["long_name" => "|ψ_{3,0} − ψ_{3,0}^{Langner}|"; lv("Viridis")],
            "elliptic_error" => [
                "long_name" => "|ψ_{3,0} − ψ_{3,0}^{BigFloat}|, elliptic at default rtol";
                lv("Viridis")
            ],
        ),
        # The root group opens as panels: the three representations above, their
        # differences and the elliptic error below (2 × 3), with linked axes and cursor.
        global_attribs = Pair{String,Any}[
            "title" => "ψ_{3,0}: elliptic formula vs Green and Langner representations";
            view_hint(
                panels = "psi30 psi30_green psi30_langner abs_diff_green abs_diff_langner elliptic_error",
                layout = [2, 3],
                link_axes = true,
                link_cursor = true,
                symmetric = true,
            )
        ],
    )
end

main_path = joinpath(OUT, "psi30.nc")
agree_path = joinpath(OUT, "psi30_agreement.nc")
isfile(main_path) || compute_main(main_path)
isfile(agree_path) || compute_agreement(agree_path)
