# λ_min, λ_2 and a21 vs E over E0 ± 1e-5 (L = 8, B10P6, R* = 1, d = 32). No cache writes.
include(joinpath(@__DIR__, "common_compute.jl"))

const RS, L, D, NAME = 1.0, 8, 32, "B10P6"
dEs = collect(range(-1e-5, 1e-5; length = 11))
zsets, nchan = basis_spec(NAME)
model0 = channel_model(L)
nL = length(model0.basis.K)
lam1 = Float64[]
lam2 = Float64[]
a21 = Float64[]
secs = Float64[]
for dE in dEs
    E = E_HE + dE
    t = @elapsed begin
        ext = exterior_form_dd(model0, E, RS, nL; zsets, nchan, lchan = L)
        p = solve_point(D, L, RS, E, ext)
    end
    push!(lam1, p.lambda)
    push!(lam2, p.lambda2)
    push!(a21, p.a21)
    push!(secs, t)
    @printf("dE=%+.1e: lambda_min=%+.6e lambda_2=%.6e a21=%.13f (%.1fs)\n", dE, p.lambda,
            p.lambda2, p.a21, t)
    flush(stdout)
end
serialize(joinpath(DATA, "energy.jls"), Dict{String,Any}(
    "Rs" => RS, "L" => L, "d" => D, "basis" => NAME, "E0" => E_HE, "dE" => dEs,
    "lambda_min" => lam1, "lambda_2" => lam2, "a21" => a21, "seconds" => secs,
))
println("wrote ", joinpath(DATA, "energy.jls"))
