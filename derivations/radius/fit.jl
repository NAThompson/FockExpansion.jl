# Fit the free constants of the Fock series to an accurate Hylleraas wavefunction and see
# how far out in ρ the truncated series agrees with it.
# Ψ = λ(Φ₀₀ + Σⱼ aⱼΦⱼ) is linear in (λ, λaⱼ) at fixed E: Φⱼ is the chain's response to the
# free constant aⱼ (coefficient of the exchange-symmetric harmonic Y_kl, k/2-l even).
include("../k4/fock_spectral.jl")
include("hylleraas.jl")
using Serialization
setprecision(BigFloat, 320)
K=parse(Int, get(ARGS, 1, "24")); n=parse(Int, get(ARGS, 2, "48"))
const Eex=-2.903724377034119598311159245194
file=joinpath(tempdir(), "hylleraas_12.jls")   # cache
h=isfile(file) ? deserialize(file) : (h=hylleraas(12, (1.6, 6.0)); serialize(file, h); h)
@printf("Hylleraas: N=%d, E-Eexact = %.1e\n", length(h.idx), Float64(h.E-big"-2.903724377034119598311159245194"))
frees=[(k, l) for k in 2:2:K for l in 0:k÷2 if iseven(k÷2-l)]
g=SGrid(Float64, n)
# angular sample points
angs=[(π*(i-0.5)/10, π*(j-0.5)/10) for i in 1:10 for j in 1:10]
# Φ values at the angular points: vals[j][(k,p)] = vector over angles
function series_vals(ψ)
    Dict(kp=>[interp(g, F, a, t) for (a, t) in angs] for (kp, F) in ψ)
end
base=series_vals(sfock(g; Z=2, E=Eex, kmax=K))
Φ=[base]
for f in frees
    v=series_vals(sfock(g; Z=2, E=Eex, kmax=K, free=Dict(f=>1.0)))
    push!(Φ, Dict(kp=>v[kp]-base[kp] for kp in keys(v)))
end
# series of each Φ at (ρ, angle index)
Sρ(F, ρ, i)=sum(ρ^k*log(ρ)^p*F[(k, p)][i] for k in 0:K for p in 0:k÷2)
function href(ρ, i)
    a, t=big.(angs[i])
    Float64(hyl_eval(h, big(ρ)*cos(a/2), big(ρ)*sin(a/2), cos(t)))
end
ρfit=parse(Float64, get(ARGS, 3, "1.0"))
ρs=collect(range(0.05, ρfit; length=20))
A=reduce(hcat, [[Sρ(F, ρ, i) for ρ in ρs for i in eachindex(angs)] for F in Φ])
b=[href(ρ, i) for ρ in ρs for i in eachindex(angs)]
c=A\b
@printf("K=%d, %d free constants, fit on ρ ≤ %.2f; λ = %.10f, a21 = %.10f, a40 = %.8f, a42 = %.8f\n",
        K, length(frees), ρfit, c[1], c[2]/c[1], c[3]/c[1], c[4]/c[1])
for ρ in (0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0)
    e=maximum(abs(sum(c[j]*Sρ(Φ[j], ρ, i) for j in eachindex(Φ))-href(ρ, i)) for i in eachindex(angs))
    s=maximum(abs(href(ρ, i)) for i in eachindex(angs))
    @printf("ρ=%.2f  max |series - Hylleraas| / max|ψ| = %.1e\n", ρ, e/s)
end
