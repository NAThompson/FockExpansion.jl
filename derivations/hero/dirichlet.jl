# Energy from the Fock series alone: the Dirichlet condition Ψ(ρ_out, Ω) = 0 on a hypersphere.
# Ψ = Φ₀ + Σⱼ aⱼΦⱼ is linear in the free constants at fixed E, so a nonzero solution exists
# where the matrix [Φⱼ(ρ_out, Ωᵢ)] (rows weighted by the S³ measure) becomes singular.
# E_D(ρ_out) is the ground-state energy of helium confined to the hypersphere and decreases
# to the free-atom energy as ρ_out grows.
# usage: julia -t 4 --project=. dirichlet.jl T n K kfree "ρ1,ρ2,..." Emin Emax nE   (scan)
#        julia -t 4 --project=. dirichlet.jl T n K kfree ρ Elo Ehi root           (root)
using DoubleFloats, Printf, LinearAlgebra
include("../k4/fock_spectral.jl")
T=ARGS[1]=="Double64" ? Double64 : Float64
n=parse(Int, ARGS[2]); K=parse(Int, ARGS[3]); kfree=parse(Int, ARGS[4])
ρs=T.(parse.(Float64, split(ARGS[5], ",")))
Emin, Emax=parse(Float64, ARGS[6]), parse(Float64, ARGS[7]); nE=ARGS[8]=="root" ? 0 : parse(Int, ARGS[8])
g=SGrid(T, n)
# angular collocation: tensor Gauss-Chebyshev nodes in (α, θ), weights √(sin²α sin θ)
m=16
angs=[(T(π)*(i-T(1)/2)/m, T(π)*(j-T(1)/2)/m) for i in 1:m for j in 1:m]
wts=[sin(a)*sqrt(sin(t)) for (a, t) in angs]
# barycentric weights for each angle: value = ax'·F·ay
function baryvecs(a, θ)
    P=T(π); β=acos(sin(a)*cos(θ))
    tx, ty=2(a-β)/P, 2(a+β-P)/P
    x=[cos(P*(j+T(1)/2)/n) for j in 0:n-1]
    bw=[(-1)^j*sin(P*(j+T(1)/2)/n) for j in 0:n-1]
    ax=bw./(tx.-x); ay=bw./(ty.-x)
    ax/sum(ax), ay/sum(ay)
end
const BV=[baryvecs(a, t) for (a, t) in angs]
frees=[(k, l) for k in 2:2:kfree for l in 0:k÷2 if iseven(k÷2-l)]
# Φ at (ρ_out, Ωᵢ) for one chain, summed on the fly
function series(E, free)
    S=zeros(T, length(angs), length(ρs))
    S[:, :].=1   # ψ₀₀ = 1 (zero for responses, subtracted below)
    sfock(g; Z=2, E, kmax=K, free, keep=2, onk=(k, ψ, t)->begin
        for p in 0:k÷2
            v=[ax'*ψ[(k, p)]*ay for (ax, ay) in BV]
            for (r, ρ) in enumerate(ρs); S[:, r].+=ρ^k*log(ρ)^p.*v; end
        end
    end)
    S
end
function smin(E)
    E=T(E)
    cols=Vector{Matrix{T}}(undef, length(frees)+1)
    Threads.@threads for j in 1:length(frees)+1
        cols[j]=j==1 ? series(E, Dict{Tuple{Int,Int},Any}()) : series(E, Dict{Tuple{Int,Int},Any}(frees[j-1]=>1))
    end
    for j in 2:length(cols); cols[j]=cols[j]-cols[1]; end
    map(eachindex(ρs)) do r
        A=reduce(hcat, [wts.*c[:, r] for c in cols])
        A=A./sqrt.(sum(abs2, A; dims=1))   # unit columns: the responses scale like ρ^k
        s=svdvals(A); Float64(s[end]/s[1])
    end
end
@printf("T=%s n=%d K=%d kfree=%d (%d free constants), %d angles\n", T, n, K, kfree, length(frees), length(angs))
if ARGS[8]=="root"
    # σmin(E) ≈ c|E-E_D| near the root: intersect the two sides of the V and shrink.
function vroot(f, a, b)
    fa, fb=f(a), f(b)
    m=(a+b)/2; fm=f(m)
    for it in 1:30
        # three points bracketing the minimum; slopes from the outer pairs
        if fa<fm || fb<fm
            error("minimum not bracketed: f($a)=$fa, f($m)=$fm, f($b)=$fb")
        end
        sl, sr=(fm-fa)/(m-a), (fb-fm)/(b-m)   # sl < 0 < sr
        # line through (a,fa) slope sl and line through (b,fb) slope sr; take the side with
        # the smaller value at m as the reference V
        E=(fb-fa+sl*a-sr*b)/(sl-sr)
        E=clamp(E, a+(m-a)/10, b-(b-m)/10)
        # a step that lands on m gives no information: probe the larger side instead
        abs(E-m)<(b-a)/1000 && (E=b-m>m-a ? m+(b-m)/3 : m-(m-a)/3)
        fE=f(E)
        @printf("  it %2d  E = %.12f  σmin/σmax = %.2e  bracket %.1e\n", it, E, fE, b-a); flush(stdout)
        if E<m
            fE<fm ? ((b, fb, m, fm)=(m, fm, E, fE)) : ((a, fa)=(E, fE))
        else
            fE<fm ? ((a, fa, m, fm)=(m, fm, E, fE)) : ((b, fb)=(E, fE))
        end
        b-a<1e-10 && break
    end
    m
end
    @printf("E_D(ρ=%.2f) = %.12f\n", Float64(ρs[1]), vroot(E->smin(E)[1], Emin, Emax))
else
    for E in range(Emin, Emax; length=nE)
        t=@elapsed s=smin(E)
        @printf("E=%.6f  σmin/σmax:", E)
        for (ρ, v) in zip(ρs, s); @printf("  ρ=%.1f %.2e", Float64(ρ), v); end
        @printf("   (%.0f s)\n", t); flush(stdout)
    end
end
