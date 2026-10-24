# Why origin (Fock) data cannot reach large R: exponential type of the coefficients
# and amplification of 1e-12 perturbations of (a, E). Generic-typed copy of cuspcoeff.
include("exterior_lib.jl")
using Printf, Random, SpecialFunctions
const OUT = @__DIR__

function cuspcoeff_generic(b, E::S, N) where {S}
    m = length(b.K); pmax = N ÷ 2 + 1
    U = S.(b.U); Dm = S.(b.D); T = S.(b.W) + 2E * I; lam = S.(b.lam)
    c = zeros(S, N + 1, pmax + 2, m, m)
    for n = 0:N
        src = zeros(S, pmax + 2, m, m)
        for k = 0:n÷2
            if n >= 1
                src[k+1, :, :] -= (2(n - 1) * U + Dm) * c[n, k+1, :, :] + 2(k + 1) * U * c[n, k+2, :, :]
            end
            n >= 2 && (src[k+1, :, :] -= T * c[n-1, k+1, :, :])
        end
        for i = 1:m
            den = n * (n + 4) - lam[i]
            if den == 0
                for k = n÷2:-1:0
                    c[n+1, k+2, i, :] = (src[k+1, i, :] - (k + 1) * (k + 2) * c[n+1, k+3, i, :]) / ((2n + 4) * (k + 1))
                end
                c[n+1, 1, i, i] = 1
            else
                for k = n÷2:-1:0
                    c[n+1, k+1, i, :] = (src[k+1, i, :] - (2n + 4) * (k + 1) * c[n+1, k+2, i, :] - (k + 1) * (k + 2) * c[n+1, k+3, i, :]) / den
                end
            end
        end
    end
    c
end

fock_value(c, a, R) = sum(R^n * log(R)^k * (c[n+1, k+1, :, :] * a) for n = 0:size(c, 1)-1 for k = 0:n÷2)

ds = isempty(ARGS) ? [2, 4] : parse.(Int, ARGS)
setprecision(BigFloat, 256)
for d in ds
    b = cuspmodel(d); z = cuspsolve(b); E = z.E
    slow, fast = pencil_roots(b, E)
    N = 170
    t = @elapsed c = cuspcoeff_generic(b, big(E), N)
    a = big.(z.a)
    @printf("d=%d E=%.12f N=%d coeff %.1fs  max|lambda|=%.4f  max slow=%.4f max fast=%.4f\n", d, E, N, t,
        maximum(abs.(vcat(slow, fast))), maximum(real(slow)), maximum(real(fast)))
    open(joinpath(OUT, "type_d$d.csv"), "w") do io
        println(io, "n,sig_a,sig_generic")
        for n = 10:10:N
            na = maximum(norm(c[n+1, k+1, :, :] * a) for k = 0:n÷2)
            ng = maximum(opnorm(Float64.(c[n+1, k+1, :, :])) for k = 0:n÷2)
            sa = Float64(exp((logfactorial(n) + log(na)) / n))
            sg = Float64(exp((logfactorial(n) + log(big(ng))) / n))
            println(io, "$n,$sa,$sg")
            n % 40 == 0 && @printf("  n=%3d (n! max_k|c_nk a|)^(1/n)=%.3f  (n! max_k|c_nk|)^(1/n)=%.3f\n", n, sa, sg)
        end
    end
    Random.seed!(1)
    da = big.(randn(length(a))); da *= big(1e-12) * norm(a) / norm(da)
    cE = cuspcoeff_generic(b, big(E) * (1 + big(1e-12)), N)
    open(joinpath(OUT, "sensitivity_d$d.csv"), "w") do io
        println(io, "R,amp_a,amp_E,phi_norm")
        for R in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0)
            Rb = big(R)
            f0 = fock_value(c, a, Rb)
            fa = fock_value(c, da, Rb)
            fE = fock_value(cE, a, Rb) - f0
            amp_a = Float64(norm(fa) / norm(f0) / 1e-12)
            amp_E = Float64(norm(fE) / norm(f0) / 1e-12)
            println(io, "$R,$amp_a,$amp_E,$(Float64(norm(f0)))")
            @printf("  R=%.1f |Phi|=%.3e  amplification(a)=%.3e  amplification(E)=%.3e\n", R, Float64(norm(f0)), amp_a, amp_E)
        end
    end
    flush(stdout)
end
