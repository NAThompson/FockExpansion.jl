# Outer formal solution of helium (L = 0) in the He⁺(1s) + e channel.
# r₁ = inner electron, R = r₂ = outer electron, ψ = Σ_l f_l(r₁, R) P_l(cos θ₁₂), and
#   R f_l = e^{-κR} Σ_k R^{σ-k} g_kl(r₁),   g_kl(r) = e^{-Zr} · polynomial(r).
# With 1/r₁₂ = Σ_λ r₁^λ/R^{λ+1} P_λ and A_{ll'λ} = (2l+1)(l λ l'; 0 0 0)², the order R^{σ-k}:
#   L_l g_kl = κ(k-1) g_{k-1,l} + ½[(σ-k+2)(σ-k+1) - l(l+1)] g_{k-2,l}
#              - Σ_{λ≥1} Σ_l' A_{ll'λ} r^λ g_{k-λ-1,l'},
#   L_l = -½(d²/dr² + (2/r)d/dr) + l(l+1)/(2r²) - Z/r + Z²/2,
# after k = 0 fixes κ = √(2(-Z²/2 - E)) (g₀₀ = e^{-Zr}) and k = 1 fixes σ = (Z-1)/κ.
# L_l r^j e^{-Zr} = [½(l(l+1) - j(j+1)) r^{j-2} + Z j r^{j-1}] e^{-Zr}: triangular, so each
# step is a downward recursion. For l = 0, L₀ annihilates e^{-Zr}: its coefficient in
# g_{k-1,0} is fixed by solvability at order k (the r¹ coefficient of g_k0 must vanish),
# exactly like the harmonic parts in Fock's recurrence.
using Printf
setprecision(BigFloat, 512)
const Z=big(2)
function threej2(l, λ, lp)            # (l λ l'; 0 0 0)²
    J=l+λ+lp
    (isodd(J) || λ>l+lp || l>λ+lp || lp>l+λ) && return big(0)//1
    g=J÷2
    f(n)=factorial(big(n))
    (f(J-2l)*f(J-2λ)*f(J-2lp)//f(J+1))*(f(g)//(f(g-l)*f(g-λ)*f(g-lp)))^2
end
A(l, lp, λ)=(2l+1)*threej2(l, λ, lp)
# solve L_l g = b (both as coefficient vectors of r^j e^{-Zr}, index j+1); a₀ = 0
function solveL(l, b)
    M=length(b)-1
    a=zeros(BigFloat, M+3)
    for m in M:-1:0
        a[m+2]=(b[m+1]-(l*(l+1)-(m+2)*(m+3))/2*a[m+3])/(Z*(m+1))
    end
    a
end
padd(p, q)=(n=max(length(p), length(q)); [get(p, i, big(0))+get(q, i, big(0)) for i in 1:n])
shift(p, λ)=vcat(zeros(BigFloat, λ), p)
function outer(E, K)
    κ=sqrt(2(-Z^2/2-E)); σ=(Z-1)/κ
    g=Dict{Tuple{Int,Int},Vector{BigFloat}}((0, 0)=>[big(1)])
    get0(k, l)=get(g, (k, l), BigFloat[])
    for k in 1:K
        for l in 0:k
            b=κ*(k-1)*get0(k-1, l)
            b=padd(b, ((σ-k+2)*(σ-k+1)-l*(l+1))/2*get0(k-2, l))
            for λ in 1:k-1, lp in max(0, l-λ):(l+λ)
                p=get0(k-λ-1, lp); isempty(p) && continue
                c=A(l, lp, λ); iszero(c) && continue
                b=padd(b, -c*shift(p, λ))
            end
            all(iszero, b) && l>0 && continue
            if l==0 && k>=2
                # fix the e^{-Zr} coefficient t of g_{k-1,0}: it adds κ(k-1)t to b₀
                a=solveL(0, b)
                t=-a[2]*Z/(κ*(k-1))
                g[(k-1, 0)][1]+=t
                b[1]+=κ*(k-1)*t
            end
            a=solveL(l, b)
            l==0 && @assert abs(a[2])<big(10)^-100
            while length(a)>1 && iszero(a[end]); pop!(a); end
            g[(k, l)]=a
        end
    end
    (; g, κ, σ)
end
evalg(p, r)=sum(c*r^(j-1) for (j, c) in enumerate(p); init=zero(r))*exp(-Z*r)
