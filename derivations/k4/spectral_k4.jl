# ψ₄₂, ψ₄₁, ψ₄₀ from the recurrence by the separated spectral solve (see sylvester30.jl),
# compared with the package's Feynman (ψ₄₁) and Green's-function (ψ₄₀) evaluators.
include("sylvester30.jl")
# Fejér weights for Chebyshev points of the first kind on [-1,1]
function fejer(n)
    θ=[π*(j+0.5)/n for j in 0:n-1]
    [2/n*(1-2sum(cos(2k*θj)/(4k^2-1) for k in 1:n÷2)) for θj in θ]
end
n=24
X, Y, MX, MY=separated_ops(n, 32.0)
wq=fejer(n).*(π/2)                                   # dX = (π/2)dx
W=[(π/2)*(cos(X[i])-cos(Y[j]))*wq[i]*wq[j] for i in 1:n, j in 1:n]   # S³ measure
nodes=[αθ(X[i], Y[j]) for i in 1:n, j in 1:n]
x1=[cos((X[i]+Y[j])/2) for i in 1:n, j in 1:n]; x2=[cos((Y[j]-X[i])/2) for i in 1:n, j in 1:n]
Hs=[4x1.^2 .- 1, x1.*x2, (3x2.^2 .- (1 .- x1.^2))./2]               # Y₄₀, Y₄₁, Y₄₂
ip(f, g)=sum(W.*f.*g)
Gram=[ip(a, b) for a in Hs, b in Hs]
purify(F)=(c=Gram\[ip(h, F) for h in Hs]; F-sum(c[l]*Hs[l] for l in 1:3))
# pure inverse of (Λ²-32): eigen-decomposition, resonant pairs dropped
eX=eigen(Matrix(MX)); eY=eigen(Matrix(MY))
VX, λX, VY, λY=real(eX.vectors), real(eX.values), real(eY.vectors), real(eY.values)
iVX, iVY=inv(VX), inv(VY)
function Lplus(H)                                   # H = (cos X - cos Y)·h on the grid
    Rt=iVX*(-H./16)*iVY'
    F=similar(Rt); nres=0
    for i in 1:n, j in 1:n
        d=λX[i]-λY[j]
        if abs(d)<1e-6
            F[i, j]=0; nres+=1
        else
            F[i, j]=Rt[i, j]/d
        end
    end
    purify(VX*F*VY'), nres
end
Z, E, a21 = 2.0, -2.9, 0.4
w=[2sin((X[i]+Y[j])/2)*sin((Y[j]-X[i])/2) for i in 1:n, j in 1:n]   # cos X - cos Y
wV=[(α=(X[i]+Y[j])/2; β=(Y[j]-X[i])/2; 2sqrt(2)*sin(α)*cos(β/2)-4Z*sin(β)*(cos(α/2)+sin(α/2))) for i in 1:n, j in 1:n]
grid(f)=[f(nodes[i, j]...) for i in 1:n, j in 1:n]
ψ20=grid((a,t)->psi20(a, t; Z, E, a21)); ψ21=grid((a,t)->psi21(a, t; Z))
ψ30=grid((a,t)->psi30(a, t; Z, E, a21, rtol=1e-14)); ψ31=grid((a,t)->psi31(a, t; Z))
ψ42=grid((a,t)->psi42(a, t; Z))
# ψ₄₁: (Λ²-32)ψ₄₁ = 24ψ₄₂ - 2Vψ₃₁ + 2Eψ₂₁  (pure part), harmonic part from ψ₄₀ solvability
ψ41p, nres=Lplus(w.*(24ψ42+2E*ψ21)-2wV.*ψ31)
@printf("resonant pairs dropped: %d\n", nres)
src40=w.*(2ψ42+2E*ψ20)-2wV.*ψ30                    # (cos X - cos Y)(2ψ₄₂ - 2Vψ₃₀ + 2Eψ₂₀)
c=Gram\[-(ip(h, 12w.*ψ41p+src40)/12) for h in Hs]./[ip(h, w) for h in Hs]  # placeholder, replaced below
# solvability: ⟨Y, 12ψ₄₁ + 2ψ₄₂ - 2Vψ₃₀ + 2Eψ₂₀⟩ = 0 with the plain S³ measure (W already has w)
ipw(f, g)=sum(W.*f.*g)
S=[ipw(h, 2ψ42+2E*ψ20)-2sum(W./w.*h.*wV.*ψ30) for h in Hs]  # ⟨Y, 2ψ₄₂ + 2Eψ₂₀ - 2Vψ₃₀⟩
c41=-(Gram\S)./12 .- Gram\[ipw(h, ψ41p) for h in Hs]
ψ41=ψ41p+sum(c41[l]*Hs[l] for l in 1:3)
ψ40, _=Lplus(w.*(12ψ41+2ψ42+2E*ψ20)-2wV.*ψ30)
# compare with the package at a few grid nodes
for (i, j) in ((5, 9), (12, 12), (18, 4), (20, 20))
    a, t=nodes[i, j]
    p41=psi41(a, t; Z, E, a21); p40=psi40(a, t; Z, E, a21)
    @printf("(α,θ)=(%.3f,%.3f): ψ₄₁ %+.2e   ψ₄₀ %+.2e  (spectral - package)\n", a, t, ψ41[i,j]-p41, ψ40[i,j]-p40)
end
