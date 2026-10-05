# Harmonic (Y₄ₗ) part of ψ₄₁, fixed by solvability of the ψ₄₀ equation:
#   (Λ²-32)ψ₄₀ = 12ψ₄₁ + 2ψ₄₂ - 2Vψ₃₀ + 2Eψ₂₀   ⟹   ⟨Y,12ψ₄₁+2ψ₄₂-2Vψ₃₀+2Eψ₂₀⟩ = 0.
# The pure part of ψ₄₁ is ⊥ Y, so c_l = ⟨Y_l, Vψ₃₀ - Eψ₂₀ - ψ₄₂⟩/(6⟨Y_l,Y_l⟩).
# Fit c_l over the monomials of (Z, E, a₂₁) that can occur.
include("s3.jl"); using FockExpansion, Printf, LinearAlgebra
Y(l, a, t) = l==0 ? 4cos(a)^2-1 : l==1 ? 4sin(a)*cos(a)*cos(t) : sin(a)^2*(3cos(t)^2-1)/2
ξs(a,t)=sqrt(2sin(π/4-a/2)^2+2sin(a)*sin(t/2)^2)
p42(a,t,Z)=(π-2)*(5π-14)/(180π^2)*Z^2*(1-2sin(a)^2*sin(t)^2)
function cl(l, Z, E, a21; n=16, levels=9)
    f(a,t)=begin
        V=1/ξs(a,t)-Z*(1/sin(a/2)+1/cos(a/2))
        Y(l,a,t)*(V*psi30(a,t; Z, E, a21, rtol=1e-12)-E*psi20(a,t; Z, E, a21)-p42(a,t,Z))
    end
    sphint(f; n, levels)/(6sphint((a,t)->Y(l,a,t)^2; n, levels))
end
mono(Z,E,a)=[1, Z, Z^2, Z^3, Z^4, E, E*Z, E*Z^2, a, a*Z, a*Z^2]
names=["1","Z","Z²","Z³","Z⁴","E","EZ","EZ²","a","aZ","aZ²"]
samples=[(Z,E,a) for Z in (0.0, 0.5, 1.0, 2.0, 3.0), (E,a) in ((0.0,0.0),(-2.0,0.0),(0.0,0.5))][:]
@printf("c₁ (should vanish) at Z=2,E=-2.9,a=0.4: %.2e\n", cl(1, 2.0, -2.9, 0.4))
for l in (0, 2)
    A=reduce(vcat, [mono(s...)' for s in samples]); b=[cl(l, s...) for s in samples]
    coef=A\b
    chk=cl(l, 1.7, -2.2, 0.31); pred=dot(mono(1.7,-2.2,0.31), coef)
    @printf("l=%d fit residual %.1e, holdout %.2e\n", l, norm(A*coef-b), chk-pred)
    for (nm, c) in zip(names, coef); @printf("   %-4s % .15f\n", nm, c); end
    chk2=cl(l, 1.7, -2.2, 0.31; n=24, levels=11); @printf("   (n,levels)=(24,11) vs (16,9) at holdout: %.1e\n", chk2-chk)
end
