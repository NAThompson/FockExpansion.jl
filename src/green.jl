# Liverts' azimuthally integrated S³ Green kernel for k=3.
# Full recurrence source S=10ψ31-2Vψ20+2Eψ10; ψ30=∫K S dΩ/(16π).
function green_kernel(α,θ,x,t)
    b=sin(α)*sin(x)*sin(θ)*sin(t)
    δ=2sin((α-x)/2)^2+2sin(α)*sin(x)*sin((θ-t)/2)^2
    den=δ+2b; s=sqrt(den); q=δ/den
    # Retain the complementary parameter near the logarithmic singularity.
    kk=q<eps(Float64) ? log(4)-log(q)/2 : ellipk(1-q)
    ee=q<eps(Float64) ? 1.0 : ellipe(1-q)
    a=δ+b
    (4kk/s-24s*ee+(16/3)*s*(4a*ee-δ*kk))/sqrt(2)
end
function green_source(x,t,Z,E,a21)
    v=sin(x)*cos(t); ξ=sqrt(cos(x)^2/(1+sin(x))+2sin(x)*sin(t/2)^2); σ=cos(x/2)+sin(x/2)
    p20=(1-2E)/12+Z^2*(1/3+sin(x)/2)+(a21-Z*(π+4)/(9π))*v
    if Z!=0; p20+=Z*(chi(x,cos(t))-σ*ξ/3); end
    p31=Z*(π-2)*ξ*(5ξ^2-6)/(36π)+Z^2*(π-2)*σ*v/(6π)
    (10p31-2*(1/ξ-2Z*σ/sin(x))*p20+2E*(ξ/2-Z*σ))*sin(x)^2*sin(t)/(16π)
end
"""Full Liverts Green-function inversion, Float64, using only lower Fock coefficients.
`method=:duffy` regularizes the kernel singularity; `:rectangular` (default) splits the integration at the target.
`rtol` controls nested adaptive quadrature, not a certified total error bound.
"""
function psi30_green(α::Real,θ::Real;Z,E,a21,rtol=1e-8,method=:rectangular)
    geometry(α,θ)
    α=min(Float64(α),π-Float64(α)); θ=Float64(θ)
    integrand(x,t)=green_source(x,t,Z,E,a21)*(green_kernel(α,θ,x,t)+green_kernel(α,θ,π-x,t))
    total=0.0
    for xe in (0.0,π/2), te in (0.0,Float64(π))
        dx=xe-α; dt=te-θ
        (dx==0 || dt==0) && continue
        if method==:rectangular
            val=quadgk(u->quadgk(v->integrand(α+dx*u,θ+dt*v)*abs(dx*dt),0.,1.;rtol,atol=rtol/100)[1],0.,1.;rtol,atol=rtol/100)[1]
        elseif method==:duffy
            function outer(s)
                u=s*s
                quadgk(v->(integrand(α+dx*u,θ+dt*u*v)+integrand(α+dx*u*v,θ+dt*u))*abs(dx*dt)*2s^3,0.,1.;rtol,atol=rtol/100)[1]
            end
            val=quadgk(outer,0.,1.;rtol,atol=rtol/100)[1]
        else
            throw(ArgumentError("method must be :duffy or :rectangular"))
        end
        total+=val
    end
    total
end
