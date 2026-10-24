include("steklov.jl")
include("exterior_dd.jl")
using Serialization
function main()
BLAS.set_num_threads(1)
println("POINT GRADIENTS")
b = harmonic_basis(14)
T = Double64
x = T.([1.2, 0.7, 0.9]); Rs = T(1)
nL = length(b.K); P = 3
V = zeros(T,1,nL*(P+1),4)
eval3d!(V,1,b,nL,nL,P,x...,Rs,2,1)
err = zero(T)
for i in (1,2,3,4,nL), p in (-1,0,2)
    z = p == -1 ? lift_vg(b,i,x...,Rs,:gauge,T(2),T(1),:cusp,T(1.5)) : channel_vg(b,i,p,x...,Rs,T(2),T(1))
    c = p == -1 ? i : nL+(i-1)*P+p+1
    err = max(err,maximum(abs.(collect(z)-vec(V[1,c,:]))))
end
println("max analytic vs AD = ",err)
println("ANALYTIC R12 BLOCK")
terms = [Term(0,0,0,2.,1.3),Term(1,2,1,2.,1.3),Term(0,0,4,4.,4.)]
kw = (;nR=4,nphi=5,edges=[0.5,2.],npanel=4)
r2 = exterior_rule_2d(1.,T;kw...)
r3 = exterior_rule_3d(1.,T;kw...,nu=12)
A = correction_block(terms,r2,1.,E_HE,2.,1.,T)
fs = [(a,c,u)->term_vg_sym(t,a,c,u,T(1)) for t in terms]
C,_ = lift_blocks(fs,Term[],r3,1.,E_HE,2.,1.,T)
println("relative block difference = ",maximum(abs.(A-C))/maximum(abs.(A)))
println("HIGH L ANGULAR GRAM: nphi=24,nu=20 versus 40,40")
for L in (4,12,14)
    bb=harmonic_basis(L)
    grams=[]
    for (np,nu) in ((24,20),(40,40))
        r=sphere_rule(np,nu);Y=harmonic_values(bb,sphere_q(r),sphere_s(r))
        push!(grams,Y'*((r.w .* exp.(2 .* cusp_U.(r.c,r.s,r.u;Z=2.,g=1.))).*Y))
    end
    println("L=",L," relative Gram discrepancy=",norm(grams[1]-grams[2],Inf)/norm(grams[2],Inf))
end
println("2S SPLIT")
for d in (16,24,28)
    model=channel_model(d;g=0.)
    ser=fock_series(model,-1.,2d+28,[1.])
    r=sphere_rule(2d+40,2d+40);Y=harmonic_values(model.basis,sphere_q(r),sphere_s(r))
    y=Y'*(r.w .* ((1 .- r.c).*exp.(r.c).*(1 .- r.s).*exp.(r.s)))
    i=findfirst(==((2,0)),model.basis.labels)
    for L in (8,12,d)
        t=copy(y);t[harmonic_degree(model.basis).>L].=0
        a=ser.F[1]\t;a./=a[1]
        @printf("d=%d L=%d seed_error=%+.6e\n",d,L,a[i]+1/32)
    end
end
println("N CHECK AND SOLVE PRECISION (d=16,L=12)")
ext=deserialize("cache/ext_L12_R1.0_E-2.9037243770341195_B10P6.jls")
m=channel_model(16)
for N in (60,80)
    s=interior_setup(m,E_HE,1.,12;N)
    z=match(s,Float64.(ext.S));pc=physical_constants(s,z.a)
    println("N=",N," a21=",pc.a21," condF=",cond(s.F))
    aa=BigFloat.(s.F)\BigFloat.(vcat(z.t,zeros(length(z.a)-s.nL)));aa./=aa[1]
    println("BigFloat solve same F/t low seed change = ",Float64(aa[2])-z.a[2])
end
end
main()
