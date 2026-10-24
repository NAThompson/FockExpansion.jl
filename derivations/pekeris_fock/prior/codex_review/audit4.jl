include("../exterior_factorial/exterior_lib.jl")
include_string(Main, join(split(first(split(read("../exterior_factorial/twoseries.jl", String), "base = Dict(")), "\n")[4:end], "\n"))
println("ambient precision=", precision(BigFloat))
for d in (2,4,6)
 b=cuspmodel(d); E=big(cuspsolve(b).E)
 setprecision(BigFloat,1024) do
  Ys=riccati_coeffs_bigE(b,E,3;prec=1024); Y=Ys[1]; U=big.(b.U); T=big.(b.W)+2E*I
  Q=Y*Y+2U*Y+T
  println("d=",d," solvent residual=",Float64(norm(Q))," precision=",precision(Y[1]))
  # Float64 spectra suffice to identify the separated branches, not certify high digits.
  A=[zeros(length(b.K),length(b.K)) I; -Float64.(T) -2b.U]
  println(" selected eigs=",sort(real(eigvals(Float64.(Y)))))
  println(" pencil eigs=",sort(real(eigvals(A))))
 end
end
b=cuspmodel(2); p=(Rs=4.0,K=200,omega=0.224,nfock=110,prec=384,rprec=1024)
E0=cuspsolve(b).E
Er,slope=secant_root(E->match(b,E;p...).f,E0)
r=match(b,Er;p...)
println("root precision=",precision(Er)," f=",Float64(r.f)," slope=",Float64(slope)," f/slope=",Float64(r.f/slope))
println("E=",Er," a21=",r.a21)
q=merge(p,(Rs=5.0,nfock=135,))
fv=match(b,Er;q...).f; Ev=Er-fv/slope; rv=match(b,Ev;q...)
h=big"1e-30"; sv=(match(b,Ev+h;q...).f-match(b,Ev-h;q...).f)/(2h)
println("R variant dE=",Float64(Ev-Er)," slope_ratio=",Float64(sv/slope)," remaining correction=",Float64(rv.f/sv))
