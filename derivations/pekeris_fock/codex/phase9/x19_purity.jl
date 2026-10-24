include("x19_atoms.jl")
setprecision(BigFloat,256)
co=coefficients41(BigFloat);πb=BigFloat(π)
println("Independent purity quadrature on the square X=2(α-β)/π, Y=2(α+β-π)/π")
results=[]
for n in (16,24,32,48)
    z,w=gaussleg(n,BigFloat);v=zeros(BigFloat,2)
    for i=1:n,j=1:n
        a=πb*(2+z[i]+z[j])/4;b=πb*(2-z[i]+z[j])/4
        x=cos(a/2);y=sin(a/2);c=sqrt(2)*sin(b/2)
        Y0=4cos(a)^2-1;Y2=3cos(b)^2/2-sin(a)^2/2
        f=form41(x,y,c,co)
        wt=w[i]*w[j]*sin(a)*sin(b)*πb^2/8
        v .+= wt*f.*[Y0,Y2]
    end
    # Sphere norms with the ignorable azimuth removed: π, π/8.
    h=[-v[1]/πb,-8v[2]/πb];push!(results,h)
    println("n=",n," ah=",h[1]," bh=",h[2]);flush(stdout)
end
hc=parse.(BigFloat,split(read(joinpath(X19DIR,"x19_claude_harmonics.txt"),String)))
println("Claude fitted coefficients: ",hc)
println("differences from independent n48 projection: ",results[end]-hc)
# Fresh points scattered through spec41, no fitting beyond the two independently projected constants.
lines=readlines(joinpath(X19DIR,"../../claude9/spec41.txt"));errs=zeros(BigFloat,2)
for i=13:37:min(length(lines),3000)
    a,t,v=parse.(BigFloat,split(lines[i]));x,y,c=geometry(a,t)
    for (j,h) in enumerate((hc,results[end]))
        errs[j]=max(errs[j],abs(form41(x,y,c,co;ah=h[1],bh=h[2])-v))
    end
end
println("held-out max errors [Claude harmonics, pure quadrature harmonics]: ",errs)
open(joinpath(X19DIR,"x19_purity_results.txt"),"w") do io
    for (n,h) in zip((16,24,32,48),results);println(io,n," ",h[1]," ",h[2]);end
    println(io,"heldout ",errs[1]," ",errs[2])
end

@assert maximum(abs.(results[4]-results[3]))<big"1e-48"
@assert errs[1]<big"1e-39"
println("PASS independent purity convergence and spectral comparison; peak_rss_bytes=",Sys.maxrss())
# Independent Feynman data have their own harmonic gauge; fit exactly those two constants on 80 points, test all others.
D=[parse.(BigFloat,split(l)) for l in readlines(joinpath(X19DIR,"../../claude9/data2.txt"))];H=zeros(BigFloat,length(D),2);rhs=zeros(BigFloat,length(D))
for (i,row) in enumerate(D)
    a,t,v=row; x,y,c=geometry(a,t);H[i,:]=[4cos(a)^2-1,sin(a)^2*(3cos(t)^2-1)/2];rhs[i]=v-form41(x,y,c,co)
end
hfit=qr(H[1:80,:])\rhs[1:80];rf=H*hfit-rhs
println("Feynman data2 only-two-harmonics fit: ",hfit," held max=",maximum(abs,rf[81:end]))
@assert maximum(abs,rf[81:end])<big"1e-42"
println("PASS purity, spectral and independent Feynman checks; peak_rss_bytes=",Sys.maxrss())
