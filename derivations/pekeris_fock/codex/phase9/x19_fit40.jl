include("x19_atoms.jl")
# Serial, bounded dense least squares: at most ~3000 x 1840 Float64 coefficients.
function mons(x,y,c,n)
    [(x^i*y^j*c^k) for j=0:1 for i=0:n for k=0:n if i+j+k<=n]
end
function fitmatrix(geo,A,n)
    nmono=length(mons(geo[1]...,n));M=zeros(length(geo),size(A,2)*nmono)
    for i=1:length(geo)
        m=mons(geo[i]...,n)
        for a=1:size(A,2),j=1:nmono;M[i,(a-1)*nmono+j]=A[i,a]*m[j];end
    end
    M
end
function dorun(geo,A,Y,n,ntr,mult,label,io;tol=1e-12)
    M=fitmatrix(geo,A,n);sc=[norm(@view M[1:ntr,j]) for j=1:size(M,2)]
    for j=1:size(M,2);sc[j]>0 && (@views M[:,j]./=sc[j]);end
    y=Y.*mult;F=svd(M[1:ntr,:];full=false);r=count(>(tol*F.S[1]),F.S)
    co=F.V[:,1:r]*((F.U[:,1:r]'*y[1:ntr,:])./F.S[1:r])
    tr=M[1:ntr,:]*co-y[1:ntr,:];te=M[ntr+1:end,:]*co-y[ntr+1:end,:]
    for j=1:3
        a=norm(tr[:,j])/norm(y[1:ntr,j]);b=norm(te[:,j])/norm(y[ntr+1:end,j]);mx=maximum(abs,te[:,j])/maximum(abs,y[ntr+1:end,j])
        println(io,join((label,n,size(M,2),r,tol,j,a,b,mx)," "));flush(io)
        @printf("%s d=%d cols=%d rank=%d cut=%.0e C%d train=%.3e held=%.3e max=%.3e\n",label,n,size(M,2),r,tol,j,a,b,mx)
    end
    flush(stdout);nothing
end
function main()
    data=[parse.(Float64,split(l)) for l in readlines(joinpath(X19DIR,"../../claude9/spec40.txt"))]
    # No point selection by fitted values; fixed compact geometry for reliable Float64 quadrature.
    ids=findall(r->0.25<r[1]<π-0.25 && 0.25<r[2]<π-0.25,data)
    geo=[geometry(data[i][1],data[i][2]) for i in ids];Y=reduce(vcat,[permutedims(data[i][4:6]) for i in ids]);N=length(ids);ntr=floor(Int,N*0.8)
    println("points ",N," training ",ntr," held-out ",N-ntr)
    z,w=gaussleg(256);ss=(z.+1)./2;ws=w./2
    A=reduce(vcat,[permutedims(vcat(order2(g...),classical(g...),panels(g...,ss,ws))) for g in geo]);@assert size(A,2)==115
    z2,w2=gaussleg(512);s2=(z2.+1)./2;w2=w2./2;err=0.0
    for i=1:N
        ref=panels(geo[i]...,s2,w2);err=max(err,maximum(abs.((A[i,24:end]-ref)./(1 .+abs.(ref)))))
    end
    println("panel 256/512 Gauss relative-scale maximum discrepancy ",err);@assert err<1e-10
    # Delimited caches preserve input ids and evaluator values; no later evaluator imports required.
    writedlm(joinpath(X19DIR,"x19_fit40_ids.txt"),ids)
    writedlm(joinpath(X19DIR,"x19_fit40_atoms.txt"),A)
    open(joinpath(X19DIR,"x19_fit40_results.txt"),"w") do io
        println(io,"# family degree columns numerical_rank svd_cut charge train_relative_l2 held_relative_l2 held_relative_max")
        for (label,sel,degs) in (("order2",1:5,(0,2,4,6)),("order2_classical",1:23,(0,1,2,3,4)),("all115",1:115,(0,1,2,3)))
            for n in degs;dorun(geo,A[:,sel],Y,n,ntr,ones(N),label,io);end
        end
        delta=[4x*x*y*y-(1-c*c)^2 for (x,y,c) in geo]
        # These include every irreducible denominator factor in x16 joint row clearings.
        clear=[x^3*y^3*c^3*(2-c*c)^2*d^2 for ((x,y,c),d) in zip(geo,delta)]
        for (name,mult) in (("all115_delta",delta),("all115_delta2",delta.^2),("all115_commonclear",clear))
            for n in (1,2,3);dorun(geo,A,Y,n,ntr,mult,name,io);end
        end
        # Rank-cut sensitivity is necessary for a highly redundant, ill-conditioned atom family.
        for tol in (1e-10,1e-14);dorun(geo,A,Y,2,ntr,ones(N),"all115_cut",io;tol=tol);end
        # Positive control: known psi41 is a degree-4 polynomial combination of order-two atoms.
        co=coefficients41(Float64);Z=reshape([form41(g...,co) for g in geo],:,1)
        dorun(geo,A[:,1:5],repeat(Z,1,3),4,ntr,ones(N),"control41",io)
    end
    println("quadrature_max=",err," atoms=",size(A)," maxrss=",Sys.maxrss())
end
main()
