# Reuse our own definitions without executing the preceding sweep.
include_string(Main,replace(read(joinpath(@__DIR__,"x19_fit40_big.jl"),String),r"\nmain\(\)\s*$"=>"\n"),joinpath(@__DIR__,"x19_fit40_big.jl"))
function paritybasis(A)
    E=hcat(A[:,1],A[:,2],A[:,4]);O=hcat(A[:,3].-BigFloat(π)/2,A[:,5])
    for (i,j,sgn) in ((6,12,-1),(7,11,-1),(8,10,-1),(9,13,-1),(15,16,1),(17,18,1))
        E=hcat(E,A[:,i]+sgn*A[:,j]);O=hcat(O,A[:,i]-sgn*A[:,j])
    end
    E=hcat(E,A[:,[14,20,21,22,23]]);O=hcat(O,A[:,19])
    for (i,j) in ((24,47),(70,93))
        E=hcat(E,A[:,i:i+22]+A[:,j:j+22]);O=hcat(O,A[:,i:i+22]-A[:,j:j+22])
    end
    @assert size(E,2)==60 && size(O,2)==55;E,O
end
function singlet()
    A=readdlm(joinpath(X19DIR,"x19_fit40_big_atoms.txt"),BigFloat)
    ids=parse.(Int,readlines(joinpath(X19DIR,"x19_fit40_big_ids.txt")));
    lines=readlines(joinpath(X19DIR,"../../claude9/spec40.txt"));data=[parse.(BigFloat,split(lines[i])) for i in ids]
    geo=[geometry(r[1],r[2]) for r in data];N=length(geo);ntr=800
    E,O=paritybasis(A)
    # Verify the assumed exchange eigenvalues using actual independent function evaluations.
    worst=big"0"
    for i in (1,13)
        x,y,c=geo[i];Ag=reshape(vcat(order2(y,x,c),classical(y,x,c),bigpanels(y,x,c)),1,:);ee,oo=paritybasis(Ag)
        worst=max(worst,maximum(abs.((ee[1,:]-E[i,:])./(1 .+abs.(E[i,:])))),maximum(abs.((oo[1,:]+O[i,:])./(1 .+abs.(O[i,:])))))
    end
    @assert worst<big"1e-51";println("PASS exchange eigenbasis scaled discrepancy=",worst)
    n=2;em=[(i,k) for i=0:n for k=0:n if i+k<=n];om=[(i,k) for i=0:n-1 for k=0:n-1 if i+k<=n-1]
    M=zeros(BigFloat,N,60length(em)+55length(om))
    for row=1:N
        x,y,c=geo[row];s=x+y;d=x-y
        for a=1:60,(j,(i,k)) in enumerate(em);M[row,(a-1)*length(em)+j]=E[row,a]*s^i*c^k;end
        off=60length(em)
        for a=1:55,(j,(i,k)) in enumerate(om);M[row,off+(a-1)*length(om)+j]=O[row,a]*d*s^i*c^k;end
    end
    @assert size(M,2)==525
    Y=reduce(vcat,[permutedims(r[4:6]) for r in data]);delta=[4x*x*y*y-(1-c*c)^2 for (x,y,c) in geo];clear=[x^3*y^3*c^3*(2-c*c)^2*d^2 for ((x,y,c),d) in zip(geo,delta)]
    # A bounded known panel-module member is a separate positive control.
    coeff=[BigFloat(mod(j*19,17)-8)/1000 for j=1:size(M,2)];control=M*coeff
    target=hcat(Y,Y.*delta,Y.*delta.^2,Y.*clear,control)
    for j=1:size(M,2);sc=norm(M[1:ntr,j]);sc>0 && (@views M[:,j]./=sc);end
    A=nothing;E=nothing;O=nothing;GC.gc();println("QR full singlet d2 dims=",size(M)," preQR peak=",Sys.maxrss());flush(stdout)
    F=qr(M[1:ntr,:],ColumnNorm());dd=abs.(diag(F.R));r=findfirst(<(dd[1]*big"1e-42"),dd);r=isnothing(r) ? length(dd) : r-1
    qq=F.Q'*target[1:ntr,:];co=F.R[1:r,1:r]\qq[1:r,:];res=M[:,F.p[1:r]]*co-target
    open(joinpath(X19DIR,"x19_fit40_singlet_results.txt"),"w") do io
        println(io,"# degree columns QR_numerical_rank rhs train_relative_l2 held_relative_l2 held_relative_max")
        for j=1:size(target,2)
            tr=norm(res[1:ntr,j])/norm(target[1:ntr,j]);te=norm(res[ntr+1:end,j])/norm(target[ntr+1:end,j]);mx=maximum(abs,res[ntr+1:end,j])/maximum(abs,target[ntr+1:end,j])
            println(io,join((n,size(M,2),r,j,tr,te,mx)," "));println("rhs=",j," train=",Float64(tr)," held=",Float64(te));flush(stdout)
            j==13 && @assert te<big"1e-35"
        end
    end
    println("PASS all115 full singlet d2; peak_rss_bytes=",Sys.maxrss())
end
singlet()
