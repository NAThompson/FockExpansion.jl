include("x19_atoms.jl")
setprecision(BigFloat,256)
const G96=gaussleg(96,BigFloat)
const G128=gaussleg(128,BigFloat)
function adapt_panel(H,k;rule=G96)
    h=sqrt(H);pts=BigFloat[0,1]
    # Resolve both the imaginary arctan poles and the near-real quartic roots.
    if h<1
        a=h
        while a<1;push!(pts,a);a*=2;end
    end
    if k<0
        s0=sqrt(-k/2)
        if s0<1
            push!(pts,s0)
            width=sqrt(H-k*k/4)/(2s0+sqrt(sqrt(H-k*k/4)))
            d=width
            while d<2
                0<s0-d<1 && push!(pts,s0-d)
                0<s0+d<1 && push!(pts,s0+d)
                d*=2
            end
        end
    end
    sort!(unique!(pts));z,w=rule;ss=BigFloat[];ws=BigFloat[]
    for i=1:length(pts)-1
        a,b=pts[i:i+1];append!(ss,(a+b)/2 .+(b-a)/2 .*z);append!(ws,(b-a)/2 .*w)
    end
    panel(H,k,ss,ws)
end
function bigpanels(x,y,c;rule=G96)
    Q=x*x+y*y;e=2Q-c*c
    vcat(adapt_panel(x*x/(y*y),(c*c-Q)/(y*y);rule),adapt_panel(y*y/(x*x),(c*c-Q)/(x*x);rule),adapt_panel(c*c/e,2(x*x-y*y)/e;rule),adapt_panel(c*c/e,2(y*y-x*x)/e;rule))
end
function monsb(x,y,c,n)
    [x^i*y^j*c^k for j=0:1 for i=0:n for k=0:n if i+j+k<=n]
end
function bigfit(geo,A,Y,n,ntr,label,io)
    nm=length(monsb(geo[1]...,n));M=zeros(BigFloat,length(geo),nm*size(A,2))
    for i=1:length(geo)
        m=monsb(geo[i]...,n)
        for a=1:size(A,2),j=1:nm;M[i,(a-1)*nm+j]=A[i,a]*m[j];end
    end
    for j=1:size(M,2);s=norm(M[1:ntr,j]);s>0 && (@views M[:,j]./=s);end
    println("QR ",label," degree=",n," dims=",size(M));flush(stdout)
    F=qr(M[1:ntr,:],ColumnNorm());d=abs.(diag(F.R));r=findfirst(<(d[1]*big"1e-42"),d);r=isnothing(r) ? length(d) : r-1
    yy=F.Q'*Y[1:ntr,:];co=F.R[1:r,1:r]\yy[1:r,:];res=M[:,F.p[1:r]]*co-Y
    for j=1:size(Y,2)
        tr=norm(res[1:ntr,j])/norm(Y[1:ntr,j]);te=norm(res[ntr+1:end,j])/norm(Y[ntr+1:end,j]);mx=maximum(abs,res[ntr+1:end,j])/maximum(abs,Y[ntr+1:end,j])
        println(io,join((label,n,size(M,2),r,j,tr,te,mx)," "));flush(io)
        println(label," d=",n," rank=",r," rhs=",j," train=",Float64(tr)," held=",Float64(te));flush(stdout)
    end
end
function main()
    lines=readlines(joinpath(X19DIR,"../../claude9/spec40.txt"));allids=parse.(Int,readlines(joinpath(X19DIR,"x19_fit40_ids.txt")))
    ids=allids[1:2:2000];N=length(ids);ntr=800
    data=[parse.(BigFloat,split(lines[i])) for i in ids];geo=[geometry(r[1],r[2]) for r in data]
    A=zeros(BigFloat,N,115);err=big"0"
    for i=1:N
        g=geo[i];p=bigpanels(g...)
        A[i,:]=vcat(order2(g...),classical(g...),p)
        if i in (1,13,257,777,1000)
            err=max(err,maximum(abs.((p-bigpanels(g...;rule=G128))./(1 .+abs.(p)))))
        end
        i%200==0 && (println("evaluated ",i," points; quadrature discrepancy ",Float64(err));flush(stdout))
    end
    @assert err<big"1e-52"
    writedlm(joinpath(X19DIR,"x19_fit40_big_atoms.txt"),A)
    writedlm(joinpath(X19DIR,"x19_fit40_big_ids.txt"),ids)
    Y=reduce(vcat,[permutedims(r[4:6]) for r in data]);delta=[4x*x*y*y-(1-c*c)^2 for (x,y,c) in geo]
    clear=[x^3*y^3*c^3*(2-c*c)^2*d^2 for ((x,y,c),d) in zip(geo,delta)]
    targets=hcat(Y,Y.*delta,Y.*delta.^2,Y.*clear)
    open(joinpath(X19DIR,"x19_fit40_big_results.txt"),"w") do io
        println(io,"# family degree columns pivoted_QR_rank rhs train_relative_l2 held_relative_l2 held_relative_max")
        for n in (0,1);bigfit(geo,A,targets,n,ntr,"all115_big",io);end
        for n in (2,4,6);bigfit(geo,A[:,1:5],targets,n,ntr,"order2_big",io);end
        co=coefficients41(BigFloat);Z=reshape([form41(g...,co) for g in geo],:,1)
        bigfit(geo,A[:,1:5],Z,4,ntr,"control41_big",io)
    end
    println("PASS quadrature err=",err," peak_rss_bytes=",Sys.maxrss())
end
main()
