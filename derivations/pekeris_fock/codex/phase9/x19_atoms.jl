# Independent evaluations of the Codex atoms; no Claude executable imports.
using LinearAlgebra, DelimitedFiles, Printf
using ClausenFunctions: cl2
BLAS.set_num_threads(1)
const X19DIR=@__DIR__
function gaussleg(n, ::Type{T}=Float64) where T
    # Newton construction works for Float64 and BigFloat, no package/precompile needed.
    xx=zeros(T,n);ww=zeros(T,n);P=T(π)
    for j=1:cld(n,2)
        z=cos(P*(T(j)-T(1)/4)/(T(n)+T(1)/2))
        dp=zero(T)
        for _=1:100
            p=one(T);p0=zero(T)
            for k=1:n
                p,p0=((2k-1)*z*p-(k-1)*p0)/k,p
            end
            dp=n*(z*p-p0)/(z*z-1); dz=p/dp;z-=dz
            abs(dz)<eps(T)*10 && break
        end
        xx[j]=-z;xx[n+1-j]=z;ww[j]=ww[n+1-j]=2/((1-z*z)*dp*dp)
    end
    xx,ww
end
function order2(x,y,c)
    Q=x*x+y*y;P=sqrt(2Q-c*c);sa=2x*y/Q;ca=(x*x-y*y)/Q;v=(Q-c*c)/Q;g=c*P/Q
    a=2atan(y/x);b=asin(clamp(v,-one(v),one(v)))
    aa=log(abs((g-ca)/(g+ca)));bb=log(abs((sa+v)/(sa-v)))
    jj=log((ca+(x+y)*c/Q)^2/((1+sa)*(g+ca)))
    dl=cl2(a-b)+cl2(oftype(a,π)+b-a)-cl2(oftype(a,π)-a-b)-cl2(a+b)
    U=b*aa+a*bb-dl+oftype(a,π)*jj
    [one(x),log((x+y+c)/sqrt(Q)),a,P/c*(b+oftype(a,π)/2),U]
end
function classical(x,y,c)
    P=oftype(x,π);al=2atan(y/x);ps=2asin(c/sqrt(2(x*x+y*y)))
    a=al-P/2;p=ps-P;wp=(a+p)/4;wm=(a-p)/4
    L(z)=iszero(z) ? zero(z) : z*log(abs(2sin(z)))+cl2(2z)/2
    T(z)=-z*log(abs(2cos(z)))+cl2(P-2z)/2
    v=typeof(x)[]
    for z in (wp,wm)
        append!(v,[L(z-P/4),L(z),L(z+P/4),T(z)])
    end
    append!(v,[L(p/2),L(al/2),L((P-al)/2),log(abs(2cos(wp))),log(abs(2cos(wm))),a,p,a*a,p*p,one(x)])
    v
end
const PNAMES=["J_A0","J_A1","J_A2","J_B2","K_A0","K_A1","Q0_A0","Q0_A1","Q0_A2","Q0_B2","Q2_A0","Q2_A1","Q2_A2","Q2_B2","U0_A0","U1_A0","U1_A1","V0_A0","V1_B1","FP_B1","E_B1","E_B2","F0"]
function panel(H,k,ss,ws)
    # 23 reachable x16 states, defining fixed-endpoint integrals.
    h=sqrt(H);out=zeros(typeof(H),23)
    for (s,w) in zip(ss,ws)
        A=atan(s);B=h*atan(s/h);T=s^4+k*s*s+H;f=w*h/sqrt(T)
        vals=(1,A,A*A/2,B*B/2,s,s*A,1/T,A/T,A*A/(2T),B*B/(2T),s*s/T,s*s*A/T,s*s*A*A/(2T),s*s*B*B/(2T),1/(1+s*s),s/(1+s*s),s*A/(1+s*s),1/(H+s*s),s*B/(H+s*s),B/s)
        for j=1:20;out[j]+=f*vals[j];end
    end
    B=h*atan(1/h);f=h/sqrt(1+k+H)
    out[21]=f*B;out[22]=f*B*B/2;out[23]=1;out
end
function panels(x,y,c,ss,ws)
    Q=x*x+y*y;e=2Q-c*c
    vcat(panel(x*x/(y*y),(c*c-Q)/(y*y),ss,ws),panel(y*y/(x*x),(c*c-Q)/(x*x),ss,ws),panel(c*c/e,2(x*x-y*y)/e,ss,ws),panel(c*c/e,2(y*y-x*x)/e,ss,ws))
end
function geometry(a,t)
    x=cos(a/2);y=sin(a/2);c=sqrt(2sin(oftype(a,π)/4-a/2)^2+2sin(a)*sin(t/2)^2);x,y,c
end
function coefficients41(::Type{T}) where T
    [(parse(Int,r[1]),parse(Int,r[2]),parse(Int,r[3]),parse(Int,r[4]),parse(T,r[5])) for r in split.(readlines(joinpath(X19DIR,"x19_closed41_coefficients.txt")))]
end
function form41(x,y,c,co;ah=zero(x),bh=zero(x))
    A=order2(x,y,c);v=sum(k*x^i*y^j*c^l*A[a] for (a,i,j,l,k) in co)
    Q=x*x+y*y
    v+ah*(4(x*x-y*y)^2-Q*Q)+bh*(3(Q-c*c)^2/2-2x*x*y*y)
end
