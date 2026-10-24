include("constrained.jl")
function corrected(path,nq)
 lines=readlines(path);E=parse(BigFloat,split(lines[1]," = ")[2])
 bas=[Tuple(parse.(Int,split(l)[1:3])) for l in lines[2:end]];c=[parse(BigFloat,split(l)[4]) for l in lines[2:end]]
 z,w=gauss(nq);I=big(0);nc=big(0)
 maxn=maximum(b[1] for b in bas)
 rad=Dict((m,p)=>moment(m,p,ee) for m=0:maxn for p=0:1)
 for i=1:nq,j=1:nq
  y=z[i];x=y*z[j];q=(1+x*x)/2;v=1-y*y/q;V=-8/(1-x*x)+1/y
  wt=64/big(pi)*y*y*(1-x*x)/(1+x*x)^3*v*w[i]*w[j];nc+=wt*v
  val=big(0)
  for (b,coef) in zip(bas,c),(n,p,t,_,_) in terms(b,x,y)
   val+=coef*t*(n==0 ? V*ee*(log(ee)+gamma-1-log(q)/2)+E*(gamma+log(ee)-log(q)/2) : n==1 ? V*(-gamma-log(ee)+log(q)/2)-E/ee : V*rad[(n-2,p)]-E*rad[(n-1,p)])
  end
  I+=wt*val
 end
 a=-I/4+2kap/8
 println("CORRECTED $(basename(path)) nq=$nq local=$(c[2]/c[1]) global=$(a/c[1]) norm=$nc")
end
if abspath(PROGRAM_FILE)==@__FILE__
 for path in ARGS;corrected(path,24);corrected(path,48);end
end
