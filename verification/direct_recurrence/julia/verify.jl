#!/usr/bin/env julia
isdefined(@__MODULE__, :ui_pass) || include("terminal_ui.jl")
include("exact_engine.jl")
const INPUT=load_inputs(joinpath(@__DIR__,"formula_inputs.json"))
function checkzero(x,label)
 if !iszero(x)
  ui_fail(label)
  error("NONZERO: $label")
 end
 get(ENV,"FOCK_CERT_VERBOSE","0")=="1" && ui_pass(label)
end
function check_panels()
 b,y,t=rv[8:10];P=1+b+(2-4y^2)*t^2+(1-b)*t^4
 N=(1+t^2)*(1+b+(1-b)*t^2)
 h=Dict(8=>-2y*t*(1-t^2)/(1+b+(1-b)*t^2),9=>4t,10=>4y*(1+b-(1-b)*t^4)/N)
 metric=((8,8,-4*(1-b^2)),(9,9,-(1-y^2)),(8,9,-2b*(1-2y^2)/y))
 drift=((8,12b),(9,-(2-5y^2)/y))
 D(q,j)=rdiff(q,j)-q*rdiff(P,j)/(2P)
 for name in ("angZ","angZ2","rotZ2")
  Q=scalar(INPUT[name*"_Q"]);S=scalar(INPUT[name*"_S"])
  manual=name=="angZ" ? (8b^2-8)/576+(b-1)*(-18b-32y^2+16)*t^2/288 : name=="angZ2" ? -(9b^2-18)/576+(b-1)*(90y^2-45)*t^2/288 : (1-b^2)/18+2*(b-1)*(2y^2-1)*t^2/9
  checkzero(Q-manual,name*" paper weight")
  R=-21Q;W=Rat(0)
  for (i,j,g) in metric
   R+=g*D(D(Q,i),j)
   W+=g*(h[i]*D(Q,j)+h[j]*D(Q,i)+Q*(rdiff(h[i],j)-h[i]*rdiff(P,j)/(2P)))/P
  end
  for (i,g) in drift;R+=g*D(Q,i);W+=g*h[i]*Q/P;end
  checkzero(R-scalar(INPUT[name*"_R"]),name*" differential action R")
  checkzero(W-scalar(INPUT[name*"_W"]),name*" differential action W")
  checkzero(R*P^2-rdiff(S,10)*P+3S*rdiff(P,10)/2,name*" integral primitive")
  H=S/P;h0=scalar(INPUT[name*"_h0"]);h2=scalar(INPUT[name*"_h2"])
  checkzero(H-h0*t-h2*t^3,name*" primitive coefficients")
  c=[scalar(INPUT[name*"_c"*string(i)]) for i=1:3];DD=1+b+(1-b)*t^2
  checkzero(W-H*h[10]/P-c[1]*t/(1+t^2)-c[2]*t/DD-c[3]*t/DD^2,name*" rational remainder")
 end
end
function image(name,i)
 getv(s)=INPUT["chart$(i)_"*s]
 b=scalar(getv("b"));y=scalar(getv("y"));P=rv[4];L=rv[6]
 coef(s)=rsubst(scalar(INPUT[name*"_"*s]),b,y)
 h0=coef("h0");h2=coef("h2");ab=getv("ab");be=getv("be")
 el=coef("c1")*P^2/16+coef("c2")*ab^2/(4*(1-b))+coef("c3")*ab/(4*getv("sinab"))
 2h2*getv("Kinf")+(h0+h2)*getv("K1")+(h0+h2)*(1+2be/P)/(2*getv("sy"))*L+el/P
end
function angular_operator(f)
 T,U=rv[1:2];r1=(1-T^2)/(1+T^2);r2=2T/(1+T^2);cp=(1-U^2)/(1+U^2);sp=2U/(1+U^2)
 a=r1^2-r2^2;v=cp^2-sp^2;sa=2r1*r2;sb=2cp*sp
 fa=fdiff(f,1);fp=fdiff(f,2)
 -4*(fdiff(fa,1)+fdiff(fp,2)-2a*v/(sa*sb)*fdiff(fa,2)+2a/sa*fa+2v/sb*fp)-21f
end
function verify()
 ui_stage("Integral identities")
 t0=time();check_panels();ui_pass("18 exact checks";seconds=time()-t0)
 ui_stage("Coefficient-table links")
 t0=time()
 for i=1:18
  checkzero(INPUT["link$(i)_left"]-INPUT["link$(i)_right"],"classical table link $i")
 end
 ui_pass("18/18 coefficient-table links";seconds=time()-t0)
 ui_stage("Recurrence sectors")
 for tag in ("Z0","Z1","Z2","Z3","E","EZ","a21","a21Z")
  t0=time()
  r=angular_operator(INPUT["formula_"*tag])-INPUT["source_"*tag]
  if tag in ("Z1","Z2")
   name=tag=="Z1" ? "angZ" : "angZ2"
   r+=rv[3]*(image(name,1)+image(name,2))
   tag=="Z2" && (r+=image("rotZ2",3)+image("rotZ2",4))
  end
  r=duplicate(r)
  checkzero(r,tag*" exact recurrence residual")
  ui_pass(tag;seconds=time()-t0)
  if tag=="Z0"
   @assert !iszero(r-1) # adding one to the source must fail
   @assert !iszero(r+angular_operator(Form(1))) # adding one to ψ must fail

  end
 end
 checkzero(angular_operator(Form(1))+21,"constant perturbation response")
 ui_pass("Altered formula and source rejected; constant response checked")
end
verify()
