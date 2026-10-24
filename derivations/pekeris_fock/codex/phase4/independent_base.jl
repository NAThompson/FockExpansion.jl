# Own phase-3 actual-letter formulas and radial moments, reduced for phase 4.
# No Claude implementation is included.
using LinearAlgebra, Printf, ForwardDiff
import ClausenFunctions
setprecision(BigFloat,256)
function ClausenFunctions.cl2(d::ForwardDiff.Dual{T}) where T
 t=ForwardDiff.value(d)
 ForwardDiff.Dual{T}(ClausenFunctions.cl2(t),-log(abs(2sin(t/2)))*ForwardDiff.partials(d))
end
function atom(g,x,y)
 g==:one && return one(x)
 q=(1+x*x)/2; al=big(pi)/2-2atan(x); v=1-y*y/q; P=sqrt(2q-y*y);be=asin(v)
 g==:alp && return -2atan(x)
 g==:Lt && return log(1+y)-log(q)/big(pi)
 g==:V && return P/y*(be+big(pi)/2)
 cl=ClausenFunctions.cl2; L5=y*P+x; t=(1-x*x)/2
 dl=cl(al-be)+cl(big(pi)+be-al)-cl(big(pi)-al-be)-cl(al+be)
 om=be*log(t*t/(L5*L5))+(be+al)*log((1-y*y)/t)+(be-al)*log((y*y-x*x)/t)-dl
 g==:U && return om+big(pi)*log((x+y)^2/L5)
 error(g)
end
function gauss(n)
 x,v=eigen(SymTridiagonal(zeros(n),[k/sqrt(4k*k-1) for k=1:n-1])); (BigFloat.((x.+1)./2),BigFloat.(v[1,:].^2))
end
const ee=big"1.704";const gamma=BigFloat(Base.MathConstants.eulergamma)
function moment(m,p,rate)
 z=sum((inv(BigFloat(j)) for j=1:m);init=big(0))-gamma-log(rate)
 a=big(pi)^2/6-sum((inv(BigFloat(j)^2) for j=1:m);init=big(0))
 z3=5sum((-1)^(j+1)/(BigFloat(j)^3*binomial(big(2j),big(j))) for j=1:150)/2
 b=-2(z3-sum((inv(BigFloat(j)^3) for j=1:m);init=big(0)))
 c=6(big(pi)^4/90-sum((inv(BigFloat(j)^4) for j=1:m);init=big(0)))
 factorial(big(m))/rate^(m+1)*(p==0 ? one(z) : p==1 ? z : p==2 ? z*z+a : p==3 ? z^3+3z*a+b : z^4+6z*z*a+3a*a+4z*b+c)
end
