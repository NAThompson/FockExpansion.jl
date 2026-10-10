using Nemo, JSON3
const PR, VX = polynomial_ring(QQ,["T","U","R","P","G","L","J","b","y","t"])
const MODR=VX[3]^2-2; const MODI=VX[7]^2+1
red(p)=last(divrem(last(divrem(p,MODR)),MODI))
struct Rat
 n::QQMPolyRingElem
 d::QQMPolyRingElem
 function Rat(n::QQMPolyRingElem,d::QQMPolyRingElem=one(PR))
  n=red(n);d=red(d);iszero(d) && error("zero denominator")
  iszero(n) && return new(zero(PR),one(PR))
  g=gcd(n,d);new(divexact(n,g),divexact(d,g))
 end
end
Rat(n::Integer)=Rat(PR(n))
Rat(n::Rational)=Rat(PR(QQ(n)))
Base.inv(a::Rat)=Rat(a.d,a.n)
Base.iszero(a::Rat)=iszero(a.n)
Base.:-(a::Rat)=Rat(-a.n,a.d)
function Base.:+(a::Rat,b::Rat)
 iszero(a) && return b;iszero(b) && return a
 g=gcd(a.d,b.d);ad=divexact(a.d,g);bd=divexact(b.d,g)
 Rat(a.n*bd+b.n*ad,ad*b.d)
end
Base.:-(a::Rat,b::Rat)=a+(-b)
Base.:*(a::Rat,b::Rat)=Rat(a.n*b.n,a.d*b.d)
Base.:/(a::Rat,b::Rat)=Rat(a.n*b.d,a.d*b.n)
Base.:^(a::Rat,k::Integer)=k<0 ? Rat(a.d^(-k),a.n^(-k)) : Rat(a.n^k,a.d^k)
for op in (:+,:-,:*,:/)
 @eval Base.$op(a::Rat,b::Union{Integer,Rational})=Base.$op(a,Rat(b))
 @eval Base.$op(a::Union{Integer,Rational},b::Rat)=Base.$op(Rat(a),b)
end
rdiff(a::Rat,j)=Rat(derivative(a.n,j)*a.d-a.n*derivative(a.d,j),a.d^2)
const rv=Rat.(VX)
# Sparse polynomial in angle, logarithm and Clausen atoms; rational coefficients.
struct Form
 terms::Dict{Tuple{Vararg{String}},Rat}
end
Form(a::Rat)=Form(iszero(a) ? Dict{Tuple{Vararg{String}},Rat}() : Dict{Tuple{Vararg{String}},Rat}(()=>a))
Form(a::Union{Integer,Rational})=Form(Rat(a))
atom(s)=Form(Dict{Tuple{Vararg{String}},Rat}((String(s),)=>Rat(1)))
Base.iszero(a::Form)=isempty(a.terms)
function put!(d,k,v)
 v=get(d,k,Rat(0))+v
 if iszero(v);delete!(d,k);else;d[k]=v;end
end
function Base.:+(a::Form,b::Form)
 d=copy(a.terms);for (k,v) in b.terms;put!(d,k,v);end;Form(d)
end
Base.:-(a::Form)=Form(Dict(k=>-v for (k,v) in a.terms))
Base.:-(a::Form,b::Form)=a+(-b)
function Base.:*(a::Form,b::Form)
 d=Dict{Tuple{Vararg{String}},Rat}()
 for (k,v) in a.terms,(l,w) in b.terms
  put!(d,Tuple(sort!([k...;l...])),v*w)
 end
 Form(d)
end
function scalar(a::Form)
 all(isempty,keys(a.terms)) || error("not a rational coefficient")
 get(a.terms,(),Rat(0))
end
Base.:/(a::Form,b::Form)=a*Form(scalar(b)^(-1))
function Base.:^(a::Form,k::Integer)
 k<0 && return Form(scalar(a)^k)
 ans=Form(1);for _=1:k;ans=ans*a;end;ans
end
for op in (:+,:-,:*,:/)
 @eval Base.$op(a::Form,b::Union{Rat,Integer,Rational})=Base.$op(a,Form(b))
 @eval Base.$op(a::Union{Rat,Integer,Rational},b::Form)=Base.$op(Form(a),b)
end
function lam(m,n,k)
 if m==n==0
  mod(k,8) in (2,6) && return Form(rv[6]/2)
  mod(k,8)==4 && return Form(rv[6])
  error("unsupported constant logarithm $k")
 end
 if all(iseven,(m,n,k));return lam(m÷2,n÷2,k÷2)+lam(m÷2,n÷2,k÷2+4);end
 if (m,n)<(0,0);m,n,k=-m,-n,-k;end
 atom("L_$(m)_$(n)_$(mod(k,8))")
end
function cl(m,n,k)
 if m==n==0;return mod(k,8)==0 ? Form(0) : atom("Clconst$(mod(k,16))");end
 sg=1
 if (m,n)<(0,0);m,n,k=-m,-n,-k;sg=-1;end
 sg*atom("C_$(m)_$(n)_$(mod(k,16))")
end
const COT=Dict{Tuple{Int,Int,Int},Rat}()
function cotq(m,n,k)
 get!(COT,(m,n,k)) do
  T,U,R,P,G,L,I=rv[1:7]
  phase=((R+R*I)/2)^mod(k,8)
  z=phase*((1+I*T)/(1-I*T))^m*((1+I*U)/(1-I*U))^n
  I*(z+1)/(z-1)
 end
end
function atomdiff(s,j)
 s=="AL" && return Form(j==1 ? 1 : 0)
 s=="PS" && return Form(j==2 ? 1 : 0)
 startswith(s,"Clconst") && return Form(0)
 kind,ms,ns,ks=split(s,"_");m,n,k=parse.(Int,(ms,ns,ks));c=(j==1 ? m : n)//4
 iszero(c) && return Form(0)
 kind=="L" && return Form(c*cotq(m,n,k))
 @assert all(iseven,(m,n,k))
 -c*lam(m÷2,n÷2,k÷2)
end
function fdiff(a::Form,j)
 out=Form(0)
 for (k,v) in a.terms
  out+=Form(Dict(k=>rdiff(v,j)*(1+rv[j]^2)/4))
  for i=1:length(k)
   rest=Tuple(k[q] for q=1:length(k) if q!=i)
   out+=Form(Dict(rest=>v))*atomdiff(k[i],j)
  end
 end
 # Some dictionary constructors retain zero entries.
 filter!(p->!iszero(p.second),out.terms);out
end
function duplicate(a::Form)
 out=Form(0)
 for (key,c) in a.terms
  term=Form(c)
  for s in key
   if startswith(s,"C_")
    m,n,k=parse.(Int,split(s,"_")[2:4])
    if max(abs(m),abs(n))>=4 && all(iseven,(m,n,k))
     term*=2cl(m÷2,n÷2,k÷2)+2cl(m÷2,n÷2,k÷2+8);continue
    end
   end
   term*=atom(s)
  end
  out+=term
 end
 out
end
function load_inputs(path)
 data=JSON3.read(read(path,String));values=Form[]
 names=Dict(string(VX[i])=>Form(rv[i]) for i=1:length(VX))
 for node in data.nodes
  op=String(node[1])
  v=if op=="symbol"
   get(names,String(node[2])) do;atom(String(node[2]));end
  elseif op=="rational"
   Form(Rat(PR(QQ(parse(BigInt,String(node[2])),parse(BigInt,String(node[3]))))))
  elseif op=="pow";values[Int(node[2])]^Int(node[3])
  elseif op=="add";foldl(+, (values[Int(i)] for i in node[2:end]);init=Form(0))
  elseif op=="mul";foldl(*, (values[Int(i)] for i in node[2:end]);init=Form(1))
  else;error(op)
  end
  push!(values,v)
 end
 Dict(String(k)=>values[Int(v)] for (k,v) in pairs(data.roots))
end
function rsubst(f::Rat,b,y)
 vals=copy(rv);vals[8]=b;vals[9]=y
 function peval(p)
  out=Rat(0)
  for (c,es) in zip(coefficients(p),exponent_vectors(p))
   term=Rat(PR(c))
   for i=1:length(es);iszero(es[i]) || (term*=vals[i]^es[i]);end
   out+=term
  end
  out
 end
 peval(f.n)/peval(f.d)
end
