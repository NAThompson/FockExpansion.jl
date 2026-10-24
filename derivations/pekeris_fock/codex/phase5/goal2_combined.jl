include("goal2_big_panel.jl")
# Coordinates a=alpha-big(pi)/2, p=p0-big(pi). Trigonometric coefficients are
# algebraic over physical distance coordinates on the chosen chart.
function classicalC(a,p,side,axis)
 C=zeros(Complex{BigFloat},18,18); pi0=big(pi)
 args=Tuple[]
 for s=(1,-1)
  w=(a+s*p)/4
  for m=(-1,0,1); push!(args,(w+m*big(pi)/4,1/4,s/4,0,m*big(pi)/4)); end
  push!(args,(w,1/4,s/4,1,0))
 end
 append!(args,[(p/2,0,1/2,0,0),((a+big(pi)/2)/2,1/2,0,0,big(pi)/4),((big(pi)/2-a)/2,-1/2,0,0,big(pi)/4)])
 for (i,(z,za,zp,typ,constant))=enumerate(args)
  fac=(axis==1 ? za : zp)*(typ==0 ? cot(z) : tan(z))
  C[i,14]=fac*za; C[i,15]=fac*zp; C[i,18]=fac*constant
 end
 for (i,s)=((12,1),(13,-1))
  w=(a+s*p)/4; C[i,18]=-tan(w)*(axis==1 ? 1/4 : s/4)
 end
 C[14,18]=axis==1 ? 1 : 0; C[15,18]=axis==2 ? 1 : 0
 C[16,14]=axis==1 ? 2 : 0; C[17,15]=axis==2 ? 2 : 0
 C
end
function classicalL(a,p,side)
 Z=2; al=a+big(pi)/2; ps=p+big(pi)
 r1=cos(al/2);r2=sin(al/2);sg=r1+r2;d=r1-r2;xi=sqrt(big(2))*sin(ps/2);eta0=sqrt(big(2))*cos(ps/2)
 out=zeros(Complex{BigFloat},1,18); idx=0
 for s=(1,-1)
  eta=s*eta0
  for m=(-1,0,1,2)
   dd,ss=d,sg
   for _=1:mod(m,4); dd,ss=-ss,dd; end
   kz=(dd*(-9ss^2+6ss*xi+8xi^2-8)-eta*(5-2xi^2))/(144big(pi))
   kz2=-5*(xi^2-1)*(dd+eta)/(32big(pi))
   nu=m==0 || m==2 ? xi^2-sg^2 : (9sg^2+xi^2-2+(m==1 ? -8 : 8)*sg*eta)
   nv=d*nu/(36big(pi))
   idx+=1; out[idx]=(m==2 ? -4 : 4)*(Z*kz+Z^2*(kz2+nv))
  end
  el=Z*(7sg*d*eta+9sg^3+8sg*xi^2-26sg+10xi^3-12xi)/72+Z^2*sg*(137d*eta+64sg^2-154xi^2+26)/288
  out[s==1 ? 12 : 13]=el
 end
 c9=-eta0*(Z*(2xi^2-5)/(9big(pi))-5Z^2*(xi^2-1)/(2big(pi)));out[9]=c9/2
 out[10]=2Z^2/(9big(pi))*r2*(8r2^2+xi^2-5)
 out[11]=2Z^2/(9big(pi))*r1*(8r1^2+xi^2-5)
 out[16]=Z*xi*(2xi^2+1)/(144big(pi))-5Z^2*xi*(xi^2-1)/(32big(pi))
 out[14]=7Z*xi*d*sg/144+Z^2*(xi*d*sg/64+d*(xi^2-1)/(9big(pi)))
 out[17]=Z^2*sg*(4sg^2+xi^2-9)/(36big(pi))
 out[15]=-Z*eta0*(5xi^2-4)/(36big(pi))+Z^2*sg*xi*eta0/(18big(pi))
 out
end

function maps(a,p)
 r1=cos(a/2+big(pi)/4);r2=sin(a/2+big(pi)/4)
 [(-sin(a),-cos(p/2),-cos(a),0,0,sin(p/2)/2),
  (sin(a),-cos(p/2),cos(a),0,0,sin(p/2)/2),
  (cos(p),-r1,0,r2/2,-sin(p),0),
  (cos(p),-r2,0,-r1/2,-sin(p),0)]
end
function combined(a,p,center,axis)
 n=332; C=spzeros(Complex{BigFloat},n,n);l=zeros(Complex{BigFloat},3,n)
 C[1:18,1:18]=classicalC(a,p,0,axis)
 l[1,18]=1;l[3,1:18]=classicalL(a,p,0)
 center_maps=maps(big(center[1]),big(center[2])); actual=maps(a,p);Z=2
 for j=1:4
  b,y,ba,ya,bp,yp=actual[j]; b0,y0=center_maps[j][1:2]
  basepoint[]=(Float64(b0),Float64(y0));side=sign(b0); sl=19+(j-1)*78:18+j*78
  db,dy=axis==1 ? (ba,ya) : (bp,yp)
  C[sl,sl]=iszero(dy) ? db*connection(b,y,side,1) : (iszero(db) ? dy*connection(b,y,side,2) : db*connection(b,y,side,1)+dy*connection(b,y,side,2))
  L=outputs(b,y,side)
  if j<=2
   c0=-sqrt(big(2))*Z*(9Z*b*b-18Z-8b*b+8)/576
   c2=-sqrt(big(2))*Z*(1-b)*(90Z*y*y-45Z-18b-32y*y+16)/288
  else
   c0=Z^2*(1-b*b)/18;c2=-2Z^2*(1-b)*(2y*y-1)/9
  end
  l[3,sl]=c0*L[1,:]+c2*L[2,:]
 end
 # Add the two old atoms not already in the classical alphabet: L and U.
 rt=sqrt(big(2));s=rt*cos(a/2);c=rt*cos(p/2);P=-rt*sin(p/2)
 sn=cos(a);cs=-sin(a);V=-cos(p);g=c*P;t=c*c;den=sn*sn-V*V;A=cs+s*c
 ds=axis==1 ? -rt*sin(a/2)/2 : 0
 dc=axis==2 ? -rt*sin(p/2)/2 : 0;dP=axis==2 ? -rt*cos(p/2)/2 : 0
 dsn=axis==1 ? -sin(a) : 0;dcs=axis==1 ? -cos(a) : 0;dV=axis==2 ? sin(p) : 0
 kg=dc/c+dP/P
 j0=2*(dcs+ds*c+s*dc)/A-dsn/(1+sn)-(g*g*kg-cs*dcs)/den
 tx=2*(sn*dV-V*dsn)/den;tb=2*t*(cs*kg-dcs)/den
 C[331,18]=(ds+dc)/(s+c)
 C[332,18]=big(pi)*j0+big(pi)/2*tx
 C[332,14]=tx;C[332,15]=-P/c*tb
 l[2,331]=-2big(pi)*V;l[2,15]=-2*c*P;l[2,332]=cs
 C,l
end
function physical_projection(a,p)
 T=zeros(Complex{BigFloat},332,240);T[1:18,1:18]=Matrix{Complex{BigFloat}}(I,18,18)
 panels=maps(a,p)
 for panel=1:4
  off=18+(panel-1)*78;coff=18+(panel-1)*55
  for fam=1:2,i=1:11,k=1:3
   col=coff+(i-1)*5+(k==1 ? 1 : (fam==1 ? k : k+2))
   T[off+(fam-1)*39+(i-1)*3+k,col]=1
  end
  b,y=panels[panel][1:2];basepoint[]=(Float64(b),Float64(y));side=sign(b)
  h,z,u=geom(complex(b),complex(y),side);e0=1/(-2*y*u);f0=1/sqrt(1+b)
  # Endpoint angle = c + aa*a + pp*p on these physical chambers.
  if panel<=2
   A=panel==1 ? (big(pi)/4,0,0) : (big(pi)/4,big(1)/2,0)
   D=panel==1 ? (big(pi)/4,big(1)/2,0) : (big(pi)/4,0,0)
  else
   A=side>0 ? (big(pi)/4,0,0) : (big(pi)/2,0,big(1)/2)
   D=side>0 ? (big(0),0,-big(1)/2) : (big(pi)/4,0,0)
  end
  for (fam,angle)=enumerate((A,D))
   c,aa,pp=angle;E=off+(fam-1)*39+34;F=off+(fam-1)*39+37
   T[E,18]=e0;T[E+1,18]=e0*c;T[E+1,14]=e0*aa;T[E+1,15]=e0*pp
   T[E+2,18]=e0*c*c/2;T[E+2,14]=e0*c*aa;T[E+2,15]=e0*c*pp
   T[E+2,16]=e0*aa*aa/2;T[E+2,17]=e0*pp*pp/2
   T[F,18]=f0
  end
 end
 T[331,239]=1;T[332,240]=1;T
end
function combined_experiment(a,p;D=8,M=32,rho=BigFloat(".04"))
 center=(a,p);a,p=big(a),big(p);inds=[(i,j) for i=0:D for j=0:D-i];L=length(inds)
 centers=[combined(complex(a),complex(p),center,k)[1] for k=1:2]
 patterns=[findall(!iszero,Matrix(c)) for c=centers]
 # Store Taylor coefficients only at structurally nonzero connection entries.
 Cs=[zeros(Complex{BigFloat},length(patterns[k]),L) for k=1:2]
 ls=zeros(Complex{BigFloat},3,332,L)
 lp=findall(!iszero,combined(complex(a),complex(p),center,1)[2])
 for v=0:M-1,w=0:M-1
  theta=2big(pi)*v/M;phi=2big(pi)*w/M
  aa=a+rho*cis(theta);pp=p+rho*cis(phi)
  c1,l=combined(aa,pp,center,1);c2,_=combined(aa,pp,center,2);cc=(c1,c2)
  for (kk,(i,j))=enumerate(inds)
   weight=cis(-i*theta-j*phi)/M^2
   for k=1:2,(e,ij)=enumerate(patterns[k]); Cs[k][e,kk]+=weight*rho*cc[k][ij]; end
   for ij=lp;ls[ij[1],ij[2],kk]+=weight*l[ij];end
  end
  if w==0 && v%5==0;println("sampling combined ",center," ",v,"/",M);flush(stdout);end
 end
 ix=Dict(v=>i for (i,v)=enumerate(inds));rows=Dict((0,0)=>ls)
 for order=0:D
  if order>0
   for i=0:order
    j=order-i;axis=i>0 ? 1 : 2;parent=i>0 ? (i-1,j) : (i,j-1)
    old=rows[parent];new=zeros(Complex{BigFloat},3,332,L)
    for (kk,(s,t))=enumerate(inds)
     s+t>D-order && continue
     nx=axis==1 ? (s+1,t) : (s,t+1)
     new[:,:,kk]=(axis==1 ? s+1 : t+1)*old[:,:,ix[nx]]
     for u=0:s,v=0:t
      ov=ix[(u,v)];cv=ix[(s-u,t-v)]
      for (e,ij)=enumerate(patterns[axis]),r=1:3
       new[r,ij[2],kk]+=old[r,ij[1],ov]*Cs[axis][e,cv]
      end
     end
    end
    rows[i,j]=new
   end
  end
  vals=[rows[(i,k-i)][:,:,1] for k=0:order for i=0:k]
  for (name,sel)=(("psi30",3:3),("1_chi_psi30",1:3))
   B=vcat([v[sel,:] for v=vals]...)
   for i=1:size(B,1);B[i,:]/=max(norm(B[i,:]),big"1e-100");end
   for (presentation,BB) in (("padded",B),("endpoint_quotient",B*physical_projection(a,p)))
    for i=1:size(BB,1);BB[i,:]/=max(norm(BB[i,:]),big"1e-100");end
    if order==D && name=="1_chi_psi30" && presentation=="endpoint_quotient"
     F=ComplexF64.(BB);Q=Vector{Vector{ComplexF64}}();selected=[]
     labels=[(src,i,k-i) for k=0:order for i=0:k for src=1:3]
     for r=1:size(F,1)
      v=copy(F[r,:]);for repeat=1:2,q in Q;v-=dot(q,v)*q;end
      if norm(v)>1e-8;push!(Q,v/norm(v));push!(selected,labels[r]);end
     end
     println((point=center,candidate_derivative_basis=selected,count=length(selected)));flush(stdout)
    end
    sv=svdvals(ComplexF64.(BB));println((point=center,order=order,output=name,presentation=presentation,ranks=[count(>(tol),sv) for tol in (1e-7,1e-9,1e-11)]));flush(stdout)
   end
  end
 end
end
if abspath(PROGRAM_FILE)==@__FILE__
 setprecision(160) do
  for point=((-.4,-1.3),(-.6,-1.1));combined_experiment(point...);end
 end
end
