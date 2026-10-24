include("goal2_panel_rank.jl")
# Coordinates a=alpha-pi/2, p=p0-pi. Trigonometric coefficients are
# algebraic over physical distance coordinates on the chosen chart.
function classicalC(a,p,side,axis)
 C=zeros(ComplexF64,18,18); pi0=pi
 args=Tuple[]
 for s=(1,-1)
  w=(a+s*p)/4
  for m=(-1,0,1); push!(args,(w+m*pi/4,1/4,s/4,0,m*pi/4)); end
  push!(args,(w,1/4,s/4,1,0))
 end
 append!(args,[(p/2,0,1/2,0,0),((a+pi/2)/2,1/2,0,0,pi/4),((pi/2-a)/2,-1/2,0,0,pi/4)])
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
 Z=2; al=a+pi/2; ps=p+pi
 r1=cos(al/2);r2=sin(al/2);sg=r1+r2;d=r1-r2;xi=sqrt(2)*sin(ps/2);eta0=sqrt(2)*cos(ps/2)
 out=zeros(ComplexF64,1,18); idx=0
 for s=(1,-1)
  eta=s*eta0
  for m=(-1,0,1,2)
   dd,ss=d,sg
   for _=1:mod(m,4); dd,ss=-ss,dd; end
   kz=(dd*(-9ss^2+6ss*xi+8xi^2-8)-eta*(5-2xi^2))/(144pi)
   kz2=-5*(xi^2-1)*(dd+eta)/(32pi)
   nu=m==0 || m==2 ? xi^2-sg^2 : (9sg^2+xi^2-2+(m==1 ? -8 : 8)*sg*eta)
   nv=d*nu/(36pi)
   idx+=1; out[idx]=(m==2 ? -4 : 4)*(Z*kz+Z^2*(kz2+nv))
  end
  el=Z*(7sg*d*eta+9sg^3+8sg*xi^2-26sg+10xi^3-12xi)/72+Z^2*sg*(137d*eta+64sg^2-154xi^2+26)/288
  out[s==1 ? 12 : 13]=el
 end
 c9=-eta0*(Z*(2xi^2-5)/(9pi)-5Z^2*(xi^2-1)/(2pi));out[9]=c9/2
 out[10]=2Z^2/(9pi)*r2*(8r2^2+xi^2-5)
 out[11]=2Z^2/(9pi)*r1*(8r1^2+xi^2-5)
 out[16]=Z*xi*(2xi^2+1)/(144pi)-5Z^2*xi*(xi^2-1)/(32pi)
 out[14]=7Z*xi*d*sg/144+Z^2*(xi*d*sg/64+d*(xi^2-1)/(9pi))
 out[17]=Z^2*sg*(4sg^2+xi^2-9)/(36pi)
 out[15]=-Z*eta0*(5xi^2-4)/(36pi)+Z^2*sg*xi*eta0/(18pi)
 out
end

function classical_experiment(b,y;D=8,M=24,rho=.035)
 basepoint[]=(b,y); side=sign(b); inds=[(i,j) for i=0:D for j=0:D-i]; L=length(inds)
 Cs=[zeros(ComplexF64,18,18,L) for _=1:2]; ls=zeros(ComplexF64,1,18,L)
 # Normalized Taylor variables (b-b0)/rho,(y-y0)/rho.
 for p=0:M-1,q=0:M-1
  t=2pi*p/M; s=2pi*q/M; bb=b+rho*cis(t); yy=y+rho*cis(s)
  C=[rho*classicalC(bb,yy,side,k) for k=1:2]; l=classicalL(bb,yy,side)
  for (v,(i,j))=enumerate(inds)
   w=cis(-i*t-j*s)/M^2
   for k=1:2; Cs[k][:,:,v] .+= w.*C[k]; end
   ls[:,:,v] .+= w.*l
  end
 end
 idx=Dict(ij=>k for (k,ij)=enumerate(inds))
 rows=Dict((0,0)=>ls)
 for order=0:D
  if order>0
   for a=0:order
    bidx=order-a; axis=a>0 ? 1 : 2; parent=a>0 ? (a-1,bidx) : (a,bidx-1)
    old=rows[parent]; new=zeros(ComplexF64,1,18,L)
    for (v,(i,j))=enumerate(inds)
     i+j>D-order && continue
     nxt=axis==1 ? (i+1,j) : (i,j+1)
     new[:,:,v]=(axis==1 ? i+1 : j+1)*old[:,:,idx[nxt]]
     for p=0:i,q=0:j
      new[:,:,v]+=old[:,:,idx[(p,q)]]*Cs[axis][:,:,idx[(i-p,j-q)]]
     end
    end
    rows[(a,bidx)]=new
   end
  end
  vals=[rows[(a,k-a)][:,:,1] for k=0:order for a=0:k]
  for (name,sel)=(("classical",1:1),)
   B=vcat([v[sel,:] for v=vals]...)
   for i=1:size(B,1); B[i,:]/=max(norm(B[i,:]),1e-300); end
   sv=svdvals(B)
   println((point=(b,y),D=D,M=M,rho=rho,order=order,output=name,ranks=[count(>(tol),sv) for tol in (1e-7,1e-9,1e-11)],tail=sv[max(1,end-2):end]))
  end
 end
end

if abspath(PROGRAM_FILE)==@__FILE__
for point=((-.4,-1.3),(-.7,-.9)); classical_experiment(point...;D=8,M=40,rho=.15); end
end
