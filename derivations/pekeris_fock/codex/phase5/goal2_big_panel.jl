# Standalone port of panel_connection.py; standard library only.
using LinearAlgebra, SparseArrays
BLAS.set_num_threads(1)
const nstate=78
const basepoint=Ref((.3,-.6))
con(a,b) = [sum(a[j+1]*b[k-j+1] for j=0:k) for k=0:2]
function jm(a,X)
 T=spzeros(Complex{BigFloat},3,3)
 for k=0:2,j=0:k; T[k+1,k-j+1]=a[j+1]; end
 T*X
end
function ji(a)
 v=Complex{BigFloat}[1/a[1],0,0]; for k=1:2; v[k+1]=-sum(a[j+1]*v[k-j+1] for j=1:k)/a[1]; end; v
end
jc(a)=Complex{BigFloat}[a,0,0]
function geom(b,y,side)
 h=sqrt((1+b)/(1-b)); b0,y0=big.(basepoint[])
 dc=complex((4y0*y0-2)^2-4*(1-b0*b0)); de0=sqrt(dc)
 de=de0*sqrt(((4y*y-2)^2-4*(1-b*b))/dc);
 v1=(2-4y0*y0+de0)/(2*(1-b0)); v2=(2-4y0*y0-de0)/(2*(1-b0))
 z1=sqrt(v1)*sqrt(((2-4y*y+de)/(2*(1-b)))/v1)
 z2=sqrt(v2)*sqrt(((2-4y*y-de)/(2*(1-b)))/v2)
 h,Complex{BigFloat}[z1,-z1,z2,-z2,im,-im,im*h,-im*h],side<0 ? h : one(h)
end
function info(b,y,side,axis,fam)
 h,z,u=geom(b,y,side); hb=h/(1-b*b)
 dz=Complex{BigFloat}[-(axis==1 ? 1-q^4 : 8y*q*q)/(4*(1-b)*q^3+2*(4y*y-2)*q) for q in z[1:4]]
 append!(dz,[0,0,axis==1 ? im*hb : 0,axis==1 ? -im*hb : 0])
 du=axis==1 && side<0 ? hb : 0
 a=[jc(-.5) for _=1:4]; tj=fam==1 ? Complex{BigFloat}[0,1,0] : jc(0); sj=fam==2 ? Complex{BigFloat}[0,1,0] : jc(0)
 append!(a,[tj/(2im),-tj/(2im),sj/(2im),-sj/(2im)])
 l0=axis==1 ? -1/(2*(1+b)) : 0; kap=jc(l0)
 for k=1:8; kap-=a[k]*(dz[k]/z[k]); end
 h,z,u,dz,du,a,l0,kap
end
function connection(b,y,side,axis)
 out=spzeros(Complex{BigFloat},78,78)
 for fam=1:2
  X=[spzeros(Complex{BigFloat},3,78) for _=1:13]
  for i=1:13,k=1:3; X[i][k,(fam-1)*39+(i-1)*3+k]=1; end
  h,z,u,dz,du,a,l0,kap=info(b,y,side,axis,fam)
  J,K=X[1:2]; Is=X[3:10]; J0,E,F0=X[11:13]
  Y=[jm(kap,J),jm(kap,K)]
  for k=1:8
   Y[1]-=dz[k]*jm(a[k],Is[k]); Y[2]-=dz[k]*jm(a[k],J+z[k]*Is[k])
  end
  Y[1]+=du*E; Y[2]+=u*du*E
  for i=1:8
   mixed=[j==i ? zero(E) : (Is[i]-Is[j])/(z[i]-z[j]) for j=1:8]
   bd=E/(u-z[i])+F0/z[i]
   for j=1:8; j!=i && (bd-=jm(a[j],mixed[j])); end
   hii=jm(ji(a[i]-jc(1)),bd); row=jm(kap,Is[i])+dz[i]*hii
   for j=1:8; row-=dz[j]*jm(a[j],j==i ? hii : mixed[j]); end
   push!(Y,row+du/(u-z[i])*E)
  end
  row=jm(kap,J0)
  for i=1:8; row+=dz[i]/z[i]*jm(a[i],J0-Is[i]); end
  push!(Y,row+du/u*E); ep=copy(kap)
  for i=1:8; ep+=a[i]*(du-dz[i])/(u-z[i]); end
  append!(Y,[jm(ep,E),l0*F0])
  for i=1:13,k=1:3; out[(fam-1)*39+(i-1)*3+k,:]=Y[i][k,:]; end
 end
 out
end
function outputs(b,y,side)
 h,z,u=geom(b,y,side); X=[[spzeros(Complex{BigFloat},3,78) for _=1:13] for _=1:2]
 for f=1:2,i=1:13,k=1:3; X[f][i][k,(f-1)*39+(i-1)*3+k]=1; end
 qs=[]
 for f=1:2
  a=info(b,y,side,1,f)[6]; V=X[f]; q=u^3*V[12]
  for i=1:8; q-=jm(a[i],z[i]*V[2]+z[i]^2*V[1]+z[i]^3*V[i+2]); end
  push!(qs,q)
 end
 A,D=X; qa,qd=qs
 l0=big(pi)/2*A[1][1,:]-2*A[1][2,:]+4/big(pi)*(A[1][3,:]-D[1][3,:])
 hd=(D[9][2,:]+D[10][2,:])/2
 l2=-big(pi)/2*qa[1,:]+2*qa[2,:]-4/big(pi)*(qa[3,:]-qd[3,:])
 l2+=4h/big(pi)*(D[2][2,:]+D[11][2,:]+4*(y*y-1)/(1-b)*hd)
 T=1+b+(4y*y-2)*u*u+(1-b)*u^4
 l2-=4T/(pi*(1-b)*u)*D[12][3,:]
 [transpose(l0); transpose(l2)]
end
function experiment(b,y;D=8,M=24,rho=.035)
 basepoint[]=(b,y); b,y=big(b),big(y); side=sign(b); inds=[(i,j) for i=0:D for j=0:D-i]; L=length(inds)
 Cs=[zeros(Complex{BigFloat},78,78,L) for _=1:2]; ls=zeros(Complex{BigFloat},2,78,L)
 patterns=[findall(!iszero,Matrix(connection(complex(b),complex(y),side,k))) for k=1:2]
 lp=findall(!iszero,Matrix(outputs(complex(b),complex(y),side)))
 # Normalized Taylor variables (b-b0)/rho,(y-y0)/rho.
 for p=0:M-1,q=0:M-1
  t=2big(pi)*p/M; s=2big(pi)*q/M; bb=b+rho*cis(t); yy=y+rho*cis(s)
  C=[rho*connection(bb,yy,side,k) for k=1:2]; l=outputs(bb,yy,side)
  if q==0 && p%5==0; println("sampling ",b," ",p,"/",M);flush(stdout);end
  for (v,(i,j))=enumerate(inds)
   w=cis(-i*t-j*s)/M^2
   for k=1:2, ij=patterns[k]; Cs[k][ij[1],ij[2],v]+=w*C[k][ij]; end
   for ij=lp; ls[ij[1],ij[2],v]+=w*l[ij]; end
  end
 end
 idx=Dict(ij=>k for (k,ij)=enumerate(inds))
 rows=Dict((0,0)=>ls)
 for order=0:D
  if order>0
   for a=0:order
    bidx=order-a; axis=a>0 ? 1 : 2; parent=a>0 ? (a-1,bidx) : (a,bidx-1)
    old=rows[parent]; new=zeros(Complex{BigFloat},2,78,L)
    for (v,(i,j))=enumerate(inds)
     i+j>D-order && continue
     nxt=axis==1 ? (i+1,j) : (i,j+1)
     new[:,:,v]=(axis==1 ? i+1 : j+1)*old[:,:,idx[nxt]]
     for p=0:i,q=0:j
      ov=idx[(p,q)]; cv=idx[(i-p,j-q)]
      for ij=patterns[axis], r=1:2
       new[r,ij[2],v]+=old[r,ij[1],ov]*Cs[axis][ij[1],ij[2],cv]
      end
     end
    end
    rows[(a,bidx)]=new
   end
  end
  vals=[rows[(a,k-a)][:,:,1] for k=0:order for a=0:k]
  for (name,sel)=(("N0",1:1),("N2",2:2),("joint",1:2))
   B=vcat([v[sel,:] for v=vals]...)
   for i=1:size(B,1); B[i,:]/=max(norm(B[i,:]),1e-300); end
   sv=svdvals(ComplexF64.(B))
   println((point=(b,y),D=D,M=M,rho=rho,order=order,output=name,ranks=[count(>(tol),sv) for tol in (1e-7,1e-9,1e-11)],tail=sv[max(1,end-2):end]))
  end
 end
end
if abspath(PROGRAM_FILE)==@__FILE__
setprecision(160) do
 for point=((.3,-.6),(-.27,-.7)); experiment(point...;D=8,M=40,rho=BigFloat(".06")); end
end
end
