"""Dependency-free exact algebra (Q coefficients); no CAS/repo imports."""
from fractions import Fraction as Q
from pathlib import Path
import json,resource
p=Path(__file__).resolve().parent
checks={}; facts={}
# Polynomial rational functions in t, stored low to high.
def clean(a):
 while len(a)>1 and not a[-1]:a.pop()
 return a
def add(a,b):return clean([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
def scale(a,c):return clean([x*c for x in a])
def mul(a,b):
 c=[Q(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return clean(c)
def der(a):return clean([i*a[i] for i in range(1,len(a))] or [Q(0)])
def radd(a,b):return (add(mul(a[0],b[1]),mul(b[0],a[1])),mul(a[1],b[1]))
def rscale(a,c):return (scale(a[0],c),a[1])
def rmul(a,b):return mul(a[0],b[0]),mul(a[1],b[1])
def rd(a):return (add(mul(der(a[0]),a[1]),scale(mul(a[0],der(a[1])),-1)),mul(a[1],a[1]))
one=[Q(1)];t=[Q(0),Q(1)];den=[Q(1),Q(0),Q(1)]
d4=mul(mul(den,den),mul(den,den))
sin8=(scale(mul(t,mul([Q(1),0,-1],[Q(1),0,-6,0,1])),8),d4)
f=radd(([Q(1,9)],one),rscale(sin8,Q(-7,45)))
f1=rd(f);f2=rd(f1);f3=rd(f2)
res=radd(radd(rmul((mul(den,den),one),f3),rmul((scale(mul(t,den),6),one),f2)),rmul(([Q(66),0,6],one),f1))
checks['guessed_ODE_annihilates_rational_fit']=all(x==0 for x in res[0])
facts['rational_fit']='1/9-(56/45)*t*(1-t^2)*(1-6*t^2+t^4)/(1+t^2)^4'
checks['fit_is_rational_so_order3_not_minimal']=len(f[0])>1 and len(f[1])>1
# 2-variable polynomials for the spherical Laplacian.
def pa(a,b):
 c=a.copy()
 for ij,x in b.items():c[ij]=c.get(ij,0)+x
 return {ij:x for ij,x in c.items() if x}
def ps(a,c):return {ij:c*x for ij,x in a.items() if c*x}
def pm(a,b):
 c={}
 for (i,j),x in a.items():
  for (k,l),y in b.items():c[(i+k,j+l)]=c.get((i+k,j+l),0)+x*y
 return c
def pd(a,idx):
 c={}
 for ij,x in a.items():
  if ij[idx]:
   kk=list(ij);kk[idx]-=1;c[tuple(kk)]=x*ij[idx]
 return c
u={(1,0):Q(1)};v={(0,1):Q(1)};I={(0,0):Q(1)};u2=pm(u,u);v2=pm(v,v);uv=pm(u,v)
H={(2,0):Q(2,9),(0,2):Q(2,9),(0,0):Q(-1,9),(1,1):Q(-14,45)}
Lap=pa(pa(pm(pa(I,ps(u2,-1)),pd(pd(H,0),0)),ps(pm(uv,pd(pd(H,0),1)),-2)),pm(pa(I,ps(v2,-1)),pd(pd(H,1),1)))
Lap=pa(Lap,pa(ps(pm(u,pd(H,0)),-3),ps(pm(v,pd(H,1)),-3)))
checks['harmonic_P40_H_zero']=not pa(Lap,ps(H,8))
fit={(0,0):Q(1,9),(1,1):Q(-14,45)}
checks['boundary_fit_minus_H_is_2delta_over9']=pa(fit,ps(H,-1))=={(0,0):Q(2,9),(2,0):Q(-2,9),(0,2):Q(-2,9)}
facts['harmonic_extension']='H=(2/9)*(u^2+v^2)-1/9-(14/45)*u*v = H40/27-7*H41/90+4*H42/27'
# SO5 type-II counterexample in basis e1,e2,f1,f2,h, <ei,fj>=deltaij.
Z=lambda:[[Q(0) for j in range(5)] for i in range(5)]
def mm(A,B):return [[sum(A[i][k]*B[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
def tr(A):return list(map(list,zip(*A)))
def ma(A,B):return [[A[i][j]+B[i][j] for j in range(5)] for i in range(5)]
G=Z();G[0][2]=G[2][0]=G[1][3]=G[3][1]=G[4][4]=Q(1)
def wedge(i,j):
 A=Z()
 for k in range(5):A[i][k]=G[j][k];A[j][k]=-G[i][k]
 return A
A=wedge(0,1);B=wedge(0,3)
checks['typeII_A_square_zero']=mm(A,A)==Z();checks['typeII_B_square_zero']=mm(B,B)==Z()
checks['typeII_A_orthogonal']=ma(mm(tr(A),G),mm(G,A))==Z()
checks['typeII_B_orthogonal']=ma(mm(tr(B),G),mm(G,B))==Z()
checks['distinct_images_double_nonzero']=mm(A,B)!=Z()
checks['distinct_images_triple_zero']=mm(B,mm(A,B))==Z()
checks['both_rank2_images_distinct']=set(i for i,row in enumerate(A) if any(row))=={0,1} and set(i for i,row in enumerate(B) if any(row))=={0,3}
facts['counterexample']={'A':[[str(x) for x in row] for row in A],'B':[[str(x) for x in row] for row in B]}
# Log recurrence signs: for l=log(-cos X), c*l'=-s and (c*l')'=-c.
# Check an explicit elementary manufactured case at X=pi, A=cos X:
# M(A*l)-l*M(A) = -2*sin X*A'-cos X*A, opposite code.
cc=([Q(1),0,-1],den);ss=(scale(t,2),den)
cr=rmul(rd(rscale(ss,-1)),(scale(den,Q(1,2)),one))
checks['positive_log_cross_zero_order_is_minus_cos']=all(x==0 for x in radd(cr,cc)[0])
facts['log_cross_term']='positive log: -2*sin(X)*A_prime-cos(X)*A; code uses the opposite, hence negative log'
# Exact indicial polynomials at the two internal divisors, from supplied exports.
import math
root=p.parents[1]
def shiftpoly(c,a):return [sum(c[i]*Q(math.comb(i,j))*a**(i-j) for i in range(j,len(c))) for j in range(len(c))]
def falling(n):
 a=[Q(1)]
 for j in range(n):a=mul(a,[Q(-j),Q(1)])
 return a
def indicial(rows,at):
 sh=[shiftpoly(row,at) for row in rows]
 vals=[next(i for i,c in enumerate(row) if c) for row in sh]
 k=min(v-j for j,v in enumerate(vals))
 a=[Q(0)]
 for j,(row,v) in enumerate(zip(sh,vals)):
  if v-j==k:a=add(a,scale(falling(j),row[v]))
 return scale(a,1/a[-1])
rowsA=[[Q(x) for x in row.split()] for row in (root/'claude10/op_e2_A.txt').read_text().splitlines()]
rowsC=[[Q(x) for x in row.split()] for row in (root/'claude10/op_eC.txt').read_text().splitlines()]
expA=[Q(0),Q(1,2),Q(1),Q(2),Q(3)]
expC=[Q(0),Q(0),Q(1),Q(1),Q(2)]
def product_roots(roots):
 a=[Q(1)]
 for x in roots:a=mul(a,[-x,Q(1)])
 return a
checks['u_minus1_exact_indicial']=indicial(rowsA,Q(-1))==product_roots(expA)
checks['v_plus1_exact_indicial']=indicial(rowsC,Q(1))==product_roots(expC)
facts['u_minus1_exponents']=[str(x) for x in expA]
facts['v_plus1_exponents']=[str(x) for x in expC]
assert all(checks.values()),checks
res={'checks':checks,'facts':facts,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(p/'x26_exact.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
