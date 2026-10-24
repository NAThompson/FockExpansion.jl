"""Exact rational interpolation via polynomial extended Euclid, with holdouts.
Degrees are minimal-total-degree among regular rational interpolants to training
samples. Held-out agreement is evidence for, not proof of, the generic function.
"""
P=1004535809

def trim(a):
 while a and not a[-1]:a.pop()
 return a

def add(a,b,sgn=1):
 c=a[:]+[0]*max(0,len(b)-len(a))
 for i,v in enumerate(b):c[i]=(c[i]+sgn*v)%P
 return trim(c)
def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%P
 return trim(c)
def divmodp(a,b):
 a=a[:];q=[0]*max(0,len(a)-len(b)+1);inv=pow(b[-1],-1,P)
 while a and len(a)>=len(b):
  k=len(a)-len(b);v=a[-1]*inv%P;q[k]=v
  for j,x in enumerate(b):a[k+j]=(a[k+j]-v*x)%P
  trim(a)
 return trim(q),a
def ev(a,x):
 r=0
 for v in reversed(a):r=(r*x+v)%P
 return r
def interp(xs,ys):
 g=[1];f=[]
 for x,y in zip(xs,ys):
  z=(y-ev(f,x))*pow(ev(g,x),-1,P)%P
  f=add(f,[z*v%P for v in g]);g=mul(g,[-x%P,1])
 return g,f
def fit(xs,ys,holdx=(),holdy=()):
 if not any(ys):
  return dict(degrees=[-1,0],numerator=[],denominator=[1],heldout=not any(holdy))
 g,f=interp(xs,ys);r0,r1=g,f;b0,b1=[],[1];candidates=[]
 while r1:
  if all(ev(b1,x) for x in xs):
   candidates.append((len(r1)+len(b1)-2,r1[:],b1[:]))
  q,r=divmodp(r0,r1);b=add(b0,mul(q,b1),-1);r0,r1=r1,r;b0,b1=b1,b
 if not candidates:return None
 _,num,den=min(candidates,key=lambda v:(v[0],max(len(v[1]),len(v[2]))))
 scale=pow(den[-1],-1,P);num=[v*scale%P for v in num];den=[v*scale%P for v in den]
 good=all(ev(den,x) and ev(num,x)==y*ev(den,x)%P for x,y in zip(holdx,holdy))
 return dict(degrees=[len(num)-1,len(den)-1],numerator=num,denominator=den,heldout=good)

if __name__=='__main__':
 xs=list(range(1,16));ys=[(x*x+3)*pow(x+19,-1,P)%P for x in xs]
 r=fit(xs[:10],ys[:10],xs[10:],ys[10:]);assert r['degrees']==[2,1] and r['heldout']
 assert fit(xs,[0]*len(xs))['degrees']==[-1,0]
 print('PASS rational interpolation manufactured cases')
