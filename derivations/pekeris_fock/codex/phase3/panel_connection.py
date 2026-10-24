"""Explicit outer-variable Pfaffian connection; no repository imports.
Taylor jets are ordinary coefficients (derivative / factorial).
Order: J,K,J_z[0:8],J_0,E,F0: 13 generating states, 65 padded jets.
"""
import mpmath as m

def add(a,b): return [a[i]+b[i] for i in range(3)]
def neg(a): return [-v for v in a]
def sub(a,b): return add(a,neg(b))
def mul(a,b): return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(3)]
def scale(a,c): return [c*v for v in a]
def const(c): return [c,m.mpf(0),m.mpf(0)]
def inv(a):
    v=[1/a[0],0,0]
    for n in (1,2): v[n]=-sum(a[j]*v[n-j] for j in range(1,n+1))/a[0]
    return v

def geometry(b,y):
    h=m.sqrt((1+b)/(1-b)); delta=m.sqrt((4*y*y-2)**2-4*(1-b*b)+0j)
    r1=m.sqrt((2-4*y*y+delta)/(2*(1-b)))
    r2=m.sqrt((2-4*y*y-delta)/(2*(1-b)))
    return h,[r1,-r1,r2,-r2,1j,-1j,1j*h,-1j*h],min(m.mpf(1),h)

def data(b,y,axis,family):
    h,z,u=geometry(b,y)
    hb=h/(1-b*b)
    dz=[]
    for q in z[:4]:
        tq=4*(1-b)*q**3+2*(4*y*y-2)*q
        dt=1-q**4 if axis==0 else 8*y*q*q
        dz.append(-dt/tq)
    dz += [0,0,1j*hb if axis==0 else 0,-1j*hb if axis==0 else 0]
    du=hb if axis==0 and b<0 else 0
    a=[const(-m.mpf('.5')) for _ in range(4)]
    tj=[0,1,0] if family=='A' else [0,0,0]
    sj=[0,1,0] if family=='D' else [0,0,0]
    a += [scale(tj,1/(2j)),scale(tj,-1/(2j)),scale(sj,1/(2j)),scale(sj,-1/(2j))]
    l0=-1/(2*(1+b)) if axis==0 else 0
    kap=const(l0)
    for ai,zi,vi in zip(a,z,dz): kap=sub(kap,scale(ai,vi/zi))
    return h,z,u,dz,du,a,l0,kap

def state(b,y,family):
    h,z,u=geometry(b,y)
    def W(q): return m.sqrt(1+b+(4*y*y-2)*q*q+(1-b)*q**4)
    f0=1/m.sqrt(1+b)
    result=[[0]*3 for _ in range(13)]
    for n in range(3):
        def f(q):
            angle=m.atan(q) if family=='A' else m.atan(q/h)
            return angle**n/(m.factorial(n)*W(q))
        result[0][n]=m.quad(f,[0,u])
        result[1][n]=m.quad(lambda q:q*f(q),[0,u])
        for k,zi in enumerate(z): result[2+k][n]=m.quad(lambda q:f(q)/(q-zi),[0,u])
        # Hadamard finite part, log(q) cutoff convention. All positive jets converge ordinarily.
        base=f0 if n==0 else 0
        result[10][n]=m.quad(lambda q:(f(q)-base)/q,[0,u])+base*m.log(u)
        result[11][n]=f(u)
        result[12][n]=f0 if n==0 else 0
    return result

def rhs(b,y,axis,family,X):
    h,z,u,dz,du,a,l0,kap=data(b,y,axis,family)
    J,K=X[:2]; Ji=X[2:10]; J0,E,F0=X[10:]
    out=[mul(kap,J),mul(kap,K)]
    for ai,zi,vi,ii in zip(a,z,dz,Ji):
        out[0]=sub(out[0],scale(mul(ai,ii),vi))
        out[1]=sub(out[1],scale(mul(ai,add(J,scale(ii,zi))),vi))
    out[0]=add(out[0],scale(E,du));out[1]=add(out[1],scale(E,u*du))
    for i in range(8):
        mixed=[None if i==j else scale(sub(Ji[i],Ji[j]),1/(z[i]-z[j])) for j in range(8)]
        boundary=add(scale(E,1/(u-z[i])),scale(F0,1/z[i]))
        for j in range(8):
            if j!=i: boundary=sub(boundary,mul(a[j],mixed[j]))
        hii=mul(boundary,inv(sub(a[i],const(1))))
        row=add(mul(kap,Ji[i]),scale(hii,dz[i]))
        for j in range(8):row=sub(row,scale(mul(a[j],hii if j==i else mixed[j]),dz[j]))
        out.append(add(row,scale(E,du/(u-z[i]))))
    row=mul(kap,J0)
    for ai,zi,vi,ii in zip(a,z,dz,Ji):row=add(row,scale(mul(ai,sub(J0,ii)),vi/zi))
    out.append(add(row,scale(E,du/u)))
    ep=kap
    for ai,zi,vi in zip(a,z,dz):ep=add(ep,scale(ai,(du-vi)/(u-zi)))
    out += [mul(ep,E),scale(F0,l0)]
    return out

def moments_from_state(b,y,XA,XD):
    h,z,u=geometry(b,y)
    def Q(X,family):
        _,_,_,_,_,a,_,_=data(b,y,0,family)
        out=scale(X[11],u**3)
        for ai,zi,ii in zip(a,z,X[2:10]):
            out=sub(out,mul(ai,add(scale(X[1],zi),add(scale(X[0],zi**2),scale(ii,zi**3)))))
        return out
    qa,qd=Q(XA,'A'),Q(XD,'D')
    ja,jd=XA[0],XD[0]
    N0=m.pi/2*ja[0]-2*ja[1]+4/m.pi*(ja[2]-jd[2])
    HD=(XD[8][1]+XD[9][1])/2
    N2=-m.pi/2*qa[0]+2*qa[1]-4/m.pi*(qa[2]-qd[2])
    N2+=4*h/m.pi*(XD[1][1]+XD[10][1]+4*(y*y-1)/(1-b)*HD)
    T=1+b+(4*y*y-2)*u*u+(1-b)*u**4
    N2-=4*T*XD[11][2]/(m.pi*(1-b)*u)
    return m.re(N0),m.re(N2)

def direct_moments(b,y):
    h,_,u=geometry(b,y)
    def f(q,n):
        A=m.pi/2-m.atan(q); D=m.atan(q/h)
        T=1+b+(4*y*y-2)*q*q+(1-b)*q**4
        return 2/m.pi*((A*A-D*D) if n==0 else ((h*D/q)**2-(q*A)**2))/m.sqrt(T)
    return tuple(m.quad(lambda q:f(q,n),[0,u]) for n in (0,2))
