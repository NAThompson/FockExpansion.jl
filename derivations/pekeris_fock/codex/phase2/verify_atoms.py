#!/usr/bin/env python3
import sys
sys.dont_write_bytecode=True
import itertools,json,random
from pathlib import Path
import mpmath as mp
import sympy as S
import order2 as O
mp.mp.dps=45
OUT=Path(__file__).resolve().parent

def cl(t):return mp.im(mp.polylog(2,mp.exp(1j*t)))

def atoms(u,v,w):
    a,b,d=2*v+w,2*u+w,2*(u+v)
    q=a*a+b*b; p=2*q-d*d
    sn=2*a*b/q; cs=(a*a-b*b)/q; vv=(q-d*d)/q
    gg=d*mp.sqrt(p)/q
    x=2*mp.atan(b/a); beta=mp.asin(vv)
    L=mp.log((a+b+d)/mp.sqrt(q))
    A=mp.log(abs((gg-cs)/(gg+cs)))
    B=mp.log(abs((sn+vv)/(sn-vv)))
    J=mp.log((cs+(a+b)*d/q)**2/((1+sn)*(gg+cs)))
    dl=cl(x-beta)+cl(mp.pi+beta-x)-cl(mp.pi-x-beta)-cl(x+beta)
    U=beta*A+x*B-dl+mp.pi*J
    return [mp.mpf(1),L,x,mp.sqrt(p)/d*(beta+mp.pi/2),U,mp.log((a+b)/mp.sqrt(q))]

def chi_source(u,v,w):
    a,b,d=2*v+w,2*u+w,2*(u+v)
    q=a*a+b*b; x=2*mp.atan(b/a)
    ss=2*a*b/q; cc=(a*a-b*b)/q; vv=(q-d*d)/q
    mu=vv/ss; gg=mp.sqrt(1-vv*vv); sig=(a+b)/mp.sqrt(q); xi=d/mp.sqrt(q)
    beta=mp.asin(vv)
    dl=cl(x-beta)+cl(mp.pi+beta-x)-cl(mp.pi-x-beta)-cl(x+beta)
    logs=beta*(2*mp.log(ss)-2*mp.log(gg+cc))
    logs+=(beta+x)*mp.log1p(mu)+(beta-x)*mp.log1p(-mu)
    h0=-(3*mp.pi+10-16*mp.catalan)/(24*mp.pi)
    return (gg+ss-2*sig*xi)/6+(-2*mp.pi*vv*mp.log(sig+xi)
        +mp.pi*cc*mp.log((cc+sig*xi)**2/(sig**2*(gg+cc)))
        +2*gg*beta+cc*(logs-dl))/(6*mp.pi)+h0*vv

def matrix_from_json(path):
    raw=json.loads(path.read_text()); out=[]
    for var in O.xyz:
        mat=S.zeros(O.rank)
        for key,value in raw[str(var)].items():
            i,j=key.strip('()').split(',')
            mat[O.names.index(i),O.names.index(j)]=S.sympify(value,locals=vars(O.old))
        out.append(mat)
    return out

def main():
    C=matrix_from_json(OUT/'connection.json')
    functions=[S.lambdify(O.xyz,c,'mpmath') for c in C]
    errors=[]
    for point in [(mp.mpf(1)/3,mp.mpf(2)/3,mp.mpf(4)/5),
                  (mp.mpf(2)/5,mp.mpf(7)/10,mp.mpf(3)/2)]:
        values=atoms(*point)
        a,b,d,q,p,sn,cs,vv,t,f,den,A=O.geometry()
        h0=-(3*mp.pi+10-16*mp.catalan)/(24*mp.pi)
        geo=S.lambdify(O.xyz,(sn,cs,vv,t,(a+b)*d/q),'mpmath')(*point)
        ss,cc,V,T,Sd=geo
        reconstructed=(ss-2*Sd)/6-V*values[1]/3+(cc*values[4]+2*T*values[3])/(6*mp.pi)+h0*V
        errors.append(abs(reconstructed-chi_source(*point)))
        for i,fun in enumerate(functions):
            cm=fun(*point)
            for j in range(O.rank):
                e=tuple(int(k==i) for k in range(3))
                actual=mp.diff(lambda *z:atoms(*z)[j],point,e)
                expected=sum(values[k]*cm[k,j] for k in range(O.rank))
                errors.append(abs(actual-expected))
    print('Two source-chi comparisons and 36 direct first derivatives; max error',mp.nstr(max(errors),8))
    assert max(errors)<mp.mpf('1e-38')
    # Cyclic rank: work in distances, so no special perimetric relation obscures rank.
    a,b,d=S.symbols('a b d',positive=True)
    coords=(a,b,d); CD=O.connection(coords)
    _,_,_,q,p,sn,cs,vv,t,*_=O.geometry(coords)
    # chi modulo a rational function and a nonzero constant factor.
    z=S.Matrix([0,-2*S.pi*vv,0,2*t,cs,0])
    jet=lambda z,i:(z.diff(coords[i])+CD[i]*z).applyfunc(S.cancel)
    vectors=[S.Matrix([1,0,0,0,0,0]),z,jet(z,0),jet(z,1)]
    vectors.extend([jet(jet(z,0),0),jet(jet(z,0),1),jet(jet(z,1),1)])
    vectors.append(S.Matrix([0,0,0,0,0,1]))
    point={a:3,b:2,d:2}
    chosen=None
    for ids in itertools.combinations(range(len(vectors)),6):
        determinant=S.factor(S.Matrix.hstack(*(vectors[k].subs(point) for k in ids)).det())
        if determinant!=0:
            chosen=(ids,determinant);break
    assert chosen
    print('Formal cyclic module rank witness at (a,b,d)=(3,2,2):',chosen)
    (OUT/'cyclic_rank.json').write_text(json.dumps({'vectors':['1','z','da z','db z','daa z','dab z','dbb z','Slog'],
        'chosen':chosen[0],'determinant':str(chosen[1]),
        'qualification':'Rank in the presented differential module; actual-function independence is proved separately in minimal_rank.md.'},indent=2)+'\n')

if __name__=='__main__':main()
