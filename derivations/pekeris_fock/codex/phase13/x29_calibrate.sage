from pathlib import Path
import time, resource, json
load('codex/phase13/x29_ftest.sage')
start=time.time(); out=Path('codex/phase13')
# p is the x25 period: K((u-1)/(u+1))/sqrt(-u*(u+1)) on -1<u<0.
# g=(1+u)p+u*p'.  t=2log(sqrt(1-u^2))+log(-u)
#                         +log(sqrt(1+u))^2+Li_2(-u).
G=primitive(E.annihilator_of_associate((1+u)+u*D))
a=-2*u/(1-u**2)+1/u; b=1/(2*(1+u))-1/u
# t'=a+b*l, l=log(1+u), l'=1/(1+u).
M=matrix(K,[[0,0,0],[1/(1+u),0,0],[a,b,0]])
v=vector(K,[0,0,1]); rows=[v]
for j in range(3): rows.append(rows[-1].derivative(u)+rows[-1]*M)
C=matrix(K,rows[:3]); rel=rows[3]*C.inverse()
T=primitive(D**3-sum(rel[j]*D**j for j in range(3)))
L=primitive(G.lclm(T))
# Exact LCLM split; Ore quotient notation quo_rem.
Q,rem=L.quo_rem(G); assert rem==0
print('G=',G,flush=True); print('T=',T,flush=True); print('L order/degree=',L.order(),max(R(c).degree() for c in L),flush=True)
print('G gauge E:',gauge(G,E),flush=True)
# Construct complete first-order factorization of Q from the rational Tate connection:
# Q is equivalent to T because the elliptic/Tate modules have no common simple factors.
print('factor Q start',flush=True)
fq=Q.factor(); print('factor Q=',fq,flush=True)
res=ftest(L,fq+[G]); print('ftest:',res,flush=True)
assert any(h[1]=='E' for h in res['hits'])
assert not gauge(Et,E); assert gauge(Et,Et) is not None
assert not ftest(T,T.factor())['hits']
# exact polynomial ansatz rank determines ALL (r,d) with a formal connection identity.
# Connection direct sum (p,p',1,l,t), observed output g+t.
B=block_diagonal_matrix(companion(E),M)
w=vector(K,[1+u,u,0,0,1]); jets=[w]
for j in range(8): jets.append(jets[-1].derivative(u)+jets[-1]*B)
# For order 5 the unique primitive minimal operator gives minimal degree directly.
assert matrix(K,jets[:5]).rank()==5
assert sum(L[j]*jets[j] for j in range(6))==0
cs=[[str(R(c)[j]) for j in range(R(c).degree()+1)] for c in L]
payload={'order':int(L.order()),'degree':int(max(R(c).degree() for c in L)), 'coefficients':cs,
 'G':str(G),'T':str(T),'L':str(L),'Q':str(Q),'factors':[str(f) for f in fq+[G]],
 'hits':[(int(i),name,str(H)) for i,name,H in res['hits']],
 'jets':[[[[str(z) for z in q.numerator().list()],[str(z) for z in q.denominator().list()]] for q in row] for row in jets], 'exact_cyclic_rank':5,'seconds':time.time()-start,'rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(out/'x29_exact.json').write_text(json.dumps(payload,indent=2,default=int)+'\n')
print('EXACT CHECKS PASS',payload['seconds'],payload['rss_bytes'],flush=True)
