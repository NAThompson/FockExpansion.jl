import sympy as s
r1,r2,c,Z=s.symbols('r1 r2 c Z');Q=r1*r1+r2*r2;S=r1+r2
PE=Z*S*(2*Q+r1*r2)/18-c*(6*Q-c*c)/72
PA=-Z*S*(Q-c*c)/2+c*(6*Q-5*c*c)/12
for name,p in [('E',PE),('a21',PA)]:
 p=s.Poly(s.expand(p),r1,r2,c)
 assert {sum(m) for m,k in p.terms()}=={3}
 assert s.expand(p.as_expr()-p.as_expr().xreplace({r1:r2,r2:r1}))==0
 print(name,'homogeneous cubic, exchange even:',p.as_expr())
print('In build3, a psi30 term has n>=3, a+b<=n-3, n<=N-2.')
print('Multiplying either cubic by s^(n-3)*x^a*y^b gives plain-space terms s^n*x^j*y^k with even j and j+k<=n; hence parameter changes preserve the finite span exactly.')
