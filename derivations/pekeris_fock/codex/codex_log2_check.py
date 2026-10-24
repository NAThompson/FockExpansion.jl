import mpmath as m
m.mp.dps=90
B=[m.mpf(1)/6,-m.mpf(1)/30,m.mpf(1)/42,-m.mpf(1)/30,m.mpf(5)/66,-m.mpf(691)/2730,m.mpf(7)/6]
for s in range(2,6):
    M=m.mpf(1000)
    z=m.fsum(m.mpf(k)**(-s) for k in range(1,1000))+M**(1-s)/(s-1)+M**(-s)/2
    rising=m.mpf(s)
    for j,b in enumerate(B,1):
        z+=b/m.factorial(2*j)*rising*M**(-s-2*j+1)
        rising*=(s+2*j-1)*(s+2*j)
    print('zeta',s,'EM-minus-reference',m.nstr(z-m.zeta(s),20))
for x in (1,3,9):
    ps=[m.polygamma(j,x) for j in range(4)]
    d0,d1,d2,d3=ps
    bells=[1,d0,d0*d0+d1,d0**3+3*d0*d1+d2,d0**4+6*d0*d0*d1+4*d0*d2+3*d1*d1+d3]
    for p in range(5):
        expected=m.diff(m.gamma,m.mpf(x),p)
        err=abs(m.gamma(x)*bells[p]-expected)
        assert err< m.mpf('1e-75')
print('PASS Bell polynomials p=0..4 at x=1,3,9 against independent Gamma derivatives')
