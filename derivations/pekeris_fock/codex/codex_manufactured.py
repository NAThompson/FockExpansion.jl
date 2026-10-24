"""Independent quadrature and Mellin checks of the finite-part normalization."""
import mpmath as m
m.mp.dps=60
h=lambda r:m.exp(-r)*((r*r-9*r)*m.log(r)+8-2*r)
print('manufactured f=R^2 log(R) exp(-R), b=1, a=0')
for d in ('1e-4','1e-6','1e-8'):
    d=m.mpf(d)
    fp=m.quad(lambda r:h(r)/r,[d,m.mpf('.01'),1,m.inf])+8*m.log(d)
    print('cutoff',m.nstr(d,4),'FP',m.nstr(fp,24),'a',m.nstr(-fp/8-m.mpf(1)/8,20))
    assert abs(fp+1)<200*d*abs(m.log(d))
for z in ('1e-4','1e-6','1e-8'):
    z=m.mpf(z)
    fp=m.diff(m.gamma,z+2)-9*m.diff(m.gamma,z+1)+8*m.gamma(z)-2*m.gamma(z+1)-8/z
    print('Mellin z',m.nstr(z,4),'FP',m.nstr(fp,24))
    assert abs(fp+1)<20*z
print('PASS: both limits approach -1, hence a=0; omitting -b/8 gives 1/8.')
