from mpmath import mp, mpf, pi, pslq, catalan, log
mp.dps = 15
vals = {"l0 Z": -0.001355480179782, "l0 Z2": 0.003059320401767, "l0 Z3": 0.000701904941965,
        "l0 EZ": 0.000839853538016, "l0 aZ": -0.006040684371402,
        "l2 Z": 0.005105394795546, "l2 Z2": -0.005367956424392, "l2 Z3": -0.001900348776693}
for k, x in vals.items():
    x = mpf(x)
    for basis, nm in (([1, 1/pi, 1/pi**2], "1,1/π,1/π²"), ([1, 1/pi, 1/pi**2, catalan/pi, log(2), catalan/pi**2, log(2)/pi], "+G,log2")):
        r = pslq([x] + basis, maxcoeff=20000, maxsteps=100000, tol=mpf(10)**-12)
        if r: print(f"{k:7s} = {nm}: {r}"); break
    else: print(f"{k:7s}: no relation found")
