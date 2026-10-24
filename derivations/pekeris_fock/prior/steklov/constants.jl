# Final constants at R* = 1 with d-tail extrapolation (increments ~ d^-5.87, error ~ d^-4.87).
include("production.jl")
p_err = 4.87
# Geometric-in-power-law tail: error(d) = C d^-p, so error(48) = (x48 - x44) / ((48/44)^p - 1).
tailpow(x44, x48) = (x48 - x44) / ((48 / 44)^p_err - 1)
names = ["a1", "a2(0,1)", "a3(2,0)", "a4(0,2)", "a5(2,1)", "a6(0,3)", "a7(4,0)", "a8(2,2)", "a9(0,4)"]
for (L, name) in ((12, "B12P6"), (12, "B10P6"), (14, "B10P6"))
    ext = get_ext(L, 1.0, E_HE, name)
    p44 = solve_point(44, L, 1.0, E_HE, ext)
    p48 = solve_point(48, L, 1.0, E_HE, ext)
    @printf("L=%d %s\n", L, name)
    @printf("  a21: d=48 %.12f  tail %+.2e  extrap %.12f\n", p48.a21, tailpow(p44.a21, p48.a21),
            p48.a21 + tailpow(p44.a21, p48.a21))
    for i = 2:9
        @printf("  %-8s d=48 %+.12f  extrap %+.12f\n", names[i], p48.a[i], p48.a[i] + tailpow(p44.a[i], p48.a[i]))
    end
    for (j, lab) in enumerate(("Y(2,0)", "Y(0,2)"))
        @printf("  <psi40,%s> d=48 %+.12f  extrap %+.12f\n", lab, p48.psi40[j],
                p48.psi40[j] + tailpow(p44.psi40[j], p48.psi40[j]))
    end
    flush(stdout)
end
