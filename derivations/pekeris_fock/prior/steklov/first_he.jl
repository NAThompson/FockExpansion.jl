include("steklov.jl")
for d in (4, 8)
    model = channel_model(d)
    setup = interior_setup(model, E_HE, 1.0, d)
    for omega in (6, 8, 10, 12)
        t = @elapsed r = run_match(model, E_HE, 1.0, d; omega, setup)
        pc = physical_constants(r.setup, r.a)
        @printf("d=%d L=%d omega=%2d nterms=%d: lambda=%+.3e lambda2=%.3e a21=%.12f a=[%s] kept=%d/%d asym=%.1e (%.1fs)\n",
                d, d, omega, r.nterms, r.lambda, r.lambda2, pc.a21,
                join([@sprintf("%.6f", x) for x in r.a[2:min(4, end)]], ", "),
                r.ext.kept, r.ext.total, r.setup.asym, t)
    end
end
