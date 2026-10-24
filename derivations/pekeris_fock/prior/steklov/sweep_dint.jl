# Exterior Steklov form once (it does not depend on d_int), then sweep d_int.
include("steklov.jl")
include("exterior_dd.jl")
using Serialization
L = parse(Int, ARGS[1])
Rs = parse(Float64, ARGS[2])
k = sqrt(-2(E_HE + 2))
zsets = [(2.0, k, 10), (4.0, 4.0, 6), (8.0, 8.0, 6)]
nchan, lchan = 4, L
file = "cache/ext_L$(L)_R$(Rs)_H1c$(lchan)P$(nchan).jls"
if isfile(file)
    ext = deserialize(file)
else
    model0 = channel_model(L)
    nL = length(model0.basis.K)
    ext = exterior_form_dd(model0, E_HE, Rs, nL; zsets, nchan, lchan, verbose = true)
    serialize(file, ext)
end
for d in parse.(Int, ARGS[3:end])
    model = channel_model(d)
    t = @elapsed setup = interior_setup(model, E_HE, Rs, L)
    r = match(setup, Float64.(ext.S))
    pc = physical_constants(setup, r.a)
    @printf("L=%d R*=%.2f d=%2d m=%3d: lambda=%+.4e a21=%.12f a3=%+.10f a4=%+.10f asym=%.1e ode=%.1e (%.1fs)\n",
            L, Rs, d, nchannels(model), r.lambda, pc.a21, r.a[3], r.a[4], setup.asym, setup.ode, t)
    flush(stdout)
end
