# Cached exterior Steklov forms and interior sweeps.
include("steklov.jl")
include("exterior_dd.jl")
using Serialization

const KAPPA = sqrt(-2(E_HE + 2))

# Exterior bases: Hylleraas sets (za, zb, omega) plus P channel corrections per trace harmonic.
function basis_spec(name)
    specs = Dict(
        "B10P4" => ([(2.0, KAPPA, 10), (4.0, 4.0, 6), (8.0, 8.0, 6)], 4),
        "B12P4" => ([(2.0, KAPPA, 12), (4.0, 4.0, 8), (8.0, 8.0, 6)], 4),
        "B10P6" => ([(2.0, KAPPA, 10), (4.0, 4.0, 6), (8.0, 8.0, 6)], 6),
        "B12P6" => ([(2.0, KAPPA, 12), (4.0, 4.0, 8), (8.0, 8.0, 6)], 6),
        "B12P6x" => ([(2.0, KAPPA, 12), (1.2, 1.2, 6), (4.0, 4.0, 8), (8.0, 8.0, 6), (16.0, 16.0, 4)], 6),
        "B14P8" => ([(2.0, KAPPA, 14), (4.0, 4.0, 8), (8.0, 8.0, 8)], 8),
        "B8P2" => ([(2.0, KAPPA, 8), (4.0, 4.0, 4), (8.0, 8.0, 4)], 2),
    )
    return specs[name]
end

# quad = (nphi, nu); the default (24, 20) cache files carry no tag.
function get_ext(L, Rs, E, name; verbose = true, quad = (24, 20))
    tag = quad == (24, 20) ? "" : "_q$(quad[1])x$(quad[2])"
    file = "cache/ext_L$(L)_R$(Rs)_E$(E)_$(name)$(tag).jls"
    isfile(file) && return deserialize(file)
    zsets, nchan = basis_spec(name)
    model0 = channel_model(L)
    nL = length(model0.basis.K)
    t = @elapsed ext = exterior_form_dd(model0, E, Rs, nL; zsets, nchan, lchan = L, verbose,
                                        nphi = quad[1], nu = quad[2])
    verbose && @printf("  exterior L=%d R*=%.2f %s: %.1fs\n", L, Rs, name, t)
    serialize(file, ext)
    return ext
end

function solve_point(d, L, Rs, E, ext; N = 2d + (Rs >= 1.5 ? 44 : 28))
    model = channel_model(d)
    setup = interior_setup(model, E, Rs, L; N = max(44, N))
    r = match(setup, Float64.(ext.S))
    pc = physical_constants(setup, r.a)
    return (; d, L, Rs, E, lambda = r.lambda, lambda2 = r.lambda2, a21 = pc.a21, a = r.a,
            psi40 = pc.psi40, asym = setup.asym, ode = setup.ode, labels = model.basis.labels)
end
