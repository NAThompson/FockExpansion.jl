# Fock seeds at R* = 1 (the constants.jl runs, which were never logged), from cached exteriors.
include(joinpath(@__DIR__, "common_compute.jl"))

configs = [(12, "B12P6"), (12, "B10P6"), (14, "B10P6")]
names = ["a21 (physical)", "a_(0,1)", "a_(2,0)", "a_(0,2)", "a_(2,1)", "a_(0,3)", "a_(4,0)",
         "a_(2,2)", "a_(0,4)", "<psi40,Y(2,0)> (physical)", "<psi40,Y(0,2)> (physical)"]
vals = Dict{String,Any}("names" => names, "configs" => ["L=$L $n" for (L, n) in configs])
for (L, name) in configs
    ext = load_ext(L, 1.0, E_HE, name)
    row(p) = vcat(p.a21, p.a[2:9], p.psi40)
    t = @elapsed begin
        p44 = solve_point(44, L, 1.0, E_HE, ext)
        p48 = solve_point(48, L, 1.0, E_HE, ext)
    end
    vals["L=$L $name d=44"] = row(p44)
    vals["L=$L $name d=48"] = row(p48)
    @printf("L=%d %s: a21(44)=%.13f a21(48)=%.13f (%.0fs)\n", L, name, p44.a21, p48.a21, t)
    flush(stdout)
end
serialize(joinpath(DATA, "constants.jls"), vals)
println("wrote ", joinpath(DATA, "constants.jls"))
