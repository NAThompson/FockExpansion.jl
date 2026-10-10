#!/usr/bin/env -S julia --startup-file=no
# A standalone exact certificate environment; no Python runtime is used.
using Pkg
Pkg.activate(joinpath(@__DIR__, "julia"))
Pkg.instantiate()
using SHA, JSON3
include(joinpath(@__DIR__, "julia", "terminal_ui.jl"))
ui_header()
started=time()
try
    include(joinpath(@__DIR__, "julia", "verify.jl"))
catch
    ui_fail("Certificate did not complete; diagnostic follows.")
    rethrow()
end
files=("verify_all.jl", "julia/exact_engine.jl", "julia/terminal_ui.jl", "julia/verify.jl", "julia/formula_inputs.json", "julia/Project.toml", "julia/Manifest.toml")
report=Dict("status"=>"passed", "julia_version"=>string(VERSION),
    "seconds"=>time()-started, "sectors"=>["Z0","Z1","Z2","Z3","E","EZ","a21","a21Z"],
    "classical_table_links"=>18, "sha256"=>Dict(f=>bytes2hex(sha256(read(joinpath(@__DIR__,f)))) for f in files))
open(joinpath(@__DIR__, "certificate_julia.json"), "w") do io
    JSON3.pretty(io,report)
end
println()
ui_pass("Certificate verified · 8/8 sectors";seconds=time()-started)
println("Report: certificate_julia.json")
