# Native styled output: Julia's IO color setting is respected, as is NO_COLOR.
ui_io()=IOContext(stdout,:color=>get(stdout,:color,false) && !haskey(ENV,"NO_COLOR"))
function ui_header()
 printstyled(ui_io(),"\nψ₃,₀ · recurrence certificate\n";bold=true,color=:cyan)
 println("Exact arithmetic · Julia / Nemo–FLINT\n")
end
function ui_stage(label)
 printstyled(ui_io(),label*"\n";bold=true);flush(stdout)
end
function ui_pass(label;seconds=nothing)
 printstyled(ui_io(),"  ✓ PASS  ";bold=true,color=:green)
 print(label)
 isnothing(seconds) || print("  ·  ",round(seconds;digits=2)," s")
 println();flush(stdout)
end
function ui_fail(label)
 printstyled(ui_io(),"\n  ✗ FAILED  ";bold=true,color=:red)
 println(label);flush(stdout)
end
