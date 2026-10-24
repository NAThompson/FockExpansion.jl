include("stencil.jl")
blocks = read_terms("classical_terms.txt")
k = sort(collect(keys(blocks)))[1]
println(k, " terms ", length(blocks[k]))
@time st = stencil(blocks[k]); @time r = recurrence(st); println(length(r))
k2 = "('sigma2P', 'g17', 'g17')"; println(k2, " terms ", length(blocks[k2]))
@time r2 = recurrence(stencil(blocks[k2])); println(length(r2))
