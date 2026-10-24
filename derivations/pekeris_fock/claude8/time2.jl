include("stencil2.jl")
blocks = read_terms("classical_terms.txt")
for k in (sort(collect(keys(blocks)))[1], "('sigma2P', 'g17', 'g17')")
    block_recurrence(blocks[k])
    t = @elapsed r = block_recurrence(blocks[k]); println(k, " offsets ", length(r), " time ", round(t; digits = 3))
end
