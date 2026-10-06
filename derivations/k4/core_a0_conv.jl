include("core_a0.jl")
for (n, lv) in ((6, 8), (8, 10), (10, 12))
    t=@elapsed v=core_a0(0.7, 1.1, 2.0; nt=n, nu=n, levels=lv)
    @printf("graded n=%d levels=%d: % .13f  (Green 0.064502533732, diff %.1e, %.1f s)\n", n, lv, v, v-0.064502533732, t)
end
