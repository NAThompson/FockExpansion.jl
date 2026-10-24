include("exterior_lib.jl")
using Printf
for d in (2, 4, 6)
    b = cuspmodel(d); E = cuspsolve(b).E
    Ys = riccati_coeffs(b, E, 45; prec = 256)
    lines = pole_scan(b, E, Ys; rmin = 0.4, dy = 0.01, ymax = 6.0)
    # local maxima in Im of max|Y| with the maximum attained off the real-axis endpoint
    for i = 2:length(lines)-1
        y, mx, at = lines[i]
        if mx > lines[i-1][2] && mx > lines[i+1][2] && at > 0.41
            @printf("d=%d candidate singularity near R = %.3f + %.2fi  (max|Y|=%.1e)\n", d, at, y, mx)
        end
    end
    flush(stdout)
end
