# TASK_codex_12 independent step 1 — saved before claude6 inspection

All branch/constant tests completed at (10,17) and (13,22), with prime
1004535809 (different from phase five/Claude's announced prime). The
10x45 boundary-relation arrays are invariant under all eight independent
radical choices. See goal2x_constant_field.md for the pi/i/sqrt2 results.

Independent line t=10+u, q=17+2u: 261 distinct split regular samples saved
in goal2x_line.json. Every sample is checked in all 240 represented
coordinates. Bivariate Taylor order 9 is used, including coefficient-row
derivatives and actual endpoint/zero-jet quotienting.

Degree tests: at 32 and 128 training samples, minimum-total-degree regular
interpolants for ALL 450 entries have balanced degree pairs (16,15) and
(64,63), respectively; NONE passes the five held-out samples. At 256
training samples, the 40 entries with basis columns 0,1,5,44 in all ten
boundary rows have (128,127); again NONE passes five holdouts. Thus no
empirically stable degree has yet been determined. The exact data give
finite-sample total-degree lower bounds 127 (all entries) and 255 (those
40 entries); the displayed balanced degrees themselves are interpolation
artifacts, NOT the true degrees. No low-degree reconstruction is claimed.

This step's independent experiments are saved, including their negative
outcome, before opening claude6. The remaining uncertainty is the actual
minimal numerator/denominator degrees; the next comparison may supply
higher-degree candidates that can be tested without changing this record.

All computation serial; sample time about 1.7 seconds after exact arithmetic
acceleration. No Claude-owned file modified. No outside-directory write
was requested by any script. No subprocess or agent parallelism.
