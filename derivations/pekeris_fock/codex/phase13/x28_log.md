# TASK 28 execution log — 2026-10-24

Read TASK_codex_28, x18/x24/x26/x27 reports and Claude's SOps implementation.
No agents, source modifications, large solves, dependency installation or commits.
All output files are phase13/x28_*; shared logs receive appended entries only.

Derived the exact zonal single-loop difference and an explicit four-term Appell
F1 expression for its rotation average. Wrote a serial independent checker.
The initial checker stopped on Julia's top-level loop scoping of gamma; an
explicit global assignment fixed the diagnosed issue. Final rerun passes:
angular/Appell error 1.33e-15, lifted branch-difference error 2.58e-17,
three-mode weak spectral errors <=4.12e-13, endpoint errors <=2.3e-15.
Peak RSS 384,122,880 bytes (<2 GB), final run about 3 seconds. No result relies
on the nonconverged point-kernel diagnostic. No ore_algebra computation is claimed.

Task A: explicit elliptic reduction of the branch-variation angular kernel DONE;
identification/evaluation of K32 itself OPEN, not passed off as established.
Task B: exact Delta G DONE, with actual-difference versus normalized-N conventions.
Task C: feasibility audit DONE; practical fallback co-signed, universal
impossibility and operator-size lower bounds not claimed.

FINAL SUMMARY: Delta G=(4z^3-3z)/(2pi sqrt(1-z^2)); its rotation average has an explicit elliptic/Appell F1 form; K32 identification and B_A annihilator remain OPEN; practical small-hours telescoping fallback co-signed with qualifications.
