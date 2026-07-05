# Task 2.2 anomaly diagnostics (analysis stage)

These JSONs are **not** part of the convergence ladder. They are short,
throwaway probes run during the Task 2.2 analysis stage to characterize the
deterministic `secondOrderLimited` CVODE failure at the finest tolerance
rungs (see the spec's "Phase 2 findings" addendum, §P2.4). No `src/` changes
were made to produce them.

| file | what it probes | result |
|---|---|---|
| `diag_strained_sol_rung5_noregrid.json` | rung-5 tolerances with `Debug(regridding=False)` | FAILED, byte-identical grid trajectory to the production failure — the `regridding` flag is only a verbose-print toggle (`src/debugUtils.h`), so this is a **determinism reproduction**, not a frozen-grid test |
| `diag_strained_sol_rung5_frozengrid.json` | rung-5 tolerances with `regridStepInterval=1e9`, `regridTimeInterval=1e30` (grid frozen after initial adaptation) | **COMPLETED**, N=100, bounded convection cost (~42 CVODE steps/global step). Under-resolved answer (peak_T/consumption_speed are off) — expected, since the point is to test integrator robustness on a non-refining grid, not accuracy |
| `diag_strained_sol_rung4p5.json` | intermediate "rung 4.5": `vtol=0.040, dvtol=0.065, gridMax=5.0e-5`, adaptation on | **COMPLETED**, N=269, converged scalars (consumption_speed 0.35915, peak_T 1558.92) — a stable grid fixed point exists just below rung 5 |

Work dir for these runs: `build/test/convergence-work-diag/` (not committed).
