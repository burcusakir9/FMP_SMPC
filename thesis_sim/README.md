# Thesis simulations: mission speed, radius scheduling and safety filters in RSC funnel chains

Everything needed to reproduce the thesis results for hypotheses H1-H3.

```matlab
cd thesis_sim
runAll          % H1, H2, H3 (~15 min with a parallel pool)
```

Each `exp_H*.m` has a `rerun` flag at the top; set it to `false` to redraw the figures from `results/H*.mat` without simulating.

## Files

| File | Purpose |
|---|---|
| `thesisConfig.m` | All shared parameters (map, chains, vessel, controller, filters) |
| `buildChain.m` | RSC funnel chain for a given seed (`RSC.m` as a function), cached in `cache/` |
| `simulateRun.m` | One run: controller, filters, disturbances, vessel model, metrics |
| `runBatch.m` | Runs a list of simulations in parallel and returns a table |
| `exp_H1.m`, `exp_H2.m`, `exp_H3.m` | One experiment per hypothesis |
| `saveFigure.m` | Saves figures to `figures/` |
| `results/` | `H*.mat`, `H*_runs.csv` (one row per run), `H*_summary.csv`, `H3_paired.csv` |
| `figures/` | All figures (PNG, 200 dpi) |

Vessel models, `getScenario.m` and `map.kml` are in the parent folder.

## Simulation setup (identical for H1-H3)

**Map and chains.** Real coastline map (`getScenario(4)`, 400 x 400 m). Five RSC funnel chains, from seeds 1-5, built with the parameters of `RSC.m`: 2 m safety margin, minimum funnel area 3 m², child center at 0.9 R_parent. Each chain has 44-66 funnels, with radii from 1.1 m to 92.5 m.

**Vessel.** `tugboat3d.m`, the 1/40 Pacific Islander tug of Erünsal (2015): 0.9 m long, 10.2 kg, two stern thrusters. Thrust per thruster is +26 N (measured) and −14.5 N (assumed). Thruster lag is 0.25 s (assumed). The reported top speed is 2 m/s; with full thrust the model reaches 3.7 m/s. RK4 integration with dt = 0.01 s.

**Funnel controller.** Durmaz et al. (2024), Eq. 33:

```
u = s(ρ) cos α,    w = Ka α + (u/ρ) sin α,    Ka = 0.3
```

The active funnel is the highest-index funnel of the chain that contains the vessel's (measured) position; if it is in none, the previous one is kept.

**Speed laws `s(ρ)`.** Any `s ≥ 0` keeps the unicycle guarantees `ρ̇ = −s cos²α ≤ 0` and `α̇ = −Ka α`.

| Name | Intermediate funnels | Goal funnel |
|---|---|---|
| Durmaz (original) | `2 Kv ρ`, `Kv = 0.05` | same |
| MS (mission speed) | `U` | `U tanh(2 Kv ρ / U)` |
| MS + radius scheduling | `U_k = min(U, R_k / T_R)`, `T_R = 1/Ka = 3.3 s` | `U_k tanh(2 Kv ρ / U_k)` |

The radius rule comes from the heading dynamics. The heading error decays with time constant `1/Ka`, during which the vessel travels about `U/Ka`. Keeping that distance within the funnel radius gives `U ≤ R_k Ka`.

**Low-level control.** PI loops on surge speed and yaw rate (gains of `lowLevelControl.m`). Their integrators are frozen while a thruster saturates (conditional integration). Thrust allocation is `F_L,R = X/2 ± N/(2d)`, saturated to the thrust limits.

**Safety filters.** Both use the barrier `b = R_active − ρ` with `k1 = k2 = 5`. Both are QPs with a heavily penalized slack, so they always return the input that violates the barrier condition least.
- **CBF (kinematic):** acts on the references `(u, w)`. It is a CBF in `u` (relative degree 1) and a HOCBF in `w` (relative degree 2). The measured sway velocity is treated as drift.
- **HOCBF (dynamic):** acts on the thrust `(F_L, F_R)`, a second-order HOCBF using the vessel model `M`, `C(ν)`, `D(ν)`.

**Disturbances.** None of them is known to the controller or the filters.
- Current: constant, added to the ground velocity.
- INS: white noise on the measured position, heading and velocities.
- Actuator: thrust efficiency of one thruster, plus thrust noise.

**Metrics.** All are computed on the true state.

| Metric | Definition |
|---|---|
| exit | Maximum distance outside the active funnel. The active funnel is chosen with the same rule, from the true position. A run "leaves a funnel" if exit > 1 cm. |
| left chain | The run was outside every funnel of the chain at some point. |
| reached | Within 2 m of the goal. |
| stuck | No new funnel entered for 300 s. |
| collision | Logged position inside an obstacle. |
| CV(u) | std(u)/mean(u) in the intermediate funnels. |
| filter active / infeasible | Share of steps where the filter changed the command, or where its condition could not be met within the input limits. |

## Results

### H1 — Mission speed keeps the kinematic guarantee and removes the stop-and-go speed profile

Results: `results/H1_summary.csv`. Figures: `H1_speed.png`, `H1_nominal.png`, `H1_current.png`.

| | Durmaz | MS 1 m/s | MS 2 m/s |
|---|---|---|---|
| Unicycle: runs leaving a funnel (also MS 3 m/s: 0 %) | 0 % | 0 % | 0 % |
| Tugboat: time to goal | 1163 s | 946 s | 505 s |
| Tugboat: CV(u) | 0.98 | 0.06 | 0.09 |
| Tugboat: runs leaving a funnel | 0 % | 0 % | 20 % (max 2.7 cm) |
| Current 0.15 m/s: runs reaching the goal | 5 % | — | 100 % |
| Current 0.30 m/s: runs reaching the goal | 0 % | — | 10 % |

**Supported.**
- On the ideal unicycle no law ever leaves a funnel, at any speed, as the theory predicts.
- On the tugboat the original law's speed swings between about 0 and 3.7 m/s at every handoff (CV ≈ 1). Mission speed keeps it nearly constant (CV < 0.1) and halves the travel time at 2 m/s.
- The original law also cannot reach the goal against a current: it slows down near every funnel center, where even 0.15 m/s of current wins.
- **Limitation:** at 0.3 m/s of current every law, including mission speed, stalls in the goal funnel, because the goal-funnel slowdown takes the speed below the current.

### H2 — Radius scheduling removes the funnel exits

Results: `results/H2_summary.csv`. Figures: `H2_sweep.png`, `H2_mechanism.png`, `H2_example.png`.

| U [m/s] | 1 | 1.5 | 2 | 2.5 | 3 | 3.5 |
|---|---|---|---|---|---|---|
| MS: runs leaving a funnel | 0 % | 0 % | 20 % | 40 % | 40 % | 60 % |
| MS: worst exit [m] | 0 | 0 | 0.03 | 1.42 | 0.50 | 0.97 |
| MS + radius scheduling: runs leaving a funnel | 0 % | 0 % | 0 % | 0 % | 0 % | 0 % |
| Time: MS / scheduled [s] | 946 / 951 | 651 / 661 | 505 / 520 | 420 / 438 | 372 / 387 | 335 / 351 |

**Supported.**
- With a constant mission speed, exits appear from 2 m/s upward. They occur only in funnels with `R_k / U < T_R` (`H2_mechanism.png`).
- Scheduling the speed by funnel radius removes every exit at every speed, on all five chains, for 1-5 % more travel time.
- Speeds above 2 m/s are beyond the tug's reported top speed but within the model's capability. They are included to show where exits start.
- No run collided with the coastline.

### H3 — Safety filters under unmodeled disturbances

Controller: radius-scheduled mission speed at U = 2 m/s. Results: `results/H3_summary.csv`, `results/H3_paired.csv`. Figures: `H3_summary.png`, `H3_activity.png`, `H3_example.png`.

| Disturbance (runs per filter) | Runs leaving a funnel (none / CBF / HOCBF) | Worst exit [m] (none / CBF / HOCBF) |
|---|---|---|
| nominal (5) | 0 / 0 / 0 % | 0 / 0 / 0 |
| current 0.15 m/s (40) | 2.5 / 2.5 / 2.5 % | 0.08 / 0.04 / 0.06 |
| current 0.30 m/s (40) | 10 / 12.5 / 10 % | **1.15 / 0.15 / 0.16** |
| INS low and high (25 + 25) | 0 / 0 / 0 % | 0 / 0 / 0 |
| actuator, left or right at 60 % (5 + 5) | 0 / 0 / 0 % | 0 / 0 / 0 |

**Partly supported.**
- With radius scheduling, INS noise and thruster faults never push the tugboat out of a funnel. The filters have nothing to correct there: they are active in < 0.01 % of steps.
- Only the strong current causes exits. Both filters then cut the worst exit about 7× (1.15 m → 0.15 m), and the mean exit 4-5×.
- They do not change how many runs leave a funnel: in the paired comparison, 2-3 of 40 runs improve and 0-1 get worse.
- In the worst case (`H3_example.png`) the HOCBF prevents the exit by stopping at the funnel boundary until the run stalls. It is safe, but it makes no progress.
- The filters' conditions were almost never unsatisfiable (≤ 0.005 % of steps), so the thrust limits are not what limits them here.

## Overall conclusion

1. Mission speed keeps the kinematic funnel guarantee, removes the stop-and-go speed profile of the original law, and is much faster (H1).
2. On the real vessel, safety comes from scheduling the cruise speed by funnel radius (`U_k = min(U, R_k Ka)`). It removes every exit in nominal conditions at a few percent travel-time cost (H2), and with it INS noise and thruster faults cause no exits either (H3).
3. Barrier-function filters are a secondary layer. They reduce the size of exits caused by a strong current, but not how often they happen, and they can trade progress for safety (H3).

## Known limitations

- **Tugboat only.** The Otter and CyberShip II models exist (`otter3d.m`, `cybership3d.m`) but are not used yet.
- **Assumed values:** the reverse thrust limit and the thruster time constant.
- **Speeds above 2 m/s** exceed the reported top speed.
- **Goal stall under a strong current:** the goal-funnel slowdown goes below the current speed, so the vessel cannot close the last metres.
- **Filters are disturbance-blind.** Neither filter models the disturbances. A robust CBF (tightening the barrier by a known disturbance bound) is the natural extension.
- **Collision check:** it tests the vessel's reference point, not its hull outline.
