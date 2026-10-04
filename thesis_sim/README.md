# Thesis simulations: radius-scheduled mission speed vs the original Durmaz2024 law, and safety filters

Everything needed to reproduce the thesis results.

```matlab
cd thesis_sim
runAll          % H1, H2, S (~10 min with a parallel pool)
```

**Run counts.** These are supporting simulations, not Monte Carlo statistics: 2 RSC chains per experiment (`cfg.seeds = 1:2` in `thesisConfig.m`). For statistics, set e.g. `cfg.seeds = 1:20`; only the number of runs changes.

**Redrawing only.** Each `exp_*.m` has a `rerun` flag at the top; set it to `false` to redraw the figures from `results/*.mat` without simulating.

**Dependencies.** The vessel models (`tugboat3d.m`, `otter3d.m`) and `getScenario.m` are in the parent folder.

## Files

| File | Purpose |
|---|---|
| `thesisConfig.m` | All shared parameters (map, chains, vessel, controller, filters) |
| `buildChain.m` | RSC funnel chain for a given seed (`RSC.m` as a function), cached in `cache/` |
| `simulateRun.m` | One run: controller, filters, disturbances, vessel model, metrics |
| `runBatch.m` | Runs a list of simulations in parallel and returns a table |
| `exp_H1.m` | H1: radius-scheduled mission speed vs the original Durmaz2024 law |
| `exp_H2.m` | H2: safety filters under unmodeled disturbances |
| `exp_S.m` | S: sensitivity of the radius-scheduled mission speed |
| `runAll.m` | Runs all three |
| `wilsonCI.m`, `saveFigure.m` | 95 % confidence intervals of shares; figure export |
| `results/` | `*.mat`, `*_runs.csv` (one row per run), `*_summary.csv`, `H2_paired.csv` |
| `figures/` | All figures (PNG, 200 dpi) |

## Simulation setup (identical for all experiments)

### Map and chains

**Map.** Tug-scale harbour, 120 x 80 m (`getScenario(6)`):
- an open approach in the south with two rocks;
- a staggered breakwater with a narrow dog-leg entrance;
- a basin with three piers and two moored hulls.

Start at (10, 10) m; goal berth in the eastern slip at (108, 68) m.

**RSC chains.** Built with the parameters of `RSC.m`: 2 m safety margin, minimum funnel area 3 m², child center at 0.9 R_parent. Seeds 1 and 2 give chains of 51 and 46 funnels, with radii from 1.0 to 12 m and path lengths of 173 and 162 m.

`buildChain.m` differs from `RSC.m` in two ways:
- It uses the exact point-to-edge distance to obstacles. `RSC.m` sampled 400 points per edge, which overestimates the distance slightly; next to the long piers that made funnels touch the buffered obstacles.
- It keeps a 1 cm tolerance on the radius.

### Vessels

**Tugboat.** `tugboat3d.m`, the 1/40 Pacific Islander tug of Erünsal (2015): 0.9 m long, 10.2 kg, two stern thrusters.
- Thrust per thruster: +26 N (measured), −14.5 N reverse (assumed).
- Thruster time constant: 0.25 s (assumed).
- Top speed: 2 m/s reported; the model reaches 3.7 m/s at full thrust.

**Otter.** `otter3d.m`, Fossen's Otter USV: 2 m, 55 kg, +120 / −67 N per propeller, top speed 3.1 m/s.

Integration is RK4 with dt = 0.01 s.

### Funnel controller

Durmaz et al. (2024), Eq. 33:

```
u = s(ρ) cos α,    w = Ka α + (u/ρ) sin α,    Ka = 0.3 s⁻¹
```

The active funnel is the highest-index funnel of the chain that contains the measured position; if it is in none, the previous one is kept.

### The two speed laws compared

| Name | Intermediate funnels | Goal funnel |
|---|---|---|
| Durmaz (original, as published) | `s = 2 Kv ρ`, `Kv = 0.05` | same |
| MSR (radius-scheduled mission speed) | `s = U_k = min(U, R_k / T_R)`, `T_R = 1/Ka = 3.3 s` | `s = max(U_k tanh(2 Kv ρ/U_k), min(U_k, 0.5 m/s))` |

- **Same kinematic guarantee.** Both keep the unicycle guarantees `ρ̇ = −s cos²α ≤ 0` and `α̇ = −Ka α`, because `s ≥ 0`.
- **What differs is the speed profile.**
  - The original law's speed is proportional to the distance to the funnel center, so it falls toward zero in every funnel.
  - MSR cruises at a constant speed inside each funnel. That speed changes only between funnels, according to their size.
  - In the goal funnel MSR slows down smoothly, but never below 0.5 m/s, so it can close the last metres against a current.

**Why `U_k ≤ R_k Ka`.**
1. Under the Durmaz law the heading error decays as `α(t) = α₀ e^{−Ka t}`, so the vessel needs about `1/Ka` seconds to turn towards a new center.
2. During that time it travels about `U/Ka` while still pointing away from the center.
3. On the unicycle this cannot increase `ρ`, because `u = s cos α`. On a real vessel, sway and the lag of the yaw-rate loop let it drift outward during the turn, by an amount that grows with the distance travelled.
4. Keeping that distance within the funnel radius, `U/Ka ≤ R_k`, gives `U_k ≤ R_k Ka = R_k / T_R`.

The sensitivity study (S) shows how much margin this rule has.

### Low-level control

- PI loops on surge speed and yaw rate use the `lowLevelControl.m` gains, scaled by surge mass and yaw inertia for the Otter.
- The integrators are frozen while a thruster saturates (conditional integration).
- Thrust allocation is `F_L,R = X/2 ± N/(2d)`, saturated to the thrust limits.

### Safety filters (H2)

All use the barrier `b = R_active − ρ` with `k1 = k2 = 5`. They are QPs with a heavily penalized slack, so they always return the input that violates the condition least.

| Filter | Acts on | Condition |
|---|---|---|
| CBF (kinematic) | references `(u, w)` | CBF in `u` (relative degree 1), HOCBF in `w` (relative degree 2); measured sway as drift |
| HOCBF (dynamic) | thrust `(F_L, F_R)` | second-order HOCBF with the vessel model `M`, `C(ν)`, `D(ν)` |
| rCBF, rHOCBF (robust) | as above | `ḃ` replaced by its worst case `ḃ − V_b`, so the conditions hold for any current up to `V_b = 0.3 m/s` |

### Disturbances

None of them is known to the controller or the filters:
- **current:** constant, added to the ground velocity;
- **INS:** white noise on the measured position, heading and velocities;
- **actuator:** one thruster at 60 % efficiency, plus 2 N thrust noise.

### Metrics

All metrics are computed on the true state.

| Metric | Definition |
|---|---|
| exit | Maximum distance outside the active funnel. The active funnel is chosen with the same rule, from the true position. A run "leaves a funnel" if exit > 1 cm. |
| reached | Within 2 m of the goal. |
| stuck | No new funnel entered for 300 s. |
| collision | Logged position inside an obstacle. None occurred in any run. |
| CV(u) | std(u)/mean(u) in the intermediate funnels. |
| filter active / infeasible | Share of steps where the filter changed the command, or where its condition could not be met within the input limits. |

## Results

### H1 — Radius-scheduled mission speed vs the original Durmaz2024 law

> **H1:** Replacing the original speed law by a mission speed scheduled by funnel radius keeps the vessel inside the funnels, like the original law, but removes its stop-and-go speed profile, shortens the travel time and lets the vessel reach the goal against a current.

Figures: `H1_trajectory.png`, `H1_speed.png`, `H1_nominal.png`, `H1_current.png`.

**No disturbance** (2 chains per entry):

| | Durmaz | MSR 1 m/s | MSR 1.5 m/s | MSR 2 m/s |
|---|---|---|---|---|
| Unicycle: runs leaving a funnel (MSR 3 m/s also 0) | 0 | 0 | — | 0 |
| Tugboat: runs leaving a funnel | 0 | 0 | 0 | 0 |
| Tugboat: time to goal | 935 s | 213 s | 182 s | 173 s |
| Tugboat: mean speed | 0.17 m/s | 0.76 m/s | 0.89 m/s | 0.95 m/s |
| Tugboat: CV(u) | 0.91 | 0.31 | 0.46 | 0.56 |
| Otter (MSR 1 / 2 / 3 m/s): time to goal | 968 s | 215 s | — | 176 s / 174 s |
| Otter: runs leaving a funnel | 0 | 0 | — | 0 |

**Tugboat against a current** (8 directions x 2 chains = 16 runs per entry):

| Current | Durmaz: reached / leaving a funnel | MSR 2 m/s: reached / leaving a funnel |
|---|---|---|
| 0.15 m/s | 0 % / 0 % | **100 %** / 0 % |
| 0.30 m/s | 0 % / 0 % | **75 %** / 25 % (max 1.10 m) |

**Supported.**
- **Safety.** Neither law leaves a funnel on the unicycle, the tugboat or the Otter.
- **Speed.** MSR reaches the goal 4.4-5.4 times faster on the tug, at 4.5-5.7 times the mean speed.
- **Speed profile.** The original law's speed falls toward zero in every funnel (CV ≈ 0.9; `H1_speed.png`). MSR holds a constant speed inside each funnel and changes it only between funnels of different size (CV 0.3-0.6). Its CV grows with U because the scheduling slows it more in the narrow entrance relative to the open water.
- **Current.** The original law never reaches the goal against a current, because its speed near every center is below the current. MSR does.
- **Limit.** Under the strong current MSR leaves a funnel in a quarter of the runs. That is the motivation for H2.

### H2 — Safety filters under unmodeled disturbances

> **H2:** Under disturbances that the controller does not know (current, INS noise, actuator faults), a barrier-function safety filter on top of MSR reduces funnel exits.

Controller: MSR at U = 2 m/s. Figures: `H2_summary.png`, `H2_activity.png`, `H2_example.png`.

| Current 0.30 m/s (16 runs per filter) | none | CBF | HOCBF | rCBF | rHOCBF |
|---|---|---|---|---|---|
| Runs leaving a funnel | 25 % | 25 % | 25 % | 19 % | 12.5 % |
| Worst exit [m] | 1.10 | 0.46 | 0.26 | **0.12** | 0.24 |
| Runs reaching the goal | 75 % | 75 % | 62.5 % | **94 %** | 56 % |
| Runs stalled | 19 % | 19 % | 31 % | 6 % | 31 % |
| Filter active [% of steps] | — | 0.3 | 8.8 | 1.0 | 20.7 |

**Partly supported.**
- **When exits happen.** With MSR, the nominal case, a 0.15 m/s current, INS noise (0.1-0.3 m, 1-3°) and a thruster at 60 % cause no exits for any controller. Only the 0.3 m/s current does.
- **Nominal filters** shrink the worst exit (1.10 → 0.46 / 0.26 m), but not how often runs leave a funnel.
- **The robust CBF** does best overall: the smallest exits (0.12 m), fewer runs leaving, and more runs reaching the goal (94 % vs 75 %).
- **The robust HOCBF** has the fewest exits, but it acts in 21 % of the steps and stalls 31 % of the runs. It is safe, at the cost of progress.
- **Cost without a disturbance.** The robust filters also act occasionally (< 0.4 % of steps) when there is no disturbance, without slowing the vessel.

### S — Sensitivity of MSR

Figure: `S_sensitivity.png`. Tugboat, MSR at U = 2 and 3 m/s, no disturbance, 2 chains per entry.

| Variant | Runs leaving a funnel | Time to goal (U = 2 / 3 m/s) |
|---|---|---|
| Reverse thrust −7.25 / −14.5 / −26 N | 0 in every case | 172-174 s / 166-168 s |
| Thruster time constant 0.1 / 0.25 / 0.5 s | 0 in every case | 171-175 s / 165-169 s |
| T_R = 0.5/Ka | 0 | **126 s / 109 s** |
| T_R = 1/Ka (default) | 0 | 173 s / 167 s |
| T_R = 2/Ka | 0 | 301 s / 301 s |

- **The assumed vessel values do not change the result.** No exits occur in any variant, and the travel time stays within a few seconds.
- **The scheduling constant is a safety/speed trade-off.** Half of `1/Ka` is still exit-free here and 27-35 % faster. With `2/Ka` every funnel is speed-limited by its radius, so U no longer matters. `T_R = 1/Ka` is therefore a conservative choice.

## Overall conclusion

1. **MSR vs the original law (H1).** The radius-scheduled mission speed keeps the vessel inside the funnels as reliably as the original Durmaz2024 law. It removes the original law's stop-and-go speed profile, reaches the goal 4-5 times faster, and keeps working against a current where the original law cannot reach the goal at all.
2. **Robustness of MSR (H1, S).** MSR is robust to the assumed thruster values. The radius rule `U_k = min(U, R_k Ka)` has margin, so `T_R` can be tuned for speed.
3. **Safety filters (H2).** They are only needed for disturbances strong enough to push the vessel out, here a 0.3 m/s current. The robust CBF is the best compromise; the robust HOCBF trades progress for safety.

## Known limitations

- **Small run counts.** These are supporting simulations (2 chains); percentages are indicative.
- **Vessels.** The tugboat is the main vessel; the Otter appears only in the H1 nominal comparison, and CyberShip II is not used.
- **Assumed values.** The tug's reverse thrust and thruster time constant are assumed; S shows they do not change the results.
- **Informal radius rule.** It is an engineering rule (derivation above, margin shown by S); there is no formal invariance proof for the vessel dynamics.
- **Collision check.** It tests the vessel's reference point, not its hull outline.
