# racing-line-lap-time-optimization

![MATLAB](https://img.shields.io/badge/MATLAB-R2021a+-orange.svg)
![Optimization Toolbox](https://img.shields.io/badge/Optimization_Toolbox-fmincon_SQP-blue.svg)
![Global Optimization](https://img.shields.io/badge/Global_Optimization-Particle_Swarm-9cf.svg)
![Curve Fitting](https://img.shields.io/badge/Curve_Fitting-makima_splines-lightgrey.svg)

> **Executive Summary:** A quasi-steady-state racing line and minimum-lap-time optimization tool built on discrete point-force vehicle dynamics. The circuit is reduced to a one-dimensional lateral-deviation problem in a curvilinear (Frenet-style) coordinate system, discretized with curvature-driven **Adaptive Mesh Refinement**, and reconstructed with boundary-safe **Modified Akima (makima)** splines. A multi-phase strategy combines a geometric **curvature** pre-optimization, a physics-based **Sequential Quadratic Programming (SQP, `fmincon`)** lap-time solver, optional **Particle Swarm Optimization (PSO)** seeding, and a final perturbation search. Demonstrated on the Hungaroring with an F1-class vehicle model, the pipeline drives a 1:29.4 geometric baseline down to a best lap of **1:15.950**.

Barnabas Szeibert · Panagiotis Sachinis — April 2026

---

## Introduction

In motorsport, extracting maximum performance from a given vehicle is as critical as designing the vehicle itself. The single variable a driver can most readily adjust is the *path* taken through a circuit — the **racing line** — chosen to maximize average velocity and therefore minimize lap time.

Finding this line is a high-dimensional problem coupling track geometry (width, corner radius, sector dependency) to vehicle dynamics (power and braking limits, aerodynamic and tire forces). This project decouples the problem into two stages: a purely **geometric** optimization that minimizes integrated curvature to find a smooth, drive-able baseline, followed by a **physical** optimization that minimizes true lap time under point-force dynamic constraints. The result is a modular MATLAB pipeline that can be pointed at any circuit for which centerline and track-width data are available.

## Modeling Methodology

### Track Representation

The track centerline is parameterized by arc length $s \in [0, S]$ and described by a Cartesian position vector $r_c(s)$. An orthogonal unit normal $n(s)$ (pointing to the **left** boundary) defines the lateral domain, bounded by the left/right width functions $W_l(s)$ and $W_r(s)$. The racing line is then a single spatially-dependent design variable $\alpha(s)$ — the lateral deviation of the center of gravity from the centerline:

$$r_{race}(s) = r_c(s) + \alpha(s)\,n(s)$$

This transformation collapses a blind two-dimensional routing problem into a **one-dimensional** optimization: the solver only chooses $\alpha$ at each station $s$, and the curvilinear domain becomes a naturally convex hyper-rectangle.

### Adaptive Node Placement

Straights need almost no spatial resolution while corners demand a lot, so a uniform node spacing is wasteful. An **Adaptive Mesh Refinement (AMR)** scheme builds a density function from the track's own curvature, pulling control points off the straights and clustering them around braking zones and apexes. Raw centerline curvature $\kappa(s)$ is Gaussian-smoothed, normalized, sharpened by a focusing exponent, and integrated to a cumulative density that is inverse-sampled into the final control grid:

$$\rho(s) = 1 + w_c\,\hat{\kappa}(s)^{p}, \qquad C(s) = \int_0^s \rho(\tau)\,d\tau$$

### Trajectory Interpolation

A continuous path is reconstructed from the discrete nodes with piecewise cubic interpolation. Standard cubic splines enforce $C^2$ continuity globally, which causes spatial **overshoot** that can violate track limits even when every control point is legal. To prevent this, the model uses **Modified Akima (`makima`)** Hermite splines, which relax $C^2$ continuity and compute nodal derivatives locally. Consequently, if the optimizer bounds the discrete points, the entire continuous path is guaranteed to stay inside the drivable surface:

$$-(W_r(s) - c_m) \le P_i(s) \le (W_l(s) - c_m) \quad \forall s \in [0, S]$$

### Physical Model

The vehicle is reduced to its fundamental dynamic constraints: a **point mass** at the center of gravity, a **constant friction coefficient** $\mu$, and **velocity-dependent aerodynamics**. Downforce and drag scale with $v^2$; the lateral grip ceiling, friction-circle coupling, and powertrain limit are then combined in a multi-pass velocity solve.

* **Aerodynamic & normal forces:** $F_{DF} = \tfrac{1}{2}\rho\,C_{LA}\,v^2$, $\;F_{Drag} = \tfrac{1}{2}\rho\,C_{DA}\,v^2$, $\;F_z = mg + F_{DF}$.
* **Lateral grip ceiling:** $v_{limit} = \sqrt{\dfrac{\mu m g}{m|\kappa| - \tfrac{1}{2}\mu\rho\,C_{LA}}}$ — corners where downforce out-scales the cornering demand become drag/engine-limited rather than grip-limited.
* **Friction circle:** remaining longitudinal force $F_{lon} = \sqrt{\max(0,\,(\mu F_z)^2 - F_{lat}^2)}$.
* **Powertrain limit:** $F_{engine} = \min(P_{max}/v,\ F_{gear\_max})$.
* **Multi-pass integration:** the final velocity at each node is the geometric minimum of the grip ceiling, a backward (braking) pass, and a forward (acceleration) pass over the closed loop.

### Vehicle Parameters

Parameters approximate a modern Formula 1 car under the 2023 technical regulations (`carParams.m`).

| Parameter | Notation | Value |
|---|---|---|
| Vehicle mass | $m$ | 798 kg |
| Peak engine power | $P_{max}$ | 735.5 kW |
| Maximum tractive force | $F_{gear\_max}$ | 10 kN |
| Minimum turning radius | $R_{min}$ | 5.0 m |
| Coefficient of lift × area | $C_{LA}$ | 4.8 m² |
| Coefficient of drag × area | $C_{DA}$ | 1.2 m² |
| Air density | $\rho$ | 1.225 kg/m³ |
| Ground friction coefficient | $\mu$ | 1.7 (set in `genTrack.m`) |

## Optimization Problems

Both problems are box-constrained nonlinear programs with linear equality constraints (negative-null form), over the design vector $u = [\alpha_{ctrl,1}, \dots, \alpha_{ctrl,N_{ctrl}}]^T$. Box constraints come directly from the track width minus a `car_margin` of 0.5 m, and the equality constraints enforce $C^0$ position and $C^1$ heading continuity across the start/finish line.

**1. Curvature optimization** minimizes integrated squared curvature with a small arc-length penalty so the optimizer does not simply chase the widest possible line:

$$f_{curv}(u) = \sum_{i=1}^{N_{ctrl}} \left(\kappa_i^2\,\Delta s_i + w_{length}\,\Delta s_i\right)$$

**2. Physical lap-time optimization** minimizes true traversal time with a steering-jerk smoothness penalty (arc length is now implicit in the time objective):

$$f_{phys}(u, v) = \sum_{i=1}^{N_{ctrl}} \frac{\Delta s_i}{v_i} + w_{smooth}\sum_{j=1}^{N_{ctrl}} \left(\Delta^2 \alpha_{ctrl,j}\right)^2$$

### Choice of Solvers

The design space is high-dimensional, non-linear, and (for the physics case) non-convex, with a broad convex-like "minimal region" containing many shallow local minima. No single method captures it efficiently, so a **multi-phase** approach is used:

1. **Curvature SQP** — `fmincon` with the SQP algorithm solves the near-convex geometric problem to produce a strong baseline line.
2. **Physics SQP** — the curvature-optimized path seeds a second SQP solve on the full physics objective, landing inside the minimal region.
3. **PSO seeding (alternative)** — a custom or built-in Particle Swarm Optimizer provides a non-geometric initial guess; used to probe the solver's sensitivity to initialization.
4. **Perturbation search** — the best physics result is nudged by constant, sinusoidal, and stochastic perturbations and re-optimized to hunt for a stronger nearby minimum.

The **Particle Swarm Optimizer** updates each candidate line from its inertia, cognitive best, and the swarm's global best, with positions clamped to the track box and velocities capped to prevent divergence:

$$v_i(k+1) = w\,v_i(k) + c_p r_p\,(p_{best,i} - u_i(k)) + c_g r_g\,(g_{best} - u_i(k))$$

## Results

All runs use the **Hungaroring** (Budapest) — a short, sector-rich circuit, ideal for keeping the reduced dimension manageable.

| Stage | Method | Lap time |
|---|---|---|
| Geometric baseline | Curvature SQP | **1:29.395** |
| Physics, curvature seed | SQP (`fmincon`) | **1:16.264** |
| Physics, PSO seed | PSO → SQP | ≈ 1 s slower |
| Final, perturbed | Sine ($\omega/4$) → SQP | **1:15.950** |

The geometric optimizer reproduces the classic **out-in-out** line, using full track width on entry, kissing the apex, and tracking out wide — but its derived velocity profile is erratic because no longitudinal coupling is applied. Introducing the physics model smooths the telemetry dramatically: minimum corner speeds rise, the car carries wider lines through long corners to exploit downforce, and there are no oscillatory (neither-accelerating-nor-braking) regions. The PSO-seeded variant is macroscopically similar but occasionally fails to fully use the track at complex corners (e.g. Turn 9), because `fmincon` is a local solver and merely polishes whatever guess it inherits. Across repeated runs the stochastic PSO seed swung lap times by roughly ±5 s, which is why the curvature seed is the primary method. Finally, small perturbations of the curvature-seeded result consistently found marginally better minima, with a one-quarter-frequency sinusoid yielding the best **1:15.950**.

> **Reproducing the figures.** Each `main_*.m` script renders the racing line (velocity- and acceleration-coloured), the along-track telemetry, and convergence/sensitivity plots. Export them into a `figs/` directory if you want to embed result images in this README.

## Limitations

The point-force reduction is a deliberate simplification, so the absolute lap time should **not** be compared quantitatively to real-world data. Notable gaps include the absence of transient dynamics, weight transfer, slip angles, thermal tire behavior, and three-dimensional inertia, plus static (velocity-only) aerodynamics. Because the solvers are gradient-based and the physics design space is non-convex, **a true global minimum cannot be guaranteed** — the multi-phase search finds a strong minimum, not a proven optimum. Results are also highly sensitive to the number of control points $N_{ctrl}$ and the smoothness/length penalty weights, so these were tuned empirically rather than via a formal per-variable sensitivity study.

Future work should attack dimensionality directly — e.g. a sector-by-sector decomposition that shrinks each sub-problem enough to make metaheuristics viable — before layering in higher-fidelity vehicle physics.

---

## Repository Structure

The top-level scripts are entry points; reusable physics and geometry helpers live under `functions/`, and circuit data under `tracks/`.

```
.
├── main_lineOpti.m          # Geometric curvature optimization (SQP baseline)
├── main_carOpti.m           # Full physics lap-time SQP, curvature-seeded (primary)
├── main_carOptiPSO.m        # Physics optimization via MATLAB built-in particleswarm
├── main_hybrid_customPSO.m  # Hybrid: custom PSO (coarse) → fmincon (fine)
├── psoopti.m                # Standalone run of the custom PSO optimizer
├── carParams.m              # F1-class vehicle parameter set
├── functions/               # Physics & geometry helpers
│   ├── genTrack.m                       # Loads a track CSV, builds bounds, curvature, mu
│   ├── genAdaptiveNodes.m               # Curvature-driven adaptive node placement (AMR)
│   ├── calcCurvatureCost.m              # Geometric objective (Eq. for f_curv)
│   ├── calcLapTimeCostDetailFriction.m  # Multi-pass physics lap-time objective
│   └── lineoptifuncs/                   # Lap-time / PSO support functions
│       ├── calcLapTimePSO.m             #   Penalised lap-time objective for PSO
│       └── pso.m                        #   Custom Particle Swarm Optimizer
└── tracks/                  # Circuit coordinate data (x, y, w_right, w_left)
    ├── Budapest.csv         # Hungaroring (default)
    ├── Montreal.csv
    ├── Suzuka.csv
    └── Zandvoort.csv
```

Track CSVs follow the TUM `racetrack-database` format: columns are `x_m, y_m, w_tr_right_m, w_tr_left_m`.

## Installation & Usage

### Requirements

* **MATLAB** (R2021a or newer recommended)
* **Optimization Toolbox** — `fmincon` (SQP)
* **Global Optimization Toolbox** — `particleswarm` (only for `main_carOptiPSO.m`)
* **Curve Fitting Toolbox** — B-spline helpers (`spmak`, `augknt`, `fnval`); `makima` is built into base MATLAB
* **Parallel Computing Toolbox** — optional, accelerates the built-in PSO run

### Running an optimization

1. Open MATLAB in the repository root (so `functions/` and `tracks/` are on the path).
2. Choose a circuit by editing the `FileName` line in `functions/genTrack.m` (default `Budapest.csv`).
3. Run the script for the method you want:

```matlab
% Geometric baseline only
main_lineOpti

% Recommended: full physics lap time, curvature-seeded
main_carOpti

% Alternatives
main_carOptiPSO        % built-in particle swarm
main_hybrid_customPSO  % coarse custom PSO -> fine fmincon
```

Key tunables sit at the top of each `main_*.m`: the control-point count `n_var`, the spline type (`'makima'` or `'bspline'`), the AMR weights (`corner_weight`, `smoothing_window`), the curvature length penalty `weight_length`, and the smoothness penalty `w_smooth`. Each script prints final metrics (line length, integrated curvature, lap time) and renders the racing line and telemetry plots.

> **Note:** the physics objective is non-convex and the gradient solver finds a strong local minimum, not a proven global optimum — expect run-to-run variation when seeding from PSO. For questions, contact the authors.
