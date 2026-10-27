# Semi-Lagrangian Fish Advection Tests

MATLAB code that generates the data and figures for the semi-Lagrangian fish
advection test cases. A Gaussian patch of fish swims up a food gradient,
either horizontally (test 1) or at 45 degrees (test 2). Each case is run
with three methods:

- **A**: semi-Lagrangian advection.
- **RA**: random semi-Lagrangian advection. The swim direction is perturbed
  with von Mises noise.
- **AD**: advection plus Laplacian diffusion (κ = 1000 m²/s).

This code replicates the data and figures used in the preprint paper: 

[A method for directed horizontal movement in Eulerian fish community models, Brzenski J., Petrik C.M. 2026, submitted](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7138904)

## Requirements

- MATLAB. Tested with R2023b.
- Statistics and Machine Learning Toolbox, for `corr`.
- Parallel Computing Toolbox (optional). The speed × Δt sweeps use `parfor`.
  Without the toolbox, `parfor` runs as an ordinary `for` loop.

The cmocean colormaps are bundled in `external/cmocean`.

## Quick start

```matlab
cd repo
run_all              % generate every .mat file, then every figure
run_all('data')      % data only
run_all('figures')   % figures only (needs the .mat files)
```

All output (`.mat`, `.png`, `.svg`, `.eps`) is written to `scripts/`. A full
run takes about 25 minutes with a 30-worker parallel pool, and most of that
is the random advection runs on the diagonal grid. The `.mat` files take about 10MB of drive space. The figures alone take about 1.5 minutes.

To run one script by itself, `cd` into `scripts/` and run it. Each script
adds `../src` to the path itself. Run the data scripts before the figure
scripts that read their files (see the tables below).

## Layout

```text
repo/
├── run_all.m        Runs the data scripts, then the figure scripts
├── src/             Model functions (semi-Lagrangian scheme, diffusion, setup)
├── scripts/         Data generation and figure scripts; output is written here
└── external/
    └── cmocean/     Colormaps (Thyng et al., MIT license)
```

## Data generation

| Script | Writes | Experiment |
| --- | --- | --- |
| `test1_Linear_A.m` | `test1_Linear_A_cfl{050,090,100,125}.mat` | Horizontal, A, CFL = 0.5, 0.9, 1.0, 1.25 |
| `test1_Linear_RA.m` | `test1_Linear_RA_cfl*.mat` | Horizontal, RA, same CFLs |
| `test1_Linear_A_D.m` | `test1_Linear_A_D_cfl*.mat` | Horizontal, AD, same CFLs |
| `test2_Diag_A.m` | `test2_Diag_A_cfl*.mat` | Diagonal, A, same CFLs |
| `test2_Diag_RA.m` | `test2_Diag_RA_cfl*.mat` | Diagonal, RA, same CFLs |
| `test2_Diag_AD.m` | `test2_Diag_AD_cfl*.mat` | Diagonal, AD, same CFLs |
| `test3_CFL_A_error.m` | `CFL_A_error.mat` | Horizontal, A, sweep of 20 speeds × 47 time steps |
| `test3_CFL_RA_error.m` | `CFL_RA_error.mat` | Horizontal, RA, same sweep |
| `test3_CFL_AD_error.m` | `CFL_AD_error.mat` | Horizontal, AD, same sweep |
| `test5_diag_CFL_A_error.m` | `CFL_Diag_A_error.mat` | Diagonal, A, same sweep |
| `test5_diag_CFL_RA_error.m` | `CFL_Diag_RA_error.mat` | Diagonal, RA, same sweep |
| `test5_diag_CFL_AD_error.m` | `CFL_Diag_AD_error.mat` | Diagonal, AD, same sweep |
| `test4_kappa_AD.m` | `test4_kappa_AD.mat` | Horizontal, AD at CFL = 1, κ = 0, 250, 500, 750, 1000 m²/s |

The sweeps cover speeds of 0.05–1 m/s and time steps of 1–12 h.

## Figures

| Script | Writes | Reads |
| --- | --- | --- |
| `horizontal_advect_plot.m` | `test1_plot1_fulla.{png,svg}` | `test1_Linear_{A,RA,A_D}_cfl*.mat` |
| `diagonal_advect_plot.m` | `test2_plot1_fulla.{png,svg}` | `test2_Diag_{A,RA,AD}_cfl*.mat` |
| `analyze_test1_plot2.m` | `test1_plot2b.{png,eps,svg}` | `test1_Linear_{A,RA,A_D}_cfl*.mat` |
| `analyze_test2_plot2.m` | `test2_plot2b.{png,eps,svg}` | `test2_Diag_{A,RA,AD}_cfl*.mat` |
| `error_bar_plot.m` | `error_bar_plot.{png,svg}` | test1 and test2 `.mat` files |
| `cfl_contour_plot.m` | `cfl_contour_plota.{png,svg}` | `CFL_{A,RA,AD}_error.mat`, `CFL_Diag_{A,RA,AD}_error.mat` |
| `analyze_mass_cons.m` | `mass_cons_{timeseries,sweeps,summary}.png` | All of the `.mat` files above |
| `cfl_shift_plot.m` | `cfl_shift_plot.{png,svg}` | None. It runs its own test 1 simulations at CFL = 0.72, 1.00 and 1.08 |
| `plot_von_Mises_dist.m` | `von_mises_distributions{,_bw}.{png,svg}` | None |

`horizontal_initial_conditions.m` (in `scripts/`) builds the test 1 initial
conditions for `horizontal_advect_plot.m`.

## Model functions (`src/`)

| Function | Role |
| --- | --- |
| `SemiLagrangianAdvectFish.m` | One time step of semi-Lagrangian fish advection (A, AD) |
| `RandomSemiLagrangianAdvectFish.m` | One time step of random semi-Lagrangian advection (RA) |
| `semiLagrangianFish.m`, `moveFishFunction.m` | Flux of fish from a cell to its four neighbours |
| `RandomSemiLagrangianFish.m`, `RandomMoveFishFunction.m` | The same, with random weights |
| `randomSpeed.m`, `decompose_angle.m`, `circ_vmrnd.m` | Von Mises perturbation of the swim direction |
| `percentMoreFoodFourWay.m`, `MoreFoodFourWay.m` | Food gradient in the four directions |
| `ApparentCurrentFull.m` | Swim speed plus ocean current |
| `laplacian2D.m` | Laplacian diffusion step (AD) |
| `makeFood.m` | Initial food and fish fields (gradient, Gaussian hump) |
| `AdvectFood.m`, `eatFoodDomain.m`, `smooth2D.m`, `zeroCornerMatrix.m` | Helpers |

## Notes

- **Random numbers.** The A and AD results are deterministic.
  - `test1_Linear_RA.m` and `test2_Diag_RA.m` do not set a seed, so their
    results depend on MATLAB's random state when they start. `run_all`
    resets the generator with `rng('default')` before every script. Each
    script then gives the same result as when it is run alone in a fresh
    MATLAB session. This is different from the images generated for the 
    paper, but, the differences should be minimal.
  - `test3_CFL_RA_error.m` and `test5_diag_CFL_RA_error.m` call `rng(42)`.
    That seeds only the MATLAB client, not the `parfor` workers, and each
    iteration draws from whichever worker runs it. Their results change
    slightly from run to run when a parallel pool is used.
- **Figure size in batch mode.** The two-column figures are 1642 × 991
  pixels. Under `matlab -batch`, MATLAB limits a figure to the size of the
  virtual screen, so it can come out a few pixels narrower.

## Third-party code

- `external/cmocean/`: cmocean colormaps by Kristen Thyng, MATLAB version
  by Chad Greene. MIT license; see `external/cmocean/LICENSE`.
- `src/circ_vmrnd.m`: from the CircStat toolbox by Philipp Berens and
  Marc J. Velasco (2009).
