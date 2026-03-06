# for_fig — Running the Pipeline on MATLAB Figure Files

This is an older version of the code that reads data directly from saved MATLAB `.fig` files rather than CSVs. It runs a single trial at a time. The math is identical to `for_csv` — only the data loading is different.

> **Heads up:** this version is hardcoded around the specific `.fig` files in this directory. If you want to use it with your own data, you'll need to adjust the peak indices and possibly the figure loading logic to match your files.

# Getting started

Open MATLAB, navigate to this folder, and run `TSD.m`.

## Setting things up

Edit the **User Inputs** block at the top of `TSD.m`:

```matlab
plot_switch        = 1;   % show plots (1) or skip them (0)
spin_switch        = 0;   % 0 = reverse spin, 1 = forward spin
fit_switch         = 1;   % refine the decay fit with a nonlinear solver (recommended)
process_data_switch = 0;  % leave this at 0

peak_index = 49;  % which peak to use as the start of the free-decay fit
                  % 49 works well for reverse spin; try 60 for forward
```

The key thing to set is `peak_index`. It needs to point to a peak that falls after the forcing has stopped — somewhere in the free-decay tail. If you're getting bad fits, this is usually the first thing to adjust.

# Input files

The two `.fig` files should be in the same folder as `TSD.m`:

| File | What's in it |
|---|---|
| `reverse_800.fig` | Reverse-spin rotation data at Re ≈ 800 |
| `forward_800.fig` | Forward-spin rotation data at Re ≈ 800 |

Each figure has two embedded line objects: the full φ(t) time series and the pre-identified peaks. The code pulls those out directly — no peak detection needed.

# Output

This version doesn't write anything to disk — results are shown as plots (if `plot_switch = 1`). Two relative L2 errors are printed to the console: one for how well the analytical decay fits the free-decay tail, and one for how well the ODE reconstruction matches the full original signal.

# What the code does

1. **Read the figure** (`read_data`) — opens the `.fig` with `openfig`, digs out the two line objects, and hands back the full time series and peak arrays. The peak at `peak_index` sets the cut point between the forced and free-decay regions.

2. **Estimate frequency** (`estimate_Omega`) — works out the damped oscillation frequency Ω from the spacing between consecutive peaks (with outlier removal).

3. **Estimate damping** (`estimate_gamma`) — fits an exponential envelope to the same-sign peaks to get γ.

4. **Get the integration constants** (`get_constant`) — solves for c₁ and c₂ analytically from the peak value and the zero-derivative condition at that peak.

5. **Refine the fit** (`fit_phi`, optional) — runs a nonlinear least-squares solver over the free-decay tail to tighten up all four parameters: γ, ω₀, c₁, c₂. Turn this off with `fit_switch = 0` if you want the initial estimates only.

6. **Build the padded signal** (`combine_data`) — stitches together the full signal: zeros on the left, real data up to the cut point, and the analytical free-decay solution on the right. This gives the FFT a clean, uniform domain to work with.

7. **Recover the torque** (`torque_solver`) — in Fourier space the ODE becomes a simple division: τ̂(k) = φ̂(k) / G(k). Invert it to get τ(t).

8. **Verify the result** (`phi_from_torque`) — plugs the recovered torque back into the ODE with `ode45` and checks how well it reproduces the original φ(t). Small errors mean the recovery worked.


# Files in this folder

| File | What it does |
|---|---|
| `TSD.m` | Entry point — single-trial run |
| `read_data.m` | Opens the `.fig` and extracts the time series and peaks |
| `estimate_Omega.m` | Estimates damped frequency from peak spacing |
| `estimate_gamma.m` | Fits an exponential to the peak envelope to get γ |
| `get_constant.m` | Solves for c₁ and c₂ analytically |
| `fit_phi.m` | Nonlinear refinement of all four parameters |
| `combine_data.m` | Builds the zero-padded franken-signal for the FFT |
| `torque_solver.m` | FFT-based inverse solve — the core of the pipeline |
| `phi_from_torque.m` | ODE45 forward solve to verify the recovered torque |
| `plot_*.m` | Plotting helpers (controlled by `plot_switch`) |
| `plot_settings.m` | Global formatting: LaTeX interpreter, font sizes, line widths |
| `reverse_800.fig` | Input data — reverse spin at Re ≈ 800 |
| `forward_800.fig` | Input data — forward spin at Re ≈ 800 |
