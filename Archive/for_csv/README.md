# for_csv — Running the Pipeline on CSV Data

This is the main working version of the code. It reads rotation data from `.csv` files, processes each trial one by one, and saves the recovered torque signal back to disk when it's done.

# Getting started

Open MATLAB, navigate to this folder, and run `main.m`. Make sure you set up as follows in the next section! Everything else is handled automatically :).

> Note: This code and plots are only written for three trails. Edits will be needed to plots if more than three plots are wanted to keep the colors/transparency consistent.

## Setting things up

Before you run it, fill in the **User Inputs** block at the top of `main.m`:

```matlab
plot_switch        = 1;   % show plots while it runs (set to 0 to skip)
fit_switch         = 1;   % refine the decay fit with a nonlinear solver (recommended)
process_data_switch = 0;  % leave this at 0

data_dir   = "/path/to/your/data/";  % where your CSV files live
spin_dir   = "forward";              % "forward" or "reverse"
Re         = "500";                  % Reynolds number label for your data
num_trials = 3;                      % how many trial files to process
t_target   = 48;                     % a time (in seconds) that falls after forcing stops
```

The most important one to get right is `t_target`. It needs to land somewhere in the free-decay region of the signal — after the sprinkler has been let go and is coasting to a stop on its own. If you're not sure where that is, just plot your raw data first and look for where the forcing ends.

## What your data files need to look like

Files should be named like this:

```
forward_500_trial1.csv
forward_500_trial2.csv
forward_500_trial3.csv
```

The pattern is `<spin_dir>_<Re>_trial<N>.csv`. Each file needs two columns: time `t` and angular position `φ`. That's it.

# What you get out

For each trial, the recovered torque is saved as a CSV:

```
<data_dir>/<spin_dir>_<Re>_trial<N>_torque_signal.csv
```

The file has two columns, `t` and `phi` (the torque values), covering the second half of the padded domain — the part that corresponds to your real measurement.

The console also prints a summary for each trial that looks like this:

```
--------------------------------------------------
  Trial 1  |  Forward  |  Re = 500
--------------------------------------------------
  Fit error    (phi_an vs. data)  :  0.1234
  Torque integral                 :  12.3456
  Recon. error  (ODE45 vs. data)  :  0.2345
--------------------------------------------------
```

Two errors to keep an eye on: how well the analytical decay fits the tail of your data, and how well the ODE reconstruction matches the full signal. Both should be small.

# The plots

Set `plot_switch = 1` and you'll get five figures. All of them use the same visual style: each trial drawn faint in its own color, with a bold mean line on top.

| Color | Style | What it is |
|---|---|---|
| Hot pink | solid | Trial 1 |
| Soft purple | dashed | Trial 2 |
| Periwinkle blue | dash-dot | Trial 3 |
| Deep violet | solid (bold) | Mean across trials |

| Figure | What it shows |
|---|---|
| 1 | Raw φ(t) data and the analytical free-decay fit for each trial. A dotted vertical line marks where each trial's fit starts; the mean start time is shown as a dashed line. |
| 2 | The recovered torque signal τ(t). |
| 3 | The full padded signal that goes into the FFT. The zero-padded region on the left is shaded light purple-grey so you can see the boundary clearly. |
| 4 | The ODE45 forward reconstruction using the recovered torque, plotted against the original data. This is our sanity check. |
| 5 | The cumulative integral of τ(t) over time. A datatip marks the peak of the mean curve. |


# What the code actually does

For each trial, the code works through these steps:

1. **Read and find peaks** (`read_data`) — loads the CSV, finds all the local maxima of |φ(t)|, and picks the peak closest to your `t_target` as the handoff point between the forced and free-decay regions.

2. **Estimate frequency** (`estimate_Omega`) — figures out the damped oscillation frequency Ω by looking at the spacing between consecutive peaks. Outliers are removed before averaging.

3. **Estimate damping** (`estimate_gamma`) — fits an exponential envelope to the same-sign peaks to get the damping coefficient γ.

4. **Get the integration constants** (`get_constant`) — uses the peak value and the fact that the derivative is zero at a peak to solve for c₁ and c₂ analytically.

5. **Refine the fit** (`fit_phi`, optional) — runs a nonlinear least-squares solver over the whole free-decay tail to sharpen up all four parameters: γ, ω₀, c₁, c₂. Skip this with `fit_switch = 0` if you want the initial estimates only.

6. **Build the padded signal** (`combine_data`) — stitches together the full signal: zeros on the left (before the measurement starts), real data up to the cut point, and the analytical free-decay solution on the right. This gives the FFT a clean, consistent-length domain to work with.

7. **Recover the torque** (`torque_solver`) — in Fourier space, the ODE is just a division: τ̂(k) = φ̂(k) / G(k). Invert that and you have τ(t).

8. **Verify the result** (`phi_from_torque`) — plugs the recovered torque back into the ODE with `ode45` and checks how closely it reproduces the original φ(t). If the errors are small, the recovery worked.


# Files in this folder

| File | What it does |
|---|---|
| `main.m` | Entry point — loops over all trials |
| `read_name.m` | Parses spin direction, Re, and trial number out of a filename |
| `read_data.m` | Reads the CSV, finds peaks, cuts the data at `t_target` |
| `estimate_Omega.m` | Estimates the damped oscillation frequency from peak spacing |
| `estimate_gamma.m` | Fits an exponential to the peak envelope to get γ |
| `get_constant.m` | Solves for c₁ and c₂ analytically from the peak condition |
| `fit_phi.m` | Nonlinear refinement of all four parameters |
| `combine_data.m` | Builds the zero-padded franken-signal for the FFT |
| `torque_solver.m` | FFT-based inverse solve — the core of the pipeline |
| `phi_from_torque.m` | ODE45 forward solve to verify the recovered torque |
| `plot_analytical_solution.m` | Figure 1 — free-decay fits across trials |
| `plot_torque.m` | Figures 2 & 5 — recovered torque and its cumulative integral |
| `plot_processed_data.m` | Figure 3 — the padded FFT input |
| `plot_phi_from_torque.m` | Figure 4 — ODE45 reconstruction vs. original data |
| `plot_settings.m` | Global formatting: LaTeX interpreter, font sizes, line widths |
