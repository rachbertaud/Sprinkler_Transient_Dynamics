# Transient Dynamics Sprinkler — `main.py`

Written by Rachel Bertaud, Colorado School of Mines, 2026.

---

## What This Code Does

Extracts the **torque signal** of a sprinkler system from measured angular position data. It works by assuming the free-decay tail of the sprinkler's angular position fits a **damped harmonic oscillator ODE**, then uses that fit to identify the system's natural frequency (ω) and damping coefficient (γ). With those parameters known, a Fourier-based method reconstructs the full torque signal. A forward ODE solve is then run to verify the result.

---

## System Overview

```
Step 0 → Set user-defined inputs (switches, data selection)
Step 1 → Read CSV data, interactively pick fit window
Step 2 → Estimate ω (natural frequency) and γ (damping coefficient)
Step 3 → Clean noise from data (optional)
Step 4 → Extract torque signal via Fourier transform
Step 5 → Forward ODE solve to verify extracted torque signal
Step 6 → Save torque signal to CSV and plot results
```

---

## Dependencies

```
numpy
scipy
matplotlib
```

Install with:
```bash
pip install numpy scipy matplotlib
```

---

## Data Format

Input data is a `.csv` file with **two columns** (header row is skipped):

| Column 0 | Column 1 |
|----------|----------|
| time (t) | angular position (φ) |

**File naming convention:** `{direction}_{Re}_trial{N}.csv`

- `direction`: `forward` or `rev`
- `Re`: Reynolds number (e.g. `500`, `1000`)
- `N`: trial number (e.g. `1`, `2`)

Example: `forward_1000_trial1.csv`

---

## User-Defined Inputs (Step 0)

At the top of `main.py`, set these variables before running:

| Variable | Type | Description |
|---|---|---|
| `plot_switch` | `0` or `1` | `1` = show all plots at end, `0` = suppress all plots |
| `write_t` | `0` or `1` | `1` = save extracted torque signal to CSV, `0` = skip saving |
| `fit_switch` | `0` or `1` | `1` = refine ω/γ estimates with a curve fit, `0` = use raw peak estimates |
| `proc_data_switch` | `0` or `1` | `1` = clean noise from pre/post-forcing data, `0` = use raw data |
| `spin_dir` | string | Spin direction — `"forward"` or `"rev"` |
| `re` | string | Reynolds number matching the filename of interest — e.g. `"1000"` |
| `trial_num` | int | Trial number matching the filename of interest — e.g. `1` |
| `data_dir` | string | **Full path** to the directory where your CSV data files live |

---

## How to Run

```bash
python main.py
```

The script is **interactive** — it will open two plot windows in sequence, each with a slider:

1. **First window:** Drag the slider to select the **start time** for the free-decay fit (should be near a peak just before the sprinkler is released). Close the window when done.
2. **Second window:** Drag the slider to select the **end time** of the fit region. Close the window when done.

After both windows are closed, the script runs automatically and prints estimated/fitted values to the terminal.

---

## Outputs

**Terminal output:**
- Estimated γ, ω, and ODE constants (c1, c2) from peak analysis
- Fitted γ, ω, c1, c2 after curve fitting (if `fit_switch = 1`)
- Path of saved signal file (if `write_t = 1`)

**Saved file (if `write_t = 1`):**

Output is written to `data_dir` with the name `{spin_dir}_{Re}_trial{N}_signal.csv` — e.g. `forward_1000_trial1_signal.csv`. If a file with that name already exists, it is deleted and replaced. The file contains two columns covering the second half of the stitched signal (the physically meaningful forcing window):

| Column | Description |
|--------|-------------|
| `t` | time |
| `torque_signal` | real part of extracted torque signal |

**Plots (if `plot_switch = 1`):**
1. **Analytical solution** — raw data (scatter) overlaid with the damped oscillator fit (line); zoomed to the fit region, with the fit start point highlighted
2. **Franken signal** — the stitched signal used for the FFT, with colour-coded segments: prepended constant (blue), original data (pink), and analytical tail if `proc_data_switch = 1` (purple)
3. **Extracted torque signal** — real part of the torque signal over the forcing window
4. **ODE verification** — forward-solve φ(t) (line) overlaid on measured data (scatter)

---

## File Structure

```
Transient_Dynamics_Sprinkler/
├── main.py                          # Main code — run this
└── Main_Functions/
    ├── dataread_funcs.py            # CSV reading, peak finding, interactive slider plots
    ├── estimate_funcs.py            # ω and γ estimation and curve fitting
    ├── process_funcs.py             # Noise removal and signal stitching ("franken" signal)
    ├── fft_funcs.py                 # Fourier-based torque extraction and forward ODE solve
    └── plot_funcs.py                # All plotting functions
```

---

## Method Summary

1. **ω estimation:** Period is measured from peak-to-peak spacing in the free-decay tail. `ω_d = 2π / T_mean`, then `ω = sqrt(ω_d² + γ²)`.

2. **γ estimation:** Peak amplitudes decay as `e^(-γt)`. An exponential is fit to normalized peak amplitudes to extract γ.

3. **ODE fit (optional):** Both ω and γ are further refined by fitting the full free-decay segment to the analytical damped oscillator solution using `scipy.optimize.curve_fit`.

4. **Signal stitching ("Franken" signal):** A constant pre-forcing segment is prepended to the data. If the data is being processes, this includes the analytical free-decay solution is appended after the forcing ends. This creates a periodic-compatible signal for the FFT.

5. **Torque extraction (FFT):** The torque signal `f(t)` is extracted by deconvolving the ODE's transfer function `G(k)` in Fourier space:
   `F(k) = Y(k) / G(k)`, where `G(k) = -1 / (k² + 2iγk - ω²)`

6. **Verification (forward ODE):** The extracted torque signal is used as forcing when solving ODE (RK45) to reconstruct φ(t). 
