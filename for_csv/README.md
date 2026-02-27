# Torque Signal Decomposition from Sprinkler Rotation Data

This code analyzes angular rotation data $\phi(t)$ of a sprinkler over time. It fits the data to extract parameters $\gamma$ (damping) and $\omega$ (angular frequency),
which are then used to recover the torque signal via Fourier transforms. This code outputs torque signals , whcih are saved as `.csv` files in a user-specified data directory.

## Usage

Edit the **User Inputs** section at the top of `main.m` before running:
```matlab
%% User Inputs

% Display plots (1 = on, 0 = off)
% Note: plots are not fully refined yet — recommended off unless you're debugging
plot_switch = 1;

% Fit data after initial parameter estimates (1 = fit, 0 = skip fit)
% Recommended: keep on (1) — it's fast and improves results marginally
fit_switch = 1;

% Process raw data (1 = process, 0 = skip)
process_data_switch = 0;

% Path to data directory
data_dir = "<path_to_data>";
```

### Data Configuration

Data files should follow the naming convention: `<spin_dir>_<Re>_trial<#>.csv`


```matlab
spin_dir   = "forward";  % Spin direction ("forward" or "reverse")
Re         = "500";      % Reynolds number
num_trials = 3;          % Number of trials to analyze
```
> **Note:** Target time for parameter fit must be AFTER forcing has ended.
% Tip: plot full_x vs. full_y to visually identify this value.
```matlab
% Target time for parameter fit — must be AFTER forcing has ended.
% Tip: plot full_x vs. full_y to visually identify this value.
t_target = 48;
```

## Output

Torque signals for each trial are written to `.csv` files in `data_dir`.
