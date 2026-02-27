# Transient Dynamics of Our Sprinkler

This project recovers the torque signal driving a sprinkler's rotation from measured angular position data φ(t). 

## What's in this repo

```
for_csv/       the main working version — reads data from .csv files, loops over trials
for_fig/       an older version that reads data from saved MATLAB .fig files
Archive/       even older scripts and figures from earlier development
Presentations/ slides and exported plots
```

**Both `for_csv` and `for_fig` run the same math. The only difference is how they load data — see each folder's README for the details. If you're starting fresh, use `for_csv`.**

## What the code does

The sprinkler is modeled as a damped harmonic oscillator. When the torque is present, it's the forcing term on the right:

$$\ddot{\phi} + 2\gamma\dot{\phi} + \omega^2 \phi = \tau(t)$$

where γ is the damping coefficient, ω is the natural frequency, and τ(t) is the torque we want to recover. 

### Step 1 — Extract γ and ω from the free-decay tail

Once the forcing stops, the sprinkler coasts freely and the ODE becomes homogeneous:

$$\ddot{\phi} + 2\gamma\dot{\phi} + \omega^2 \phi = 0$$

The tail of the signal is treated as a clean damped sinusoid. We can fit it to extract γ and ω.

### Step 2 — Recover τ(t) using Fourier transforms

With γ and ω, we go back to the full forced ODE and solve for τ(t) spectrally.

The ODE can be written as a convolution in the time domain:

$$\phi(t) = \int G(t - s)\,\tau(s)\,ds$$

where $G$ is the Green's function of the system. In Fourier space, convolution becomes multiplication:

$$\hat{\phi} = \hat{G} \cdot \hat{\tau}$$

To find $\hat{G}$, we take the Fourier transform of the ODE directly (where $k$ is the frequency variable):

$$-k^2\hat{\phi} - 2i\gamma k\hat{\phi} + \omega^2\hat{\phi} = \hat{\tau}$$

$$(-k^2 - 2i\gamma k + \omega^2)\hat{\phi} = \hat{\tau}$$

$$\hat{\phi} = \left(\frac{-1}{k^2 + 2i\gamma k - \omega^2}\right)\hat{\tau} = \hat{G}\cdot\hat{\tau}$$

So the Green's function in Fourier space is:

$$\hat{G} = \frac{-1}{k^2 + 2i\gamma k - \omega^2}$$

Since we know $\hat{G}$ and we can compute $\hat{\phi}$ from our data, we can solve directly for the torque:

$$\hat{\tau} = \frac{\hat{\phi}}{\hat{G}}$$

$$\tau = \text{IFFT}\!\left(\frac{\text{FFT}(\phi)}{\hat{G}}\right)$$

### Step 3 — Verify

To check the result, we plug τ(t) back into the ODE and solve it forward with `ode45`. If the reconstructed φ(t) matches the original data, the recovery worked.


