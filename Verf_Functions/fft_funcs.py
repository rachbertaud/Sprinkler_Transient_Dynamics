import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

def torque_solver(N_f, gamma, w_0, L, franken_y):
    def G_hat(k):
        return -1.0 / (k**2 + 2.0*1j*gamma*k - w_0**2)
    
    kk = (2*np.pi / L) * np.concatenate([
        np.arange(0, N_f//2),
        [0],
        np.arange(-N_f//2 + 1, 0)
    ])
    
    signal = np.fft.ifft(np.fft.fft(franken_y) / G_hat(-kk))
    
    return signal

def phi_from_torque_new(t_seg, signal, gamma, w_0):
    # interpolant for the forcing
    para = PchipInterpolator(t_seg, signal, extrapolate=False)

    def dydt(t, y):
        p = para(t) if t_seg[0] <= t <= t_seg[-1] else 0.0
        return [y[1], -(2*gamma)*y[1] - (w_0**2)*y[0] + p]

    y0 = [0.0, 0.0]
    sol = solve_ivp(dydt, [t_seg[0], t_seg[-1]], y0, t_eval=t_seg, method='RK45')

    phi_gen = sol.y[0]
    return phi_gen

def phi_from_torque(N_f, franken_t, signal, gamma, w_0):
    signal_chunk = np.real(signal[N_f//2:])  # signal from 0 to end of full_t

    t_seg = franken_t[N_f//2:]

    # interpolant for the forcing
    para = PchipInterpolator(t_seg, signal_chunk, extrapolate=False)

    def dydt(t, y):
        p = para(t) if t_seg[0] <= t <= t_seg[-1] else 0.0
        return [y[1], -(2*gamma)*y[1] - (w_0**2)*y[0] + p]

    y0 = [0.0, 0.0]
    sol = solve_ivp(dydt, [t_seg[0], t_seg[-1]], y0, t_eval=t_seg, method='RK45')

    phi_gen = sol.y[0]
    return phi_gen
