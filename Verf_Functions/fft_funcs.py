import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator, interp1d


def torque_solver(N_f2, gamma, w_0, L, franken_y):
    def G_hat(k):
        return -1.0 / (k**2 + 2.0*1j*gamma*k - w_0**2)
    
    kk = (2*np.pi / L) * np.concatenate([
        np.arange(0, N_f2),
        [0],
        np.arange(-N_f2 + 1, 0)
    ])
    
    signal = np.real(np.fft.ifft(np.fft.fft(franken_y) / G_hat(-kk)))
    
    return signal

# def phi_from_torque_new(t_seg, signal, gamma, w_0):
#     # interpolant for the forcing
#     para = PchipInterpolator(t_seg, signal, extrapolate=False)

#     def dydt(t, y):
#         p = para(t) if t_seg[0] <= t <= t_seg[-1] else 0.0
#         return [y[1], -(2*gamma)*y[1] - (w_0**2)*y[0] + p]

#     y0 = [0.0, 0.0]
#     sol = solve_ivp(dydt, [t_seg[0], t_seg[-1]], y0, t_eval=t_seg, method='RK45')

#     phi_gen = sol.y[0]
#     return phi_gen
#
def phi_from_torque_new(t_seg, signal, gamma, w_0):
    para = interp1d(t_seg, signal, kind='nearest', fill_value=0.0, bounds_error=False)
    def dydt(t, y):
        p = float(para(t))
        return [y[1], -(2*gamma)*y[1] - (w_0**2)*y[0] + p]
    y0 = [0.0, 0.0]
    sol = solve_ivp(dydt, [t_seg[0], t_seg[-1]], y0, t_eval=t_seg, method='RK45', rtol=1e-8, atol=1e-10, max_step=np.min(np.diff(t_seg)))

    # print(f"sol.success: {sol.success}")
    # print(f"sol.message: {sol.message}")
    # print(f"sol.y shape: {sol.y.shape}")
    # print(f"signal range: {signal.min():.4f} to {signal.max():.4f}")
    # print(sol.y[0].max())
    # print(f"t_seg range: {t_seg[0]:.4f} to {t_seg[-1]:.4f}")
    return sol.y[0]

def phi_from_torque(N_f2, franken_t, signal, gamma, w_0):
    signal_chunk = signal[N_f2:]  # signal from 0 to end of full_t

    t_seg = franken_t[N_f2:]

    # interpolant for the forcing
    para = PchipInterpolator(t_seg, signal_chunk, extrapolate=False)

    def dydt(t, y):
        p = para(t) if t_seg[0] <= t <= t_seg[-1] else 0.0
        return [y[1], -(2*gamma)*y[1] - (w_0**2)*y[0] + p]

    y0 = [0.0, 0.0]
    sol = solve_ivp(dydt, [t_seg[0], t_seg[-1]], y0, t_eval=t_seg, method='RK45')

    phi_gen = sol.y[0]
    return phi_gen
