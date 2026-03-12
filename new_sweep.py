"""
sweep.py - Runs the verf_of_main pipeline for all combinations of
           fit_switch and proc_data_switch, captures all error outputs,
           and writes them to sweep_results.txt.

The user picks t_target and t_insert once interactively at the start;
those values are then reused for every combination automatically.
"""

import io
import sys
import os
import traceback
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from itertools import product
import csv

sys.path.append(os.path.join(os.path.dirname(__file__), 'Verf_Functions'))

from construct_funcs import square_wave, error
from dataread_funcs import plot_data, fit_segments
from estimate_funcs import est_omega_d, est_gamma, get_constants, fit_phi
from process_funcs import combine_data, remove_noise, denoise_from_quiet_region
from fft_funcs import torque_solver, phi_from_torque, phi_from_torque_new
# PLOT FUNCTIONS
from plot_funcs import plot_analytical, plot_franken, plot_phi_gen, plot_torque


# ---- Fixed Parameters (same as verf_of_main.py) ----
gamma_exact = 0.548
omega_exact = 6.24
print("gamma exact: ", gamma_exact)
print("omega_exact", omega_exact)

res = 1
N_start     = int(2048)
tau_mag     = 50
start_sq    = 11
end_sq      = 23
count = 0

N_exact = 2048
t_exact = np.linspace(0, 60, N_exact)
tau_exact = square_wave(tau_mag, t_exact, start_sq, end_sq, 1)
#plt.plot(t_exact, tau_exact)
phi_exact = phi_from_torque_new(t_exact, tau_exact, gamma_exact, omega_exact)
t_int_exact = np.trapezoid(tau_exact, x=t_exact)
print("exact torque integral: ", t_int_exact)

plt.plot(t_exact, phi_exact, color='hotpink')
plt.show()

plt.plot(t_exact,tau_exact, color='hotpink')
plt.show()

fname = f"Results_N_{N_start}_Gam_{gamma_exact:.3f}_Omega_{omega_exact:.3f}.csv"

outpath = os.path.join(os.path.dirname(__file__), fname)

if os.path.exists(outpath):
    os.remove(outpath)

with open(outpath, 'w', newline='') as csvfile:
    wrtr = csv.writer(csvfile)
    wrtr.writerow(['Beginning Sweep!'])
    wrtr.writerow(['exact t', *t_exact])
    wrtr.writerow(['exact phi', *phi_exact])
    wrtr.writerow(['exact tau', *tau_exact])
    wrtr.writerow(['exact int', t_int_exact])


for res in range(1,2):
    # Build base data
    N = N_start * res
    # for downsampling tau

    indices = np.linspace(0, N_exact-1, N, dtype=int)
    t         = np.linspace(0, 60, N)

    # tau_copy = tau_exact.copy()
    # phi_copy = phi_exact.copy()
    tau = tau_exact[indices]
    phi = phi_exact[indices]
    # plt.plot(t, tau)
    # plt.show()
    
    dt = 60/N
    
    # t_check   = np.linspace(0, end_t, N_exact)
    # plt.plot(t_check,tau, '-o')
    # plt.plot(t, tau_exact, '-o')
    # plt.show()
    # 
    # plt.plot(t, phi)
    # plt.show()



    # # Show setup plots (same as verf_of_main.py section 1)
    # plt.plot(t, tau, '-o',color='pink')
    # plt.title("User Defined Square Wave")
    # plt.show()

    # plt.plot(t, phi_exact,'-o', color='pink')
    # plt.title("User Defined Phi Data")
    # plt.show()

    if(count == 0):
        # User picks t_target and t_insert once — reused for all combinations
        # t_target = plot_data(t_exact, phi_exact, 1)
        # t_insert = plot_data(t_exact, phi_exact, 0)
        t_target = 28.52
        t_insert = 36.02
        count += 1


    with open(outpath, 'a', newline='') as csvfile:
        wrtr = csv.writer(csvfile)
        wrtr.writerow([*['*'*2048]])
        wrtr.writerow(['STARTING VALUES'])
        wrtr.writerow(['dt', 'gamma', 'omega', 'tau mag.','tau start time', 'tau end time', 't target'])
        wrtr.writerow([dt, gamma_exact, omega_exact, tau_mag, start_sq, end_sq, t_target])
        wrtr.writerow([*['-'*2048]])

        

            
    print('-' * 60)

    header = (
        f"Sweep Starting!\n"
        f"dt = {dt:.4f} s,  NRES = {N}"
    )
    print(header)


    with open(outpath, 'a', newline='') as csvfile:
        wrt = csv.writer(csvfile)
        wrt.writerow(['t', *t])
        wrt.writerow(['phi', *phi])
        wrt.writerow(['tau', *tau])
        wrt.writerow([*['-'*2048]])


        for fit_switch, proc_data_switch in product([0, 1], [0, 1]):
            wrt.writerow(['fit switch', fit_switch, 'proc_data_switch', proc_data_switch])

            print('*' * 30)
            print('fit switch :', fit_switch, ', proc switch :', proc_data_switch)
            # Fresh data copies each iteration (remove_noise mutates in-place)
            full_t = t.copy()
            full_y = phi.copy()
            # plt.plot(full_t, full_y)
            # plt.show()

            # Peaks
            loc, _ = find_peaks(np.abs(full_y))
            t_peaks = full_t[loc]
            y_peaks = full_y[loc]

            # Segment selection (replaces fit_segments called after plot_data)
            index, peak_index, insert_index, t_seg, y_seg = fit_segments(
                full_t, full_y, t_peaks, t_target, t_insert
            )

            # Section 2 - Estimate omega and gamma
            omega_d   = est_omega_d(peak_index, t_insert, t_peaks)
            gamma_est = est_gamma(peak_index, t_peaks, y_peaks)
            #print(gamma_est)
            omega_est = np.sqrt(omega_d**2 + gamma_est**2)
            c1_est, c2_est = get_constants(gamma_est, omega_est, full_y[index])

            if fit_switch == 1:
                gamma, omega, c1, c2 = fit_phi(
                    t_seg, y_seg, gamma_est, omega_est, c1_est, c2_est, full_t[index]
                )
            else:
                gamma, omega, c1, c2 = gamma_est, omega_est, c1_est, c2_est

            print("gamma :", gamma)
            print("omega : ", omega)
            err_w = error(omega, omega_exact, 0, t)
            err_g = error(gamma, gamma_exact, 0 ,t)
            # Capture loop-local values in defaults so the closure is safe
            def phi_an(t_val, _g=gamma, _w=omega, _c1=c1, _c2=c2, _t0=t_peaks[peak_index]):
                wd = np.sqrt(_w**2 - _g**2)
                return np.exp(-_g * (t_val - _t0)) * (
                    _c1 * np.cos(wd * (t_val - _t0)) + _c2 * np.sin(wd * (t_val - _t0))
                )

            wrt.writerow(['phi_an', *phi_an(t[index:])])

            # plt.plot(t[index:], phi_an(t[index:]))
            # plt.plot(t[index:], phi_exact[index:])
            err_an = error(phi_an(t[index:]), phi[index:], 1, t[index:])
            #print(f"Error of analytical case with noise:     {err_an:.6f}")

            #print("------------------------------")



            # Section 3 - Process / combine data
            if proc_data_switch == 1:
                full_y = remove_noise(full_t, full_y, threshold=1)

            N_f = 2 * N_start * res
            N_f2 = N_start * res
            franken_t, franken_y, t_end = combine_data(full_t, full_y, index, phi_an, proc_data_switch, N_f, N_f2)

            # Section 4 - FFT / torque extraction

            L      = franken_t[-1] - franken_t[0]
            signal = torque_solver(N_f2, gamma, omega, L, franken_y)

        
            # plt.plot(t, signal[N_f2:],'-')
            # plt.plot(t, tau, color='pink')
            # plt.show()
            # plt.clf()


            wrt.writerow(['torque', *signal[N_f2:]])
            # print('sig', len(signal))
            # print('t', len(franken_t))
            # calculates torque signal integral
            torque_integral = np.trapezoid(signal[N_f2:], x=franken_t[N_f2:])
            err_int = error(torque_integral, t_int_exact, 0, t)

            print("torque int :", torque_integral)
           # print(f"Error of torque int:       {err_int:.6f}")

            # plt.plot(t, tau, color='purple')
           # print("------------------------------")

            err_sig = error(signal[N_f2:], tau, 1, t)

            # plt.plot(t, signal[N_f2:],'-')
            # plt.plot(t, tau, color='pink')
            # plt.show()

            # print("------------------------------")

            # # Section 5 - Forward ODE solve
            phi_gen = phi_from_torque(N_f2, franken_t, signal, gamma, omega)

            wrt.writerow(['generated phi', *phi_gen])

    
            err_for = error(phi_gen, phi, 1, t)
            # print(f"Error of forward case with noise:     {err_for:.6f}")

            wrt.writerow([ 'err omega', 'err gamma', 'err integral', 'err an/phi', 'err tau/syn', 'err for/phi'])
            wrt.writerow([err_w, err_g, err_int, err_an, err_sig, err_for])
            wrt.writerow([*['-'*2048]])

            new_t   = franken_t[N_f2:]
            new_y   = franken_y[N_f2:]
            new_sig = signal[N_f2:]

            plot_analytical(full_t, full_y, phi_an(full_t), index)
            plot_franken(franken_t, franken_y, full_t, index, t_end)
            plot_torque(franken_t, signal)
            plot_phi_gen(new_t, new_y, phi_gen)
            plt.show()

            full_t      = None
            full_y      = None
            franken_t   = None
            franken_y   = None
            signal      = None
            phi_gen     = None
            t_seg       = None
            y_seg       = None




        print(f"Results written to: {outpath}")
        print('-' * 60)


