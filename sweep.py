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
from process_funcs import combine_data, remove_noise
from fft_funcs import torque_solver, phi_from_torque, phi_from_torque_new


# ---- Fixed Parameters (same as verf_of_main.py) ----
gamma_exact = 0.345
omega_exact = 7.412
print("gamma exact: ", gamma_exact)
print("omega_exact", omega_exact)

res = 1
N_start     = int(256)
tau_mag     = 200
start_sq    = 11
end_sq      = 23
count = 0

fname = f"Results_Gam_{gamma_exact:.3f}_Omega_{omega_exact:.3f}.csv"

outpath = os.path.join(os.path.dirname(__file__), fname)


if os.path.exists(outpath):
    os.remove(outpath)


for res in range(16,20,2):
    # Build base data
    N = N_start * res
    t         = np.linspace(0, 60, N)
    tau_exact = square_wave(tau_mag, t, start_sq, end_sq, 0)
    phi_exact = phi_from_torque_new(t, tau_exact, gamma_exact, omega_exact)
    t_int_exact = (end_sq - start_sq)*tau_mag
    print("exact torque integral: ", t_int_exact)
    dt = 60/N

    tau = square_wave(tau_mag, t, start_sq, end_sq, 1)
    phi = phi_from_torque_new(t, tau, gamma_exact, omega_exact)


    # # Show setup plots (same as verf_of_main.py section 1)
    # plt.plot(t, tau, '-o',color='pink')
    # plt.title("User Defined Square Wave")
    # plt.show()

    # plt.plot(t, phi_exact,'-o', color='pink')
    # plt.title("User Defined Phi Data")
    # plt.show()

    if(count == 0):
        # User picks t_target and t_insert once — reused for all combinations
        t_target = plot_data(t, phi, 1)
        t_insert = plot_data(t, phi, 0)
        count += 1


    with open(outpath, 'w', newline='') as csvfile:
        wrtr = csv.writer(csvfile)
        wrtr.writerow(['STARTING VALUES'])
        wrtr.writerow(['dt', 'gamma', 'omega', 'tau mag.','tau start time', 'tau end time', 't target'])
        wrtr.writerow([dt, gamma_exact, omega_exact, tau_mag, start_sq, end_sq, t_target])
        wrtr.writerow([])

    savename =  f"Results_ttarg_{t_target:.4f}_dt_{dt}_TauM_{tau_mag}.csv"

    savepath = os.path.join(os.path.dirname(__file__), savename)

    if os.path.exists(savepath):
        os.remove(savepath)
            
    print('-' * 60)

    header = (
        f"Sweep Starting!\n"
        f"dt = {dt:.4f} s,  NRES = {N}"
    )
    print(header)


    with open(outpath, 'w', newline='') as csvfile:
        wrt = csv.writer(csvfile)
        wrt.writerow(['t', *t])
        wrt.writerow(['phi with noise', *phi])
        wrt.writerow(['phi exact', *phi_exact])
        wrt.writerow(['tau with noise', *tau])
        wrt.writerow(['tau exact', *tau_exact])

        for fit_switch, proc_data_switch in product([0, 1], [0, 1]):
            wrt.writerow([*['-'*60]])
            wrt.writerow(['fit switch', fit_switch, 'proc_data_switch', proc_data_switch])

            print('*' * 30)
            print('fit switch :', fit_switch, ', proc switch :', proc_data_switch)
            # Fresh data copies each iteration (remove_noise mutates in-place)
            full_t = t.copy()
            full_y = phi.copy()

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
        
            err_an_c = error(phi_an(t[index:]), phi_exact[index:], 1, t[index:])
            #print(f"Error of analytical case with noise:     {err_an:.6f}")

            err_an = error(phi_an(t[index:]), phi[index:], 1, t[index:])
            #print(f"Error of analytical case with no noise:  {err_an:.6f}")

            #print("------------------------------")



            # Section 3 - Process / combine data
            if proc_data_switch == 1:
                full_y = remove_noise(full_t, full_y, threshold=1)

            N_f = 2 * N_start * res
            franken_t, franken_y, t_end = combine_data(full_t, full_y, index, phi_an, proc_data_switch, N_f)

            # Section 4 - FFT / torque extraction
    
            L      = franken_t[-1] - franken_t[0]
            signal = torque_solver(N_f, gamma, omega, L, franken_y)

            wrt.writerow(['torque', *signal[N_f//2:]])

            # calculates torque signal integral
            torque_integral = np.trapezoid(signal, franken_t)
            err_int = error(torque_integral, t_int_exact, 0, t)

            print("torque int :", torque_integral)
           # print(f"Error of torque int:       {err_int:.6f}")

    
           # print("------------------------------")

            err_sig = error(signal[N_f // 2:], tau, 1, t)
            #print(f"Error of signal with noise:           {err_sig:.6f}")

            err_sig_c = error(signal[N_f // 2:], tau_exact, 1, t)
            #print(f"Error of signal with no noise:        {err_sig:.6f}")

    
            # print("------------------------------")

            # # Section 5 - Forward ODE solve
            phi_gen = phi_from_torque(N_f, franken_t, signal, gamma, omega)

            wrt.writerow(['generated phi', *phi_gen[N_f//2:]])
    
            err_for = error(phi_gen, phi, 1, t)
            # print(f"Error of forward case with noise:     {err_for:.6f}")

            err_for_c = error(phi_gen, phi_exact, 1, t)
            # print(f"Error of forward case with no noise:  {err_for:.6f}")

    
            with open(outpath, 'a', newline='') as  csvfile:
                wrtr = csv.writer(csvfile)
                wrtr.writerow(['fit switch', 'proc switch', 'err omega', 'err gamma', 'err an/noise', 'error an/clean', 'err tau/noise', 'err tau/clean', 'err for/noise', 'err for/clean'])
                wrtr.writerow([fit_switch, proc_data_switch, err_w, err_g, err_int, err_an, err_an_c, err_sig, err_sig_c, err_for, err_for_c ])



        print(f"Results written to: {outpath} and {savepath}")
        print('-' * 60)


