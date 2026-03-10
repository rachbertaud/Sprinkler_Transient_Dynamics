'''
    Written by Rachel Bertaud during PhD at Colroado School of Mines, 2026

    This code assumes we know what omega and gamma are, then goes through the process of main
    to verify our code.
    
    This code:
        0. Creates angular position data that we produce from a noisy square wave and known
            values of omega and gamma.
        1. Reads in angular position data in csv format and request fit time from user
        2. Estimates omega (natural frequency) and gamma (damping coefficent)
        3. Cleans noise from the data before and after forcing
        4. Uses Fourier tranforms to produce the torque signal for the data
        5. Does a forward ODE solve using the extracted torque signal
        6. Compares the results of the code to the known values.
        
    All sections of this code are clearly labelled as above for
    process trasparency.
'''

# SECTION ZERO - USER DEFINED INPUTS
###################################################################################################

# switches plots on (1) or off (0)
plot_switch = 1

# saves extracted torque signal (1) or does not (0)
write_t = 1

# fit data after producing estimates of gamma dn omega (1) or no fit (0)
fit_switch = 1

# processes the data before forcing and after forcing to remove noise (1)
# or uses raw data (0)
proc_data_switch = 0

gamma = 0.345
omega = 7.412
N = int(2048)
tau_mag = 100


# DEFINE DEPENDENCIES AND UDFs
import os
import matplotlib.pyplot as plt
import numpy as np
import sys

# enters path where the function files are 
sys.path.append(os.path.join(os.path.dirname(__file__), 'Verf_Functions'))

# MAKE TAU AND PHI FUNCS
from construct_funcs import square_wave

# ESTIMATE FUNCS
from estimate_funcs import est_omega_d, est_gamma, get_constants, fit_phi

# PROCESS FUNCS
from process_funcs import combine_data, remove_noise

# PLOT FUNCTIONS
from plot_funcs import plot_analytical, plot_franken, plot_phi_gen, plot_torque

from fft_funcs import torque_solver, phi_from_torque


# SECTION ONE - READ ANGULAR POSITION DATA AND REQUEST FIT LOCATION
###################################################################################################

if(fit_switch == 1):
    print("Data is being fit!")
if(proc_data_switch == 1):
    print("Data is being processed!")

t = np.linspace(0,N,N)*0.033

tau = square_wave(tau_mag, t, 11, 23)

plt.plot(t, tau)
plt.show()



# read data for x,y data, get size of data, and find peaks of data
full_t, full_y, N, t_peaks, y_peaks = read_data(fname)

#plt.plot(full_t, full_y)
#plt.scatter(t_peaks, y_peaks, s=200, facecolors='none',
#               edgecolors='limegreen', linewidths=2.5)
#plt.show()

# lets user look at plot and define fitting target
t_target = plot_data(full_t, full_y, 1)

# lets user look at plot and define where we will input the data
t_insert = plot_data(full_t, full_y, 0)

# returns segment of t and y after UD time for fitting target
index, peak_index, insert_index, t_seg, y_seg = fit_segments(full_t, full_y, t_peaks, t_target, t_insert)


# SECTION TWO - ESTIMATE OMEGA AND GAMMA
###################################################################################################

# estimates period from data
omega_d = est_omega_d(peak_index, t_insert, t_peaks)

# estimates gamma from data
gamma_est = est_gamma(peak_index, t_peaks, y_peaks)


# estimate omega from period and gamma
omega_est = np.sqrt((omega_d*omega_d) + (gamma_est*gamma_est))

# gets constants for ODE using ours estimates for gamma and omega
# given fit t value
c1_est, c2_est = get_constants(gamma_est, omega_est, full_y[index])

print("------------------------------")
print("-----ESTIMATES FROM PEAKS-----")
print("------------------------------")
print("Estimate of gamma: ", gamma_est)
print("Estimate of omega: ", omega_est)
print("Estimate of c1: ", c1_est)
print("Estimate of c2: ", c2_est)

# if user wants a fit...
if(fit_switch == 1):
    # fit the data
    gamma, omega, c1, c2 = fit_phi(t_seg, y_seg, gamma_est, omega_est, c1_est, c2_est, full_t[index])
    print("------------------------------")
    print("-------VALUES AFTER FIT-------")
    print("------------------------------")
    print("Gamma: ", gamma)
    print("Omega: ", omega)
    print("c1: ", c1)
    print("c2: ", c2)
    print("------------------------------")
    

else:
    # set end values as estimates
    gamma = gamma_est
    omega = omega_est
    c1 = c1_est
    c2 = c2_est
    
# define the analytical solution in terms of new found values
def phi_an(t):
    t0 = t_peaks[peak_index]
    wd = np.sqrt(omega**2 - gamma**2)
    return np.exp(-gamma * (t - t0)) * (c1 * np.cos(wd * (t - t0)) + c2 * np.sin(wd * (t - t0)))

#error = np.linalg.norm(phi_an(full_t[index:insert_index]) - full_y[index:insert_index]) / np.linalg.norm(full_y[index:insert_index])
error = np.linalg.norm(phi_an(full_t[index:]) - full_y[index:], 2) / np.linalg.norm(full_y[index:], 2)

error_full = np.abs((phi_an(full_t[index:]) - full_y[index:])) / (np.abs(full_y[index:]) + 1e-4)

#plt.plot(full_t[index:], error_full, color='purple')
#plt.title("Relative Absolute Error on Tail")
#plt.show()


#print("Error of Analytical Fit: ", error)



# SECTION THREE - CLEANS NOISE FROM DATA
###################################################################################################

if(proc_data_switch == 1):
    full_y = remove_noise(full_t, full_y, threshold=1)

franken_t, franken_y, t_end = combine_data(full_t, full_y, index, phi_an, proc_data_switch)


# SECTION FOUR - FOURIER TRANSFORM
###################################################################################################

N_f = 2*2048
L = franken_t[-1] - franken_t[0]
h = L/N_f

# outputs extracted torque signal
signal = torque_solver(N_f, gamma, omega, L, franken_y)

# calculates torque signal integral
torque_integral = np.trapezoid(signal, franken_t)

# going forwards
new_t   = franken_t[N_f//2:]
new_y   = franken_y[N_f//2:]
new_sig = signal[N_f//2:]

# SECTION FIVE - DO FORWARD PROBLEM WITH TORQUE SIGNAL
###################################################################################################

phi_gen = phi_from_torque(N_f, franken_t, signal, gamma, omega)

#error_gen = np.linalg.norm(phi_gen[:index] - new_y[:index], 2) / np.linalg.norm(new_y[:index], 2)
#print("Error of Transient Phase Recovery: ", error_gen)

#error_two = np.linalg.norm(phi_gen - new_y, 2) / np.linalg.norm(new_y, 2)
#print("Error for full phi recovery :", error_two)

#wtf = np.abs(phi_gen[:index] - new_y[:index])/(np.abs(new_y[:index]) + 1e-4)

#plt.plot(new_t[:index], wtf)
#plt.show()


# SECTION SIX - SAVE SIGNAL AND PLOT RESULTS
###################################################################################################

if(write_t == 1):
    out_fname = data_name + "_signal.csv"
    out_path = os.path.join(data_dir, out_fname)
    if os.path.exists(out_path):
        os.remove(out_path)
        print("Existing file removed: ", out_path)
    np.savetxt(out_path, np.column_stack([franken_t[N_f//2:], np.real(signal[N_f//2:])]),
           delimiter=',', header='t,torque_signal', comments='')
    print("Signal saved to: ", out_path)

if(plot_switch == 1):
    plot_analytical(full_t, full_y, phi_an(full_t), index)
    plot_franken(franken_t, franken_y, full_t, index, t_end)
    plot_torque(franken_t, signal)
    plot_phi_gen(new_t, new_y, phi_gen)
    plt.show()
