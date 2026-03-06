'''
    Written by Rachel Bertaud during PhD at Colroado School of Mines, 2026
    
    To analyze the transient behavior of the sprinkler, it is imperitive to know
    the natural frequency (omega) and the damping coefficent (gammma) of the system.
    By assuming that the tail of the measured angular position data of the sprinkler
    fits the ODE for a damped, hormoic oscillator, we have found a method to find
    those values, and use them to extract the torque signal of the system.
    
    This code:
        0. Read user defined inputs that control outputs
        1. Reads in angular position data in csv format and request fit time from user
        2. Estimates omega (natural frequency) and gamma (damping coefficent)
        3. Cleans noise from the data before and after forcing
        4. Uses Fourier tranforms to produce the torque signal for the data
        5. Does a forward ODE solve using the extracted torquue signal to verify the results

    All sections of this code are clearly labelled as step 0, 1, 2, 3, 4, or 5 for
    process trasparency.
'''

# SECTION ZERO - USER DEFINED INPUTS
###################################################################################################

# switches plots on (1) or off (0)
plot_switch = 1

# fit data after producing estimates of gamma dn omega (1) or no fit (0)
fit_switch = 1

# processes the data before forcing and after forcing to remove noise (1)
# or uses raw data (0)
proc_data_switch = 1

# define the spin direction of data to use
# reads from data file name, i.e. for "forward_500_trail1" put "forward" here
spin_dir = "forward"

# define reynolds number of data to use
# reads from data file name, i.e. for "forward_500_trail1" put "500" here 
re = "500"

# define trail number  of data to use
# reads from data file name, i.e. for "forward_500_trail1" put "1" here 
trial_num = 1

# define where data is stored on local machine
data_dir = "/Users/rachelbertaud/code/Sprinkler_Data/"


# DEFINE DEPENDENCIES AND UDFs
import os
import matplotlib.pyplot as plt
import numpy as np

# DATA READ FUNCS
from dataread_funcs import read_name, read_data, plot_data, fit_segments

# ESTIMATE FUNCS
from estimate_funcs import est_period, est_gamma, get_constants, fit_phi

# PROCESS FUNCS
from process_funcs import combine_data, remove_noise

# PLOT FUNCTIONS
from plot_funcs import plot_analytical, plot_franken, plot_phi_gen

from fft_funcs import torque_solver, phi_from_torque


# SECTION ONE - READ ANGULAR POSITION DATA AND REQUEST FIT LOCATION
###################################################################################################

os.chdir(data_dir)
fname = spin_dir + "_" + re + "_trial" + str(trial_num) + ".csv"
# fname = data_name

# define values from user define inputs
spin_switch, re, trial = read_name(fname)

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
period_est = est_period(peak_index, t_insert, t_peaks)

# estimates gamma from data
gamma_est = est_gamma(peak_index, t_peaks, y_peaks)
print("Estimate of gamma: ", gamma_est)

# estimate omega from period and gamma
omega_est = np.sqrt((period_est*period_est) + (gamma_est*gamma_est))
print("Estimate of omega: ", omega_est)

# gets constants for ODE using ours estimates for gamma and omega
# given fit t value
c1_est, c2_est = get_constants(gamma_est, omega_est, full_y[index])
print("Estimate of c1: ", c1_est)
print("Estimate of c2: ", c2_est)
# if user wants a fit...
if(fit_switch == 1):
    # fit the data
    gamma, omega, c1, c2 = fit_phi(t_seg, y_seg, gamma_est, omega_est, c1_est, c2_est, full_t[index])
else:
    # set end values as estimates
    gamma = gamma_est
    omega = omega_est
    c1 = c1_est
    c2 = c2_est

print("Gamma final: ", gamma)
print("Omega final: ", omega)
print("c1 final: ", c1)
print("c2 final: ", c2)
# define the analytical solution in terms of new found values
def phi_an(t):
    t0 = t_peaks[peak_index]
    wd = np.sqrt(omega**2 - gamma**2)
    return np.exp(-gamma * (t - t0)) * (c1 * np.cos(wd * (t - t0)) + c2 * np.sin(wd * (t - t0)))

error = np.linalg.norm(phi_an(full_t[index:insert_index]) - full_y[index:insert_index]) / np.linalg.norm(full_y[index:insert_index])
print("Error of Analytical Fit: ", error)

if(plot_switch == 1):
    plot_analytical(full_t, full_y, phi_an(full_t), index)


# SECTION THREE - CLEANS NOISE FROM DATA
###################################################################################################

if(proc_data_switch == 1):
    franken_y = remove_noise(full_t, full_y, threshold=1)

franken_t, franken_y, t_end = combine_data(full_t, full_y, index, phi_an)

if(plot_switch == 1):
    plot_franken(franken_t, franken_y, full_t, index, t_end)

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

phi_gen = phi_from_torque(N_f, franken_t, signal, gamma, omega)

if(plot_switch == 1):
    plot_phi_gen(full_t, full_y, phi_gen)
