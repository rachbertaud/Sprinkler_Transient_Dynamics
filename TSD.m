clear all 


%% User Inputs
plot_switch = 1; %switches on plots (1) or no plots (0)
spin_switch = 0; %reads reverse (0) or forward (1) data
fit_switch = 1; %fits data using estimates (1) or no fit after estimates (0)
process_data_switch = 0; %process data (1) or not (0) (no option to not process data right now)



%% Main Code

if(spin_switch == 0)
    peak_index = 49; %At what t index for peaks we want the inverse solve
else
    peak_index = 60;
end 
% 49 is good for reverse, 60 is good for forward

%reads data from figure and outputs:
% full data (full_x, full_y) 
% peak data (x_peaks, y_peaks)
% size of full data (N) 
% the index (index) where the peak index of interest (peak_index) occurs in full data set
% the full data set at and after the peak index of interest (x_seg, y_seg) 
[full_x, full_y, x_peaks, y_peaks, N, index, x_seg, y_seg] = read_data(spin_switch, peak_index);

%Estimates period of function (Omega_est) after peak index using x_peaks 
Omega_est = estimate_Omega(peak_index, x_peaks);

%estimate gamma (gamma) with exponential curve using peaks after peak index
gamma = estimate_gamma(peak_index, x_peaks, y_peaks);

%estimate omega_0 (w_0) using period (Omega_est) and gamma (gamma) est.
w_0 = sqrt(Omega_est^2 + gamma^2);

%Solve for constants using function below and our estimate gamma and w_0
[c1, c2] = get_constant(gamma, w_0, x_peaks(peak_index), y_peaks(peak_index)); %Use solution at time t (user determines time t using index) where a peak occurs

if(fit_switch == 1)
    [gamma, w_0, c1, c2] = fit_phi(x_seg, y_seg, gamma, w_0, c1, c2); %runs fit
end 

% Makes function using derived stuff
phi_an = @(x) (exp(-gamma.*x)).*(c1.*cos(sqrt(w_0^2 - gamma^2).*x) + c2.*sin(sqrt(w_0^2 - gamma^2).*x)); %define the solution in terms of c_1 and c_2

%Prints relative error between analytical solution and data
error = norm(phi_an(full_x(index:end)) - full_y(index:end))/norm(full_y(index:end));
fprintf("norm two error between data and analytical fit: %.5f \n", error)

[franken_x, franken_y] = combine_data(full_x, full_y, index, phi_an);

if(process_data_switch == 1)
    [franken_x, franken_y] = process_data(franken_x, franken_y, spin_switch); %creates processed version of full_x and full_y, franken_x and franken_y respectively
end 


N_f = 2*2048;     %Domain definitions so that everything is the nice size when we go into fourier space
L = franken_x(end)*2;%Determines length of domain
h = L/N_f; %determines spacing of dx to get number of points we want on domain
x = h*(1:N_f)'-L/2;  %Defining x (or t) domain for function

[signal] = torque_solver(N_f, gamma, w_0, L, franken_y); %produces torque using FFT to IFFT

torque_intergal = trapz(franken_x, signal);
fprintf("Intergral of Torque Signal is approximately %.2f \n", torque_intergal)



%%


fprintf("Cummulative Intergral of Torque Signal is approximately %.2f \n", torque_intergal)

% Going backwards
phi_gen = phi_from_torque(N_f, full_x, signal, gamma, w_0); %generates phi from analytical torque using ODE45

error1 = norm(phi_gen(1:index) - full_y(1:index))/norm(full_y(1:index));
fprintf("norm two error of original data and output using solved torque: %.5f \n", error1)



if(plot_switch == 1)
    plot_settings(plot_switch)
    plot_analytical_solution(full_x, full_y, peak_index, x_peaks, y_peaks, phi_an, spin_switch) %plots analytical solution
    plot_processed_data(franken_x, franken_y) %plots the processed data 
    plot_torque(franken_x, signal, N_f) %plots data from torque
    plot_phi_from_torque(full_x, full_y, phi_gen, spin_switch) %plots generated phi
end





