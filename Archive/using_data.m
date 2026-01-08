fs = 20;  % set your preferred font size here

% plotting settings
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', fs);
set(groot, 'defaultTextFontSize', fs);
set(groot, 'defaultLegendFontSize', 30);

% Optional: line width and marker size defaults
set(groot, 'defaultLineLineWidth', 3);
set(groot, 'defaultLineMarkerSize', 15);

%% Main Code

y_0 = 0; %IC for y
y_prime_0 = 0.00001; %IC for y'

peak_index = 55; %At what t index for peaks we want the inverse solve

index_full = 0;

N = length(full_X);

for i = 1:1:N
    if(full_X(i) == x_peaks(peak_index))
        index = i;
    end 
end




%%

close all
clf

% [c1, c2] = get_constant(gamma, w_0, x_peaks(index), phi_gen(i_peaks(index))); %Use solution at time t (user determines time t using index) where a peak occurs
% [c1, c2] = get_constant(gamma, w_0, x_peaks(peak_index), y_peaks(peak_index)); %Use solution at time t (user determines time t using index) where a peak occurs

%%
% Fitting for gamma and omega_0
% akes our generated solution, and fits the data to solve for omega and gamma.
t_seg = full_X(index:end);   %going moderately far into the system     
y_seg = full_y(index:end); 

%%
% ft = fittype(@(gamma, w0, x) ...
%     exp(-gamma.*x) .* ( ...
%         c1.*cos( sqrt(w0.^2 - gamma.^2).*x ) + ...
%         c2.*sin( sqrt(w0.^2 - gamma.^2).*x ) ), ...
%         'independent', 'x', ...
%         'coefficients', {'gamma','w0'});
% Above defines the fit we want for gamma and w0 using c_1 and c_2 we
% derived above

ft = fittype(@(gamma,w0,c1,c2,x) ...
    exp(-gamma.*x) .* ( ...
        c1.*cos( sqrt(w0.^2 - gamma.^2).*x ) + ...
        c2.*sin( sqrt(w0.^2 - gamma.^2).*x ) ), ...
    'independent','x', ...
    'coefficients', {'gamma','w0','c1','c2'});



opts = fitoptions(ft);
opts.StartPoint = [0.2, 10.0, 5.0, -5.0];    % [gamma w0 c1 c2] [0.2, 10.0, 5.0, -5.0]; 



mdl = fit(t_seg', y_seg', ft, opts);

% [fit_1, gof] = fit(t_seg', y_seg', ft, 'StartPoint', [gamma_guess, w_0_guess]); %Calling the actual fit, with user inputted guesses
%%
gamma = mdl.gamma;
w_0    = mdl.w0;
c1    = mdl.c1;
c2    = mdl.c2;


%%

% Makes function using derived stuff
phi_an = @(x) (exp(-gamma.*x)).*(c1.*cos(sqrt(w_0^2 - gamma^2).*x) + c2.*sin(sqrt(w_0^2 - gamma^2).*x)); %define the solution in terms of c_1 and c_2

%Plotting for visual validation
figure(1)
hold on
plot(full_X, full_y, 'm-')
plot(full_X, phi_an(full_X),'b')
plot(x_peaks(peak_index), y_peaks(peak_index), 'rx')
% plot(full_X(index), full_y(index), 'go')
ylim([-4,4])
legend('Data from Experiment', 'Analytical Fit', 'Point Used for Fit')
hold off

error = norm(phi_an(full_X(index:end)) - full_y(index:end))/norm(full_y(index:end));
fprintf("norm two error of data and analytical fit: %d \n", error)

%%
spot = full_X(end);
dx = mean(diff(full_X));
x_end = spot + dx:dx:(spot+(dx*(2048-length(full_X))));
x_neg = -(2048*dx):dx:0 - dx;
y_neg = zeros(1,2048);
franken_x = [x_neg, full_X, x_end];
franken_y = [y_neg, full_y, phi_an(x_end)];

figure(2)
hold on 
plot(franken_x, franken_y)
legend('Combined Data for FFT')
hold off
%% 

N = 2*2048;     %Domain definitions so that everything is the same size when we go into fourier space
L = franken_x(end)*2;
h = L/N; 
x = h*(1:N)'-L/2;  %Defining x (or t) domain for function


G_hat = @(k) -1.0./(k.^2 + 2.0*1i*gamma.*k - w_0^2);

kk = ((2*pi)/L)*[0:N/2-1 0 -N/2+1:-1]; 

signal = ifft(fft(franken_y)./G_hat(-kk));  

%%

for i = 1:1:N
    if(abs(franken_x(i) - 12.1234) < 5e-5)
        index = i;
    end 
end


%%

%Plots output
figure(3)
hold on
plot(franken_x, signal)
plot(franken_x(index), franken_y(index), 'ro')

legend('Inverse Solve')
title("Solving for Tau with Inverse FFT")
xlim([0,40])
hold off

figure(4)
hold on
plot(franken_x, franken_y)
plot(franken_x(index), franken_y(index), 'ro')
xlim([0,40])
legend('Combined Data')
hold off
