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

peak_index = 55; %At what t index for peaks we want the inverse solve

N = length(full_X);

%determine which index of the full data set our peak index for inverse solve occurs at
for i = 1:1:N
    if(full_X(i) == x_peaks(peak_index))
        index = i;
    end 
end


close all
clf

% Pulls out data after our cutting index for inverse solve
t_seg = full_X(index:end);        
y_seg = full_y(index:end); 

%% we want to fit for gamma and w0, but can guess at what b is for the fit,
%so here is where we will make an educated guess at b
aoi = x_peaks(peak_index:end);
end_peaks = length(aoi);
period_est = zeros(1, end_peaks - 2);
j = 1;
for i = peak_index:1:(length(x_peaks) - 2)
    period_est(j) = x_peaks(i + 2) - x_peaks(i); %find out the period after our picked solve index
    j = j + 1;
end

%% Can we guess for gamma too?? maybe.
x_gamma = [];
y_gamma = [];
for i = peak_index:2:(length(x_peaks))
    x_gamma = [x_gamma x_peaks(i)] ;
    y_gamma = [y_gamma y_peaks(i)] ;
end

x_gamma = x_gamma - x_gamma(1); %start from 0 might help?
hold on 
plot(x_peaks, y_peaks)
plot(x_gamma, y_gamma)
hold off

% can we normalize y_gamma?
y_gamma = y_gamma./y_gamma(1);


exp_fit_type = fittype("exp(-gamma*x)",...
    dependent='y', independent='x',...
    coefficients='gamma');
exp_fit = fit(x_gamma', y_gamma', exp_fit_type)

%%
b_est = (2*pi)/mean(period_est); %estimate the period using the mean of the period

t_0 = full_X(index);


%below suggests fit form for gamma, w0, c1, and c2
myfittype = fittype("exp(-gamma*x)*c1*cos( b*x ) + exp(-gamma*x)*c2*sin( b*x )", ... %note here that b = sqrt(w0^2 + gamma^2)
    dependent='y', independent='x',...
    coefficients=["gamma" "b" "c1" "c2"]);



%%
% Defines some starting guesses for gamma, b, c1, and c2
opts = fitoptions(myfittype);
opts.StartPoint = [0.05, b_est, 0, 0];    % [gamma w0 c1 c2] [0.2, 10.0, 5.0, -5.0]; 

% actually runs the fit
fit_guess = fit(t_seg', y_seg', myfittype, opts)

%%
plot(t_seg, y_seg)
%%
%pulling out guesses for gamma, w_0, c1, and c2
gamma = fit_guess.gamma;
b    = fit_guess.b;
c1    = fit_guess.c1;
c2    = fit_guess.c2;

w_0 = sqrt(b^2 + gamma^2);


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
