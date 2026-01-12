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

clc
clear all

%%Extract Data

fig = openfig('reverse_800.fig'); %open figure containing data from NYU 

%
dataObjs = findobj(fig, 'Type', 'line'); %find the line of data

xData = get(dataObjs, 'XData'); %get all x data from the line
yData = get(dataObjs, 'YData'); %get all y data from the line


%pull out peaks data (orange circles)
x_peaks = xData{1}; 
y_peaks = yData{1};

%pull out the phi data (blue curve)
full_X = xData{2};
full_y = yData{2};

clf 
close all
%%

% peak index will define what peak we are choosing for the solve of coeff.
peak_index = 49; 

N = length(full_X); %Pull out full length of data!

%determine which index of the full data set our peak index for inverse solve occurs at
for i = 1:1:N
    if(full_X(i) == x_peaks(peak_index))
        index = i; %index is the index point in the full data set that the peak defined by "peak_data" exists at
    end 
end

t_seg = full_X(index:end); %getting all x data after the peak of interest       
y_seg = full_y(index:end); %getting all y data after the peak of interest 

%% Fitting for initial guesses of gamma and w0

aoi = x_peaks(peak_index:end); % this is the peak area of interest, since this is all the points at and after our determined cutting point
end_peaks = length(aoi); %determine the length of the area of interest
period_est = zeros(1, end_peaks - 2); %define zero vector to fill with our period estimate

j = 1; %indexing variable
for i = peak_index:1:(length(x_peaks) - 2) %for all the points after our area of interest
    period_est(j) = x_peaks(i + 2) - x_peaks(i); %find out the period after our picked solve index
    j = j + 1;
end

b_est = (2*pi)/mean(period_est); %estimate the period using the mean of the period data found above


x_gamma = []; %to fill with gamma data
y_gamma = []; %to fill with gamma data
%we know we are starting at a peak, so the period occurs every 2 peaks. so
%starting from our current peak and stepping in twos... 
for i = peak_index:2:(length(x_peaks))
    x_gamma = [x_gamma x_peaks(i)] ;
    y_gamma = [y_gamma y_peaks(i)] ;
end

%%confirms we got the data we wanted
% hold on 
% plot(x_peaks, y_peaks)
% plot(x_gamma, y_gamma)
% hold off

x_gamma = x_gamma - x_gamma(1); %start from 0 might help the fit (it does)

% normalizing y_gamma (check this with brennan)
y_gamma = y_gamma./y_gamma(1);

%define fit type for the expontential
exp_fit_type = fittype("exp(-gamma*x)",...
    dependent='y', independent='x',...
    coefficients='gamma');

options = fitoptions(exp_fit_type);
options.StartPoint = 0.4;

exp_fit = fit(x_gamma', y_gamma', exp_fit_type, options);

%Pull/solve for gamma/w0 from fits
gamma = exp_fit.gamma;
w_0 = sqrt(b_est^2 + gamma^2);

%Solve for constants using function below and our estimate gamma and w_0

[c1, c2] = get_constant(gamma, w_0, x_peaks(peak_index), y_peaks(peak_index)); %Use solution at time t (user determines time t using index) where a peak occurs

%%Fitting with our Guesses

ODE_fit_type = fittype(@(gamma, w0, c1, c2, x) ...
        exp(-gamma.*x) .* ( ...
        c1.*cos( sqrt(w0.^2 - gamma.^2).*x ) + ...
        c2.*sin( sqrt(w0.^2 - gamma.^2).*x ) ), ...
        'independent', 'x', 'dependent', 'y', ...
        'coefficients', {'gamma','w0', 'c1', 'c2'});

options = fitoptions(ODE_fit_type);
options.StartPoint = [gamma, w_0, c1, c2];

ODE_fit = fit(t_seg', y_seg', ODE_fit_type, options); %Calling the actual fit, with user inputted guesses


gamma = ODE_fit.gamma; %Pulling out the actual fitted values
w_0    = ODE_fit.w0;
c1 = ODE_fit.c1;
c2 = ODE_fit.c2;

% Makes function using derived stuff
phi_an = @(x) (exp(-gamma.*x)).*(c1.*cos(sqrt(w_0^2 - gamma^2).*x) + c2.*sin(sqrt(w_0^2 - gamma^2).*x)); %define the solution in terms of c_1 and c_2

%Plotting for visual validation
f = figure(1);
theme(f,"light");
hold on
d1 = plot(full_X(1:2:(3*N/4 - 80)), full_y(1:2:(3*N/4 - 80)), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255], 'DisplayName','Original Sprinkler Data');
plot(full_X((3*N/4 - 79):6:end), full_y((3*N/4 - 79):6:end), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255])
d2 = plot(full_X, phi_an(full_X),'-', 'Color',[227/255,159/255,246/255], 'DisplayName', 'Analytical Fit of Sprinkler Data');
d3 = plot(x_peaks(peak_index), y_peaks(peak_index), 'o','Color', [92/255,188/255,99/255], 'DisplayName','$\phi(t_0)$ Used for Fit');
% plot(full_X(index), full_y(index), 'go')
ylim([-2,2])
xlim([20,40])
xlabel("t")
ylabel("$$\phi(t)$$")
legend([d1 d2 d3])
hold off

error = norm(phi_an(full_X(index:end)) - full_y(index:end))/norm(full_y(index:end));
fprintf("norm two error of data and analytical fit: %d \n", error)
%% Inverse Solve

full_y = full_y - full_y(1); %cleaning out some data

spot = full_X(index); %determine x value at peak so we can transplant our estimate of tail in for smoothness
dx = mean(diff(full_X));
x_end = (spot + dx):dx:(spot+(dx*(2048-length(full_X(1:index))))); %define the end of positive x to be spot + dx to the dx*(number of points we want - how many points we have)
x_neg = -(2048*dx):dx:(0 - dx);
y_neg = zeros(1,2048);
franken_x = [x_neg, full_X(1:index), x_end]; %put our new x domain together
franken_y = [y_neg, full_y(1:index), phi_an(x_end)]; %put our new y domain together


%%

N = 2*2048;     %Domain definitions so that everything is the same size when we go into fourier space
L = franken_x(end)*2;
h = L/N; 
x = h*(1:N)'-L/2;  %Defining x (or t) domain for function

%trying to remove some noise
for i = N/2:1:(2*N/3)
    if(franken_x(i) < 8)
        franken_y(i) = 0;
    end 
end 

%smooth data for fourier transform
franken_y_check = smoothdata(franken_y, "gaussian", 5);

error_smooth = norm(franken_y_check - franken_y)/norm(franken_y);
fprintf("norm 2 error from smoothing y: %d \n", error_smooth)
franken_y = franken_y_check;

G_hat = @(k) -1.0./(k.^2 + 2.0*1i*gamma.*k - w_0^2); %defining greens function

kk = ((2*pi)/L)*[0:N/2-1 0 -N/2+1:-1]; %defining domain for FFT


signal = ifft(fft(franken_y)./G_hat(-kk));  %Do the FFT and IFFT to output the torque


%picking an x point to plot at (the point of interest here is 12.1234)
for i = 1:1:N
    if(abs(franken_x(i) - 12.1234) < 5e-5)
        index_plot = i;
    end 
end


%Plots output
f = figure(2);
theme(f,"light");
hold on
plot(franken_x, signal, 'Color',[227/255,159/255,246/255])
legend('Inverse Solve Solution')
title("Solving for Torque Function with Inverse FFT")
xlabel("t")
ylabel("$$\tau(t)$$")
yticks(-100:20:100)
xlim([0,40])
hold off

f = figure(3);
theme(f,"light");
hold on
plot(franken_x, franken_y,'Color',[227/255,159/255,246/255])
xlim([0,40])
xlabel("t")
ylabel("$$\phi(t)$$")
legend('Combined Data')
title("Plot of Combined Sprinkler Data Used in FFT")
hold off

%% Going backwards
signal_chunk = signal(N/2:(N/2 + length(full_X) - 1));

y0 = [0; 0]; %IC for y and y'

para = @(t) interp1(full_X, signal_chunk, t, 'pchip', 0);

dydt = @(t,y)[y(2); -(2*gamma)*y(2) - (w_0)^2*y(1) + para(t)]; %Actual ODE is here.


[t,y] = ode45(dydt, full_X, y0); %ODE45 solver call, generating phi data

phi_gen = y(:,1)'; %generated phi data

N = length(full_X);

f = figure(4);
theme(f,"light");
hold on
%plot((full_X), phi_gen,'-', 'Color',[227/255,159/255,246/255][55/255,55/255,55/255], 'DisplayName','Sprinkler Data Using Computed Torque Signal')
d1 = plot(full_X(1:20:(N/4)), full_y(1:20:(N/4)), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255], 'DisplayName','Original Sprinkler Data');
plot(full_X((N/4 + 1):2:(N/3 + 10)), full_y((N/4 + 1):2:(N/3 + 10)), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255])
plot(full_X((N/3 + 11):5:(N/2 + 35)), full_y((N/3 + 11):5:(N/2 + 35)), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255])
plot(full_X((N/2 + 36):2:(3*N/4)), full_y((N/2 + 36):2:(3*N/4)), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255])
plot(full_X((3*N/4 + 1):6:end), full_y((3*N/4 + 1):6:end), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255])
d2 = plot((full_X), phi_gen,'-', 'Color',[227/255,159/255,246/255], 'DisplayName','Sprinkler Data Using Computed Torque Signal');
legend([d1 d2])
title("ODE45 Soltuion Using Comptued Torque")
xlabel("t")
ylabel("$$\phi(t)$$")
hold off


error1 = norm(phi_gen(1:index) - full_y(1:index))/norm(full_y(1:index));
fprintf("norm two error of original data and output using solved torque: %d \n", error1)


%General Solution Constant Solves for ODE Analytical

function [c_1, c_2] = get_constant(gamma, w_0, t, phi_t)

gam = sqrt(w_0^2 - gamma^2);

e_up = (exp(gamma*t));
cos_up = cos(gam*t);
sin_up = sin(gam*t);

a = -gamma*e_up*cos_up - gam*e_up*sin_up; 
b = -gamma*e_up*sin_up + gam*e_up*cos_up; 

c_2 = (phi_t*e_up)/(sin_up - (b/a)*cos_up);
c_1 = -(c_2*b)/a;

end




