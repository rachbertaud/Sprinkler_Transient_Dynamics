fs = 20;  % set your preferred font size here

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

clc
clear all

y_0 = 0; %IC for y
y_prime_0 = 0.00001; %IC for y'

gamma = 0.1; %Gamma, used to create phi data
w_0 = 10; %Omega_0, used to create phi data

gamma_guess = 0.11; %Gamma, used to create phi data
w_0_guess = 10.1; %Omega_0, used to create phi data

index = 26; %At what t index for peaks we want the inverse solve


N = 2*2048;     %Domain definitions so that everything is the same size when we go into fourier space
L = 200;
h = L/N; 
x = h*(1:N)'-L/2;  %Defining x (or t) domain for function


%% Generating phi data
[t,y] = ode45(@sprk, x, [y_0;y_prime_0]); %ODE45 solver call, generating phi data

phi_gen = y(:,1); %generated phi data

%% Finding Peaks for Inverse Parameter Solving
% this section goes through the generated phi data and determines where the
% peaks exist. we want these for the inverse solve, since we know what the
% first derivative of phi at these points is (y'(peak) = 0). 
x_peaks = [];
i_peaks = [];
for i = (N/2 + 1:1:N - 1)
    if ((abs(phi_gen(i+1) - phi_gen(i-1))/(2*h)) < 8e-2)
        x_peaks = [x_peaks; x(i)];
        i_peaks = [i_peaks; i];
    end
end


 

%% Plotting Phi for Validation after Inverse Solve
% this section takes a user defined peak index, solves the ode numerically
% for c_1 and c_2 and sees how well the numerical estimate for phi is. 
close all
clf

[c1, c2] = get_constant(gamma, w_0, x_peaks(index), phi_gen(i_peaks(index))); %Use solution at time t (user determines time t using index) where a peak occurs

gam = sqrt(w_0^2 - gamma^2);
phi_an = @(x) (exp(-gamma.*x)).*(c1.*cos(sqrt(w_0^2 - gamma^2).*x) + c2.*sin(sqrt(w_0^2 - gamma^2).*x)); %define the solution in terms of c_1 and c_2

% Indexing for plotting phi_an after the t it was determined at
i_index = i_peaks(index);
phi_an_seg = phi_an(t(i_index:end));
t_seg = t(i_index:end); 
phi_gen_seg = phi_gen(i_index:end);

% Indexing for all the time before phi_an was solved for
t_beg = t(N/2:i_index);
phi_an_beg = phi_an(t_beg);
phi_gen_beg = phi_gen(N/2:i_index);

%Plotting for visual validation
figure(1)
hold on
plot(t_seg,phi_an_seg', 'bx-')
plot(t_seg, phi_gen_seg,'ro-')
plot(t_beg, phi_an_beg, 'c-x')
plot(t_beg, phi_gen_beg, 'm-o')
legend('An. sol after t_0', 'ode45 sol after t_0', 'An. sol b4 t_0', 'ode45sol b4 t_0')
hold off

%Calculating the error between the analytical and generated outputs.
%In other words, how much error is there in the estimation of phi_gen using a
%solve for phi_an at a peak
error = norm(phi_an_seg - phi_gen_seg)/norm(phi_gen_seg);
fprintf("Error in Phi Estimation: %f \n", error)

%% Fitting for gamma and omega_0
% This section takes our generated solution, and fits the data to solve for
% omega and gamma.
t_seg = t(N/2 + N/3:end);   %going moderately far into the system     
y_seg = phi_gen(N/2 + N/3:end); 

ft = fittype(@(gamma, w0, x) ...
    exp(-gamma.*x) .* ( ...
        c1.*cos( sqrt(w0.^2 - gamma.^2).*x ) + ...
        c2.*sin( sqrt(w0.^2 - gamma.^2).*x ) ), ...
        'independent', 'x', ...
        'coefficients', {'gamma','w0'});
% Above defines the fit we want for gamma and w0 using c_1 and c_2 we
% derived above


[fit_1, gof] = fit(t_seg, y_seg, ft, 'StartPoint', [gamma_guess, w_0_guess]); %Calling the actual fit, with user inputted guesses


gamma_hat = fit_1.gamma; %Pulling out the actual fitted values
w0_hat    = fit_1.w0;

%Constructing G_hat using the fitted values
G_hat = @(k) -1.0./(k.^2 + 2.0*1i*gamma_hat.*k - w0_hat^2);

% Different forcing funtion options
gauss = @(t,t_0,sigma) (1.0/(sigma*sqrt(2.0*pi)))*exp(-0.5*((t-t_0)/sigma).^2);

tau = @(t) (gauss(t,2.5,0.1) + gauss(t,1,0.25)); %tau = @(t) gauss(t,2.5,1.0) + gauss(t,8.5,0.5);

para = @(t) (-(60).*(t-2).^2 + 60).*(t>=1 & t<= 3);  % para = @(t) (-(30).*(t-2).^2 + 30).*(t>=1 & t<= 3); 


kk = ((2*pi)/L)*[0:N/2-1 0 -N/2+1:-1]'; 
%This defines domain in Fourier space for greens functions


%Takes our data, and solves for tau using our solution for G_hat
%constructed with estimates for gamma dn omega_0
signal = ifft(fft(phi_gen)./G_hat(-kk));  

%Plots output
figure(2)
hold on
plot(t, para(t), 'x')
plot(t, signal)
legend('Actual Tau','Inverse Solve')
xlim([0,5])
hold off

error_inv = norm(signal - para(t))/norm(para(t));
fprintf("Error in Tau Estimation: %f \n", error_inv)


%% Functions 

% to solve for ode45

function dydt = sprk(t,y)
%forcing torque functions, RHS 
gamma = 0.1; %ODE Parameters
w_0 = 10; %ODE Parameters
%Gaussian 
gausDist = @(x,mu,sig)(1./(sig.*sqrt(2*pi))).*exp((-(x-mu).^2)./(2.*sig.^2));

tau_sing = gausDist(t, 2.5,1);

tau = gausDist(t, 2.5,0.1);

tau1 = gausDist(t, 1,0.25);

tau_mix = (tau + tau1);

%Right Tailed Heavy
pdr = makedist('Stable','alpha',0.5,'beta',0,'gam',1,'delta',3); %delta is where it is centered
pdfr = pdf(pdr,t);

%Left Tailed heavy
pdl = makedist('Stable','alpha',0.3,'beta',-1,'gam',1,'delta',3); %delta is where it is centered
pdfl = pdf(pdl,t);

%parabola 
para_nb = @(t) -(60).*(t-2).^2 + 60;

para = @(t) (-(60).*(t-2).^2 + 60).*(t>=1 & t<= 3); 
% ode45 function to solve
dydt = [y(2); -(2*gamma)*y(2) - (w_0)^2*y(1) + para(t)]; %Actual ODE is here.
end



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
