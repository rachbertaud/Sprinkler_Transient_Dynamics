%% Second Meeting
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
%%

y_0 = 0; %IC for y
y_prime_0 = 0.00001; %IC for y'

N = 2*2048;     %Domain definitions so that everything is the same size when we go into fourier space
L = 200;
h = L/N;   %Downsample the other way pls
x = h*(1:N)'-L/2;  %Defining x (or t) domain for function

gamma = 0.1;
w_0 = 10;


[t,y] = ode45(@sprk, x, [y_0;y_prime_0]); %ODE45 solver call

[t1,y1] = ode45(@sprk1, x, [y_0;y_prime_0]); %ODE45 solver call for different forcing function
ode_sol = y(:,1);


%%
%Plotting solution to ode45 where y(:,1) is y and y(2,:) is y_prime
figure(1)
hold on
plot(t,ode_sol,'-om', 'MarkerEdgeColor','m')%,t,y(:,2),'-o')
% plot(t1,y1(:,1),'-ob', 'MarkerEdgeColor','b')%,t,y(:,2),'-o')
xlim([0, 20]);
ylim([-1,1]);
title('Solution of Transient Sprinkler Dynamics with ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('y_1')%, 'y1_1')
hold off


%% From brennan

G_hat = @(k) -1.0./(k.^2 + 2.0*1i*gamma.*k - w_0^2);

%Forcing functions
gauss = @(t,t_0,sigma) (1.0/(sigma*sqrt(2.0*pi)))*exp(-0.5*((t-t_0)/sigma).^2);

%tau = @(t) gauss(t,2.5,1.0) + gauss(t,8.5,0.5);

tau = @(t) (gauss(t,2.5,0.1) + gauss(t,1,0.25));

% para = @(t) (-(30).*(t-2).^2 + 30).*(t>=1 & t<= 3); 

para = @(t) (-(60).*(t-2).^2 + 60).*(t>=1 & t<= 3); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk = ((2*pi)/L)*[0:N/2-1 0 -N/2+1:-1]'; 
%This defines domain in Fourier space for greens functions???

%Verified phi produced by ode45 (forward)
% signal = ifft(G_hat(-kk).*fft(tau(x)));  

%Takes phi produced by ode45 to solve for tau (inverse)
signal = ifft(fft(y(:,1))./G_hat(-kk));   

%% Plots everything

clf
plot(x,signal,'xm') %Ploting output of Foruier Transform
hold on
plot(x, y(:,1),"*") %Plots output of ode45
plot(x,para(t),'-') %Plots tau (forcing function)
xlim([0, 10]) %defining limits for x
legend('Inverse Solve for \tau(t)','ode45 solve for \phi(t)', 'Real \tau(t)')
xlabel('time (s)')
ylabel('Spring deflection')
hold off


%% Downsampling


%% Testing Noise

kk = ((2*pi)/L)*[0:N/2-1 0 -N/2+1:-1]';
G_hat = @(k) -1.0./(k.^2 + 2.0*1i*gamma.*k - w_0^2);
para = @(t) (-(60).*(t-2).^2 + 60).*(t>=1 & t<= 3); 

i_int = 1:1:10; %3*10^-1, loop thru explicit noise 50%, 20%, 10%, 1%, .1%
legend_input = {};
input = zeros(length(t),length(i_int));
signal = zeros(length(t),length(i_int));
for i = i_int
    noise = (10^(-i)) * rand(size(t));
    input(:,i) = (y(:,1) + noise)';
    signal(:,i) = (ifft(fft(input(:,i))./(G_hat(-kk))))'; 
end 

%%
legend_input = {};
plot_on = 0;
clf
hold on

counter = 1;

i_sol = 1;

for i = i_int
    if (norm(signal(:,i) - para(t))/norm(para(t)))<10^(-i_sol)
        if plot_on == 0
            figure(1)
            hold on
            plot(x, para(t), 'o')
            legend_input{end+1} = sprintf('Solution \\tau(t)');
            plot(x,signal(:,i)) %Ploting output of Foruier Transform
            legend_input{end+1} = sprintf('Result with \\phi(t) + Noise of $10^{-%d}$', i);
            plot_on = 1;
            % 
            % plot(x,para(t),'-',LineWidth=2) %Plots tau (forcing function)
            xlim([0 6]) %defining limits for x
            ylim([-1 65])
            
            xlabel('time (s)')
            ylabel('Spring deflection')
            title(sprintf('Solution to Inverse Problem Showing Noise Upper Bound for Error of $10^{-%d}$ Between Result and Solution', i_sol))
            L = legend(legend_input{:});
            L.FontSize = 25;
            hold off
            counter = counter + 1;
        end 
    end
end




%% Functions to solve for ode45

% function y = para(x) 
% if x < 1 || x > 3
%     y = 0;
% else 
%    y = -(30).*(x-2).^2 + 30;
% end 
% end 

function dydt = sprk(t,y)
% forcing torque functions, RHS 
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


%%
%ode45 function to solve
dydt = [y(2); -(2*gamma)*y(2) - (w_0)^2*y(1) + para(t)]; %Actual ODE is here.
end


function dydt = sprk1(t,y)
% forcing torque functions, RHS 
gamma = 0.1; %ODE Parameters
w_0 = 5; %ODE Parameters
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

para = @(t) (-(30).*(t-2).^2 + 30).*(t>=1 & t<= 3); 

%ode45 function to solve
dydt = [y(2); -(2*gamma)*y(2) - (w_0)^2*y(1) + para(t)]; %Actual ODE is here.
end

%% General Solution to ODE Analytical

function [c_1, c_2] = get_constant(gamma, w_o, t, phi_t)

gam = sqrt(1 - gamma^2);

e_up = (exp(-gamma*w_o*t));
cos_up = cos(w_o*gam*t);
sin_up = sin(w_o*gam*t);

tan_up = sin_up/cos_up;


a = -gamma*w_o*e_up*cos_up - w_o*gam*e_up*sin_up;
b = cos_up;
c = -gamma*w_o*e_up*sin_up + w_o*gam*e_up*cos_up;
d = phi_t/e_up;
e = sin_up;


c_2 = d/(e - ((c*b)/a));

c_1 = (c_2*w_o*gam - tan_up*c_2*gamma*w_o)/(gamma*w_o + tan_up*w_o*gam);


end

%% Testing Distributions 

t_new = 0:0.1:5;

%Gaus Dist.
gausDist = @(x,mu,sig)(1./(sig.*sqrt(2*pi))).*exp((-(x-mu).^2)./(2.*sig.^2));


tau = gausDist(t_new, 2.5,0.1);

tau1 = gausDist(t_new, 1,0.25);

tau_mix = tau + tau1;



%Right tailed Dist. 
pd3 = makedist('Stable','alpha',0.3,'beta',-1,'gam',1,'delta',3); %delta is where it is centered
pdf3 = pdf(pd3,t_new);


%% First Meeting
% 
% alpha = 1/8; 
% w0 = 1;
% 
% wr = w0/sqrt(1 - alpha^2);
% 
% e = exp(1);
% 
% % phi_add1 = @(t) sin(2*t);
% % phi_add2 = @(t) -sin(2*t);
% % phi = @(t) (e.^(-alpha*w0*t)).*sin(wr*t);
% % 
% p = pi;
% phi_1 = @(t) 0.2*sin(t);
% phi_2 = @(t) -sin(t-pi)*(e.^(-alpha*w0*t));
% t_end = 1000;
% t_step = 0.01;
% % t = 0:.01:50; 
% hold on
% figure(1)
% p = 2500*pi/5103;
% figure(1)
% y = [];
% for t = 0:.01:t_end
%     if t <= 2*p 
%         y = [y ,  0.1 * phi(t)]; 
%     elseif t <= 4*p
%         y = [y, 2 * phi(t)];
%     else
%         y = [y, phi(t)];
%     end
% end
% t = 0:t_step:t_end;
% % plot(t,y)
% phi_f = fft(y);
% 
% gausDist = @(x,mu,sig)(1./(sig.*sqrt(2*pi))).*exp((-(x-mu).^2)./(2.*sig.^2));
% tau = gausDist(t, 1,0.25);
% plot(t,tau)
% tau_f = fft(tau);
% 
% G_f = phi_f./tau_f;
% 
% G = ifft(G_f);
% 
% plot(G)
% 
% hold off
% 
% 
