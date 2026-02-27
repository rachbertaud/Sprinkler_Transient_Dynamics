function [c_1, c_2] = get_constant(gamma, w_0, t_peak, phi_t)
    gam = sqrt(w_0^2 - gamma^2);
    
    % In shifted coordinates tau = t - t_peak:
    % phi(tau) = e^(-gamma*tau) * (c1*cos(gam*tau) + c2*sin(gam*tau))
    % At tau=0: phi(0) = c1 = phi_t
    % At tau=0: phi'(0) = -gamma*c1 + gam*c2 = 0 => c2 = (gamma/gam)*c1
    c_1 = phi_t;
    c_2 = (gamma/gam) * phi_t;
end
% function [c_1, c_2] = get_constant(gamma, w_0, t, phi_t)
% 
%     gam = sqrt(w_0^2 - gamma^2);
% 
%     e_up = (exp(gamma*t));
%     cos_up = cos(gam*t);
%     sin_up = sin(gam*t);
% 
%     a = -gamma*e_up*cos_up - gam*e_up*sin_up; 
%     b = -gamma*e_up*sin_up + gam*e_up*cos_up; 
%     fprintf("e_up: %.4f, cos_up: %.4f, sin_up: %.4f\n", e_up, cos_up, sin_up)
%     fprintf("a: %.4f, b: %.4f\n", a, b)
%     fprintf("denominator: %.6f\n", sin_up - (b/a)*cos_up)
% 
%     c_2 = (phi_t*e_up)/(sin_up - (b/a)*cos_up);
%     c_1 = -(c_2*b)/a;
% 
% end