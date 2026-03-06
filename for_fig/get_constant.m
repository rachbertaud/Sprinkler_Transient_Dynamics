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