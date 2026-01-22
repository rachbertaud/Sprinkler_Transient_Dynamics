function phi_gen = phi_from_torque(N_f, full_x, signal, gamma, w_0)
    signal_chunk = signal(N_f/2:(N_f/2 + length(full_x) - 1)); %signal from 0 to end of full_x
    
    y0 = [0; 0]; %IC for y and y'
    
    para = @(t) interp1(full_x, signal_chunk, t, 'pchip', 0);
    
    dydt = @(t,y)[y(2); -(2*gamma)*y(2) - (w_0)^2*y(1) + para(t)]; %Actual ODE is here.
    
    
    [t,y] = ode45(dydt, full_x, y0); %ODE45 solver call, generating phi data
    
    phi_gen = y(:,1)'; %generated phi data
end 
