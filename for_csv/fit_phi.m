function [gamma, w_0, c1, c2] = fit_phi(x_seg, y_seg, gamma, w_0, c1, c2)
    %%Fitting with our Guesses
    
    ODE_fit_type = fittype(@(gamma, w0, c1, c2, x) ...
            exp(-gamma.*x) .* ( ...
            c1.*cos( sqrt(w0.^2 - gamma.^2).*x ) + ...
            c2.*sin( sqrt(w0.^2 - gamma.^2).*x ) ), ...
            'independent', 'x', 'dependent', 'y', ...
            'coefficients', {'gamma','w0', 'c1', 'c2'});
    
    options = fitoptions(ODE_fit_type);
    options.StartPoint = [gamma, w_0, c1, c2];
    
    ODE_fit = fit(x_seg', y_seg', ODE_fit_type, options); %Calling the actual fit, with user inputted guesses
    
    
    gamma = ODE_fit.gamma; %Pulling out the actual fitted values
    w_0    = ODE_fit.w0;
    c1 = ODE_fit.c1;
    c2 = ODE_fit.c2;
end 
