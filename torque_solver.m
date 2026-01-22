function signal = torque_solver(N_f, gamma, w_0, L, franken_y)

    G_hat = @(k) -1.0./(k.^2 + 2.0*1i*gamma.*k - w_0^2);

    kk = ((2*pi)/L)*[0:N_f/2-1 0 -N_f/2+1:-1]; 


    signal = ifft(fft(franken_y)./G_hat(-kk));  
end 