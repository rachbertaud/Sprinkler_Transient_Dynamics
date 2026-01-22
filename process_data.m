function [franken_x, franken_y] = process_data(index, full_x, full_y, phi_an, spin_switch)
    
    N = 2*2048; 
    
    full_y = full_y - full_y(1); %cleaning out some data
    
    spot = full_x(index); %determine x value at peak so we can transplant our estimate of tail in for smoothness
    dx = mean(diff(full_x));
    x_end = (spot + dx):dx:(spot+(dx*(2048-length(full_x(1:index))))); %define the end of positive x to be spot + dx to the dx*(number of points we want - how many points we have)
    x_neg = -(2048*dx):dx:(0 - dx);
    y_neg = zeros(1,2048);
    franken_x = [x_neg, full_x(1:index), x_end]; %put our new x domain together
    franken_y = [y_neg, full_y(1:index), phi_an(x_end)]; %put our new y domain together

    if(spin_switch == 1) %if forward
          %trying to remove some noise
        for i = N/2:1:(2*N/3)
            if(franken_x(i) < 12)
                franken_y(i) = 0;
            end 
        end 

        %smooth data for fourier transform
        franken_y_check = smoothdata(franken_y, "gaussian", 10);
    else %if reverse
              %trying to remove some noise
        for i = N/2:1:(2*N/3)
            if(franken_x(i) < 8)
                franken_y(i) = 0;
            end 
        end 

        %smooth data for fourier transform
        franken_y_check = smoothdata(franken_y, "gaussian", 5);
    end 
    
    
    
    error_smooth = norm(franken_y_check - franken_y)/norm(franken_y);
    fprintf("norm 2 error from smoothing y: %d \n", error_smooth)
    franken_y = franken_y_check;
end