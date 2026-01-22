function [franken_x, franken_y] = combine_data(full_x, full_y, index, phi_an)
    full_y = full_y - full_y(1); %cleaning out some data
    
    spot = full_x(index); %determine x value at peak so we can transplant our estimate of tail in for smoothness
    dx = mean(diff(full_x));
    x_end = (spot + dx):dx:(spot+(dx*(2048-length(full_x(1:index))))); %define the end of positive x to be spot + dx to the dx*(number of points we want - how many points we have)
    x_neg = -(2048*dx):dx:(0 - dx);
    y_neg = zeros(1,2048);
    franken_x = [x_neg, full_x(1:index), x_end]; %put our new x domain together
    franken_y = [y_neg, full_y(1:index), phi_an(x_end)]; %put our new y domain together
end 