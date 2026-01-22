function [franken_x, franken_y] = process_data(franken_x, franken_y, spin_switch)
    
    N = 2*2048; 

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
    fprintf("norm 2 error from smoothing y: %.5f \n", error_smooth)
    franken_y = franken_y_check;
end