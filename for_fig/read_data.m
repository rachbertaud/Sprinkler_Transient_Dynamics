function [full_x, full_y, x_peaks, y_peaks, N, index, t_seg, y_seg] = read_data(spin_switch, peak_index)
    if(spin_switch == 0)
        fig = openfig('reverse_800.fig'); %open figure
    elseif(spin_switch == 1)
        fig = openfig('forward_800.fig'); %open figure
    else
        fprintf("Error: Invalid input for pin_switch. Choose 0 for reverse spin or 1 for forward spin")
    end 

    dataObjs = findobj(fig, 'Type', 'line'); %find the line
    
    xData = get(dataObjs, 'XData'); %get xData from the line
    yData = get(dataObjs, 'YData'); %get yData from the line
    
    close all 
    
    %pull out peaks data
    x_peaks = xData{1}; 
    y_peaks = yData{1};
    
    %pull out all data
    full_x = xData{2};
    full_y = yData{2};

    N = length(full_x); %Pull out full length of data!


    %determine which index of the full data set our peak index for inverse solve occurs at
    for i = 1:1:N
        if(full_x(i) == x_peaks(peak_index))
            index = i;
        end 
    end

    % Pulls out data after our cutting index for inverse solve
    t_seg = full_x(index:end);        
    y_seg = full_y(index:end); 

end 

