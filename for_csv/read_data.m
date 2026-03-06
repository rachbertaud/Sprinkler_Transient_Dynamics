function [full_x, full_y, x_peaks, y_peaks, N, index, peak_index, t_seg, y_seg] = read_data(fname, t_target)
    data = readmatrix(fname);

    %pull out all data
    full_x = data(:,1)';
    full_y = data(:,2)';
    % plot(full_x, full_y)

    N = length(full_x); %Pull out full length of data!

    [~, loc] = findpeaks(abs(full_y));

    
    %pull out peaks data
    peak_N = length(loc);
    x_peaks = zeros(1, peak_N);
    y_peaks = zeros(1, peak_N);

    for j = 1:1:peak_N
        idx = loc(j);
        x_peaks(j) = full_x(idx);
        y_peaks(j) = full_y(idx);
    end

    
    [~, idx] = min(abs(x_peaks - t_target));
    peak_index = idx;
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

