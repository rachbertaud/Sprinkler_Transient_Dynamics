function Omega_est = estimate_Omega(peak_index, x_peaks, y_peaks)

    aoi = x_peaks(peak_index:end); % this is the peak area of interest, since this is all the points at and after our determined cutting point
    end_peaks = length(aoi); %determine the length of the area of interest
    period_est = zeros(1, end_peaks - 2); %define zzero vector to fill with our period estimate
    j = 1; %indexing variable

    
    for i = peak_index:1:(length(x_peaks) - 2) %for all the points after our area of interest
        period_est(j) = x_peaks(i + 2) - x_peaks(i); %find out the period after our picked solve index
        j = j + 1;
        % hold on
        % plot(x_peaks, y_peaks, 'o-')
        % plot(x_peaks(i + 2), y_peaks(i+2), 'og')
        % plot(x_peaks(i), y_peaks(i), 'om')
        % legend("full","i+2", "i")
        % xlim([x_peaks(peak_index), 80])
        % ylim([-0.4, 0.4])
        % pause(1)
        % clf
    end
    
    period_est = rmoutliers(period_est ,"quartiles");
    Omega_est = (2*pi)/mean(period_est); %estimate the period using the mean of the period data found above
end