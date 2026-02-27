function gamma = estimate_gamma(peak_index, x_peaks, y_peaks)
    x_gamma = []; %to fill with gamma data
    y_gamma = []; %to fill with gamma data
    %we know we are starting at a peak, so the period occurs every 2 peaks. so
    %starting from our current peak and stepping in twos... 
    for i = peak_index:2:(length(x_peaks))
        x_gamma = [x_gamma x_peaks(i)] ;
        y_gamma = [y_gamma y_peaks(i)] ;
    end
    
    
    x_gamma = x_gamma - x_gamma(1); %start from 0 might help the fit (it does)
    
    % normalizing y_gamma (check this with brennan)
    y_gamma = y_gamma./y_gamma(1);


    
    
    %define fit type for the expontential
    exp_fit_type = fittype("exp(-gamma*x)",...
        dependent='y', independent='x',...
        coefficients='gamma');
    options = fitoptions(exp_fit_type);
    options.StartPoint = 0.5;
    exp_fit = fit(x_gamma', y_gamma', exp_fit_type, options);
    
    %Pull/solve data from fits
    gamma = exp_fit.gamma;
end 