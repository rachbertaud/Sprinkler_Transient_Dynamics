function plot_analytical_solution(full_x, full_y, peak_index, x_peaks, y_peaks, phi_an, spin_switch)
    N = length(full_x);
    if(spin_switch == 1)%Plotting for visual validation
        f = figure(1);
        theme(f,"light");
        hold on
        d1 = plot(full_x(1:2:(round(3*N/4) - 80)), full_y(1:2:(round(round(3*N/4)) - 80)), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255], 'DisplayName','Original Sprinkler Data');
        plot(full_x((round(3*N/4) - 79):6:end), full_y((round(round(3*N/4)) - 79):6:end), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255])
        d2 = plot(full_x, phi_an(full_x),'-', 'Color',[227/255,159/255,246/255], 'DisplayName', 'Analytical Fit of Sprinkler Data');
        d3 = plot(x_peaks(peak_index), y_peaks(peak_index), 'o','Color', [92/255,188/255,99/255], 'DisplayName','$\phi(t_0)$ Used for Fit');
        % plot(full_X(index), full_y(index), 'go')
        ylim([-2,2])
        xlim([20,40])
        xlabel("t")
        ylabel("$$\phi(t)$$")
        legend([d1 d2 d3])
        hold off
    else 
        %Plotting for visual validation
        f = figure(1);
        theme(f,"light");
        hold on
        d1 = plot(full_x(1:2:(round(3*N/4) - 80)), full_y(1:2:(round(3*N/4) - 80)), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255], 'DisplayName','Original Sprinkler Data');
        plot(full_x((round(3*N/4) - 79):6:end), full_y((round(3*N/4) - 79):6:end), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255])
        d2 = plot(full_x, phi_an(full_x),'-', 'Color',[227/255,159/255,246/255], 'DisplayName', 'Analytical Fit of Sprinkler Data');
        d3 = plot(x_peaks(peak_index), y_peaks(peak_index), 'o','Color', [92/255,188/255,99/255], 'DisplayName','$\phi(t_0)$ Used for Fit');
        % plot(full_X(index), full_y(index), 'go')
        ylim([-2,2])
        xlim([20,40])
        xlabel("t")
        ylabel("$$\phi(t)$$")
        legend([d1 d2 d3])
        hold off
    end 
end