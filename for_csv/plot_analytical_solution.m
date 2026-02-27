function plot_analytical_solution(full_x, full_y, peak_index, x_peaks, y_peaks, phi_an, spin_switch, trial, Re)
    
    figure(1)
    subplot(2,2,trial)
    hold on
    plot(full_x, full_y, 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255], 'DisplayName','Original Sprinkler Data')
    plot(x_peaks(peak_index), y_peaks(peak_index),'o','Color', [92/255,188/255,99/255], 'DisplayName','$\phi(t_0)$ Used for Fit')
    plot(full_x, phi_an(full_x), '-', 'Color',[227/255,159/255,246/255], 'DisplayName', 'Analytical Fit of Sprinkler Data')
    ylim([-2,2])
    xlabel("t")
    ylabel("$$\phi(t)$$")
    if spin_switch == 1
        title(sprintf("Forward, Re = %d, Trial %d", Re, trial))
    else
        title(sprintf("Reverse, Re = %d,  Trial %d", Re, trial))
    end
    legend()
    hold off
end



% function plot_analytical_solution(full_x, full_y, peak_index, x_peaks, y_peaks, phi_an, spin_switch)
% 
%     figure(1)
%     hold on
%     plot(full_x, full_y, '-o')
%     plot(x_peaks(peak_index), y_peaks(peak_index), 'om', 'MarkerSize',40)
%     plot(full_x, phi_an(full_x))
%     ylim([-2,2])
%     xlabel("t")
%     ylabel("$$\phi(t)$$")
%     legend("OG Data", "Peak Used", "Analytical Solution")
%     hold off
% end