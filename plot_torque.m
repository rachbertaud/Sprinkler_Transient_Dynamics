function plot_torque(franken_x, signal)
    %Plots output
    f = figure(2);
    theme(f,"light");
    hold on
    plot(franken_x, signal, 'Color',[227/255,159/255,246/255])
    legend('Inverse Solve Solution')
    title("Solving for Torque Function with Inverse FFT")
    xlabel("t")
    ylabel("$$\tau(t)$$")
    yticks((round(min(signal), -1) - 50):20:(round(max(signal), -1) + 50))
    xlim([0,40])
    hold off
end 