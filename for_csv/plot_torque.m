function plot_torque(franken_x, signal, N_f)
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

    cumtrapz_output = cumtrapz(franken_x(N_f/2:end), signal(N_f/2:end));

    max_trap = max(cumtrapz_output);
    index_of_max = find(cumtrapz_output==max_trap);

    f = figure(5);
    theme(f,"light");
    hold on
    p = plot(franken_x(N_f/2:end), cumtrapz_output,'Color',[227/255,159/255,246/255] );
    datatip(p,franken_x(N_f/2 + index_of_max), max_trap, 'Location', 'northwest'); 
    xlabel("$$t$$")
    title('Cumulative Integral of Torque Signal')
    ylabel('$$\int_0^t \tau(t)$$')
    xlim([0,70])
    ylim([min(cumtrapz_output - 10), max(cumtrapz_output + 200)])
    hold off
end 