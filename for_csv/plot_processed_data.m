function plot_processed_data(franken_x, franken_y)
    f = figure(3);
    theme(f,"light");
    hold on
    plot(franken_x, franken_y,'Color',[227/255,159/255,246/255])
    xlim([0,40])
    xlabel("t")
    ylabel("$$\phi(t)$$")
    legend('Combined Data')
    title("Plot of Combined Sprinkler Data Used in FFT")
    hold off
end 