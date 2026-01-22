function plot_phi_from_torque(full_x, full_y, phi_gen, spin_switch)
%you might be wondering wtf all these random N indexing is for, and its
%simply for visual readability of the plots. I could probably make this
%nicer, and i will, but not today. Sorry folks. spin_switch is a variable
%for this plot because the picking and choosing of plotting indicies was
%slightly different for the different spins. Doth is lyfe. 
    N = length(full_x);
    if(spin_switch == 1)
        f = figure(4);
        theme(f,"light");
        hold on
        %plot((full_x), phi_gen,'-', 'Color',[227/255,159/255,246/255][55/255,55/255,55/255], 'DisplayName','Sprinkler Data Using Computed Torque Signal')
        d1 = plot(full_x(1:20:(N/3)), full_y(1:20:(N/3)), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255], 'DisplayName','Original Sprinkler Data');
        plot(full_x((N/3 + 1):2:(N/3 + 60)), full_y((N/3 + 1):2:(N/3 + 60)), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255])
        plot(full_x((N/3 + 61):5:(N/2 + 110)), full_y((N/3 + 61):5:(N/2 + 110)), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255])
        plot(full_x((N/2 + 111):2:round(round(3*N/4))), full_y((N/2 + 111):2:(round(round(3*N/4)))), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255])
        plot(full_x((round(round(3*N/4)) + 1):6:N), full_y((round(round(3*N/4)) + 1):6:N), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255])
        d2 = plot((full_x), phi_gen,'-', 'Color',[227/255,159/255,246/255], 'DisplayName','Sprinkler Data Using Computed Torque Signal');
        legend([d1 d2])
        title("ODE45 Soltuion Using Comptued Torque")
        xlabel("t")
        ylabel("$$\phi(t)$$")
        hold off
    else
        f = figure(4);
        theme(f,"light");
        hold on
        %plot((full_x), phi_gen,'-', 'Color',[227/255,159/255,246/255][55/255,55/255,55/255], 'DisplayName','Sprinkler Data Using Computed Torque Signal')
        d1 = plot(full_x(1:20:round(N/4)), full_y(1:20:round(N/4)), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255], 'DisplayName','Original Sprinkler Data');
        plot(full_x((round(N/4) + 1):2:(round(N/3) + 10)), full_y((round(N/4) + 1):2:(round(N/3) + 10)), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255])
        plot(full_x((round(N/3) + 11):5:(round(N/2) + 35)), full_y((round(N/3) + 11):5:(round(N/2) + 35)), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255])
        plot(full_x((round(N/2) + 36):2:(round(3*N/4))), full_y((round(N/2) + 36):2:(round(3*N/4))), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255])
        plot(full_x((round(3*N/4) + 1):6:end), full_y((round(3*N/4) + 1):6:end), 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255])
        d2 = plot((full_x), phi_gen,'-', 'Color',[227/255,159/255,246/255], 'DisplayName','Sprinkler Data Using Computed Torque Signal');
        legend([d1 d2])
        title("ODE45 Soltuion Using Comptued Torque")
        xlabel("t")
        ylabel("$$\phi(t)$$")
        hold off
    end 
end 
