function plot_phi_from_torque(full_x, full_y, phi_gen, spin_switch)
%you might be wondering wtf all these random N indexing is for, and its
%simply for visual readability of the plots. I could probably make this
%nicer, and i will, but not today. Sorry folks. spin_switch is a variable
%for this plot because the picking and choosing of plotting indicies was
%slightly different for the different spins. Doth is lyfe. 

    f = figure(4);
    theme(f,"light");
    hold on
    %plot((full_x), phi_gen,'-', 'Color',[227/255,159/255,246/255][55/255,55/255,55/255], 'DisplayName','Sprinkler Data Using Computed Torque Signal')
    plot(full_x, full_y, 'square', 'MarkerSize', 5, 'Color',[252/255,70/255,170/255], 'DisplayName','Original Sprinkler Data')
    plot(full_x, phi_gen,'-', 'Color',[227/255,159/255,246/255], 'DisplayName','Sprinkler Data Using Computed Torque Signal')
    title("ODE45 Soltuion Using Comptued Torque")
    xlabel("t")
    ylabel("$$\phi(t)$$")
    legend()
    hold off
end 
