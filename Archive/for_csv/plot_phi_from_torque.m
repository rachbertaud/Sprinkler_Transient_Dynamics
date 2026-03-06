function plot_phi_from_torque(full_x, full_y, phi_gen, spin_switch, trial_num, num_trials)
    persistent x_ref phi_gen_all

    trial_rgb    = {[252/255, 70/255, 170/255], ...   % hot pink
                    [227/255, 159/255, 246/255], ...  % soft purple
                    [100/255, 149/255, 237/255]};     % periwinkle blue
    trial_styles = {'-', '--', '-.'};
    gen_alpha    = 0.5;
    data_alpha   = 0.18;
    avg_color    = [110/255, 50/255, 180/255];  % deep violet

    % Reset accumulators on first trial
    if trial_num == 1
        x_ref        = full_x;
        phi_gen_all  = zeros(num_trials, length(full_x));
    end

    % Interpolate both signals onto reference grid
    gen_interp  = interp1(full_x, phi_gen, x_ref, 'pchip', NaN);
    data_interp = interp1(full_x, full_y,  x_ref, 'pchip', NaN);

    phi_gen_all(trial_num, :) = gen_interp;

    f = figure(4);
    theme(f, "light");
    hold on

    if trial_num == 1
        yline(0, '--k', 'LineWidth', 0.8, 'HandleVisibility', 'off')
    end

    % Original data — very faint background, no legend entry
    plot(x_ref, data_interp, trial_styles{trial_num}, 'LineWidth', 0.8, ...
        'Color', [252/255, 70/255, 170/255, data_alpha], 'HandleVisibility', 'off')

    % ODE45 reconstruction — faint, per-trial color and style
    plot(x_ref, gen_interp, trial_styles{trial_num}, 'LineWidth', 1.8, ...
        'Color', [trial_rgb{trial_num}, gen_alpha], ...
        'DisplayName', sprintf('Trial %d', trial_num))

    if trial_num == num_trials
        phi_gen_avg = mean(phi_gen_all, 1, 'omitnan');
        plot(x_ref, phi_gen_avg, '-', 'Color', avg_color, 'LineWidth', 2.5, ...
            'DisplayName', 'Mean')
        title('ODE45 Solution Using Recovered Torque')
        xlabel('$t$')
        ylabel('$\phi(t)$')
        xlim([x_ref(1), x_ref(end)])
        legend('Location', 'northeast')
    end

    hold off
end
