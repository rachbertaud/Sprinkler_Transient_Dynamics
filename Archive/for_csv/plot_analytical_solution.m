function plot_analytical_solution(full_x, full_y, peak_index, x_peaks, y_peaks, phi_an, spin_switch, trial_num, Re, num_trials)
    persistent x_ref phi_all t0_vals

    trial_rgb    = {[252/255, 70/255, 170/255], ...   % hot pink
                    [227/255, 159/255, 246/255], ...  % soft purple
                    [100/255, 149/255, 237/255]};     % periwinkle blue
    trial_styles = {'-', '--', '-.'};
    fit_alpha    = 0.45;
    data_alpha   = 0.18;
    avg_color    = [110/255, 50/255, 180/255];  % deep violet

    % Reset accumulators on first trial
    if trial_num == 1
        x_ref   = full_x;
        phi_all = zeros(num_trials, length(full_x));
        t0_vals = zeros(1, num_trials);
    end

    % Evaluate fit on reference grid (phi_an is a function handle, works at any x)
    phi_vals  = phi_an(x_ref);
    data_vals = interp1(full_x, full_y, x_ref, 'pchip', NaN);

    phi_all(trial_num, :) = phi_vals;
    t0_vals(trial_num)    = x_peaks(peak_index);

    f = figure(1);
    if trial_num == 1
        theme(f, "light");
    end
    hold on

    % Raw data — very faint background, no legend entry
    plot(x_ref, data_vals, trial_styles{trial_num}, 'LineWidth', 0.8, ...
        'Color', [trial_rgb{trial_num}, data_alpha], 'HandleVisibility', 'off')

    % Analytical fit — faint, per-trial color and style
    plot(x_ref, phi_vals, trial_styles{trial_num}, 'LineWidth', 1.8, ...
        'Color', [trial_rgb{trial_num}, fit_alpha], ...
        'DisplayName', sprintf('Trial %d', trial_num))

    % Cut-point marker — faint dotted line, no legend entry
    xline(x_peaks(peak_index), ':', 'Color', trial_rgb{trial_num}, ...
        'LineWidth', 1.2, 'Alpha', fit_alpha, 'HandleVisibility', 'off')

    if trial_num == num_trials
        phi_avg = mean(phi_all, 1);
        plot(x_ref, phi_avg, '-', 'Color', avg_color, 'LineWidth', 2.5, ...
            'DisplayName', 'Mean')
        xline(mean(t0_vals), '--', 'Color', avg_color, ...
            'LineWidth', 1.5, 'HandleVisibility', 'off')
        yline(0, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off')

        all_vals = phi_all(:);
        ymargin  = 0.1 * (max(all_vals) - min(all_vals));
        ylim([min(all_vals) - ymargin, max(all_vals) + ymargin])

        xlabel('$t$')
        ylabel('$\phi(t)$')
        if spin_switch == 1
            title(sprintf('Forward, Re = %d, Analytical Fits', Re))
        else
            title(sprintf('Reverse, Re = %d, Analytical Fits', Re))
        end
        legend('Location', 'northeast')
    end

    hold off
end
