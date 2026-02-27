function plot_torque(franken_x, signal, N_f, trial_num, num_trials)
    persistent x_ref sig_all cum_all

    trial_rgb    = {[252/255, 70/255, 170/255], ...   % hot pink
                    [227/255, 159/255, 246/255], ...  % soft purple
                    [100/255, 149/255, 237/255]};     % periwinkle blue
    trial_styles = {'-', '--', '-.'};
    line_alpha   = 0.4;
    avg_color    = [110/255, 50/255, 180/255];  % deep violet

    x_pos   = franken_x(N_f/2:end);
    sig_pos = signal(N_f/2:end);

    % Reset accumulators on first trial
    if trial_num == 1
        x_ref   = x_pos;
        sig_all = zeros(num_trials, length(x_pos));
        cum_all = zeros(num_trials, length(x_pos));
    end

    % Interpolate onto reference grid
    sig_interp = interp1(x_pos, sig_pos, x_ref, 'pchip', 0);
    cum_interp = cumtrapz(x_ref, sig_interp);

    sig_all(trial_num, :) = sig_interp;
    cum_all(trial_num, :) = cum_interp;

    % --- Torque signal ---
    f = figure(2);
    theme(f, "light");
    hold on

    if trial_num == 1
        yline(0, '--k', 'LineWidth', 0.8, 'HandleVisibility', 'off')
    end

    plot(x_ref, sig_interp, trial_styles{trial_num}, 'LineWidth', 1.5, ...
        'Color', [trial_rgb{trial_num}, line_alpha], ...
        'DisplayName', sprintf('Trial %d', trial_num))

    if trial_num == num_trials
        sig_avg = mean(sig_all, 1);
        plot(x_ref, sig_avg, '-', 'Color', avg_color, 'LineWidth', 2.5, ...
            'DisplayName', 'Mean')
        title('Recovered Torque Signal (Inverse FFT)')
        xlabel('$t$')
        ylabel('$\tau(t)$')

        % Clip y to 99th percentile across all trials to suppress FFT edge spikes

        legend('Location', 'northeast')
    end

    hold off

    % --- Cumulative integral ---
    f = figure(5);
    theme(f, "light");
    hold on

    if trial_num == 1
        yline(0, '--k', 'LineWidth', 0.8, 'HandleVisibility', 'off')
    end

    plot(x_ref, cum_interp, trial_styles{trial_num}, 'LineWidth', 1.5, ...
        'Color', [trial_rgb{trial_num}, line_alpha], ...
        'DisplayName', sprintf('Trial %d', trial_num))

    if trial_num == num_trials
        cum_avg = mean(cum_all, 1);
        [max_trap, idx_max] = max(cum_avg);
        p = plot(x_ref, cum_avg, '-', 'Color', avg_color, 'LineWidth', 2.5, ...
            'DisplayName', 'Mean');
        datatip(p, x_ref(idx_max), max_trap, 'Location', 'northwest');
        xlabel('$t$')
        title('Cumulative Integral of Torque Signal')
        ylabel('$\int_0^t \tau \, dt$')
        legend('Location', 'northeast')
    end

    hold off
end
