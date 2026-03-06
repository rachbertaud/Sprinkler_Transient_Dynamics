function plot_processed_data(franken_x, franken_y, trial_num, num_trials)
    persistent x_ref y_pos_all

    trial_rgb    = {[252/255, 70/255, 170/255], ...   % hot pink
                    [227/255, 159/255, 246/255], ...  % soft purple
                    [100/255, 149/255, 237/255]};     % periwinkle blue
    trial_styles = {'-', '--', '-.'};
    line_alpha   = 0.4;
    avg_color    = [110/255, 50/255, 180/255];  % deep violet

    pos_mask = franken_x >= 0;
    x_pos    = franken_x(pos_mask);
    y_pos    = franken_y(pos_mask);

    % Reset accumulators on first trial
    if trial_num == 1
        x_ref     = x_pos;
        y_pos_all = zeros(num_trials, length(x_pos));
    end

    % Interpolate positive region onto reference grid
    y_interp = interp1(x_pos, y_pos, x_ref, 'pchip', NaN);
    y_pos_all(trial_num, :) = y_interp;

    f = figure(3);
    theme(f, "light");
    hold on

    % Zero-padded region — light shaded patch drawn first so it sits behind all data
    if trial_num == 1
        patch([franken_x(1), 0, 0, franken_x(1)], [-1e4, -1e4, 1e4, 1e4], ...
            [0.91, 0.89, 0.96], 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
            'DisplayName', 'Zero-Padded Region')
        xline(0, '-', 'Color', [0.60, 0.50, 0.70], 'LineWidth', 1.5, 'HandleVisibility', 'off')
        yline(0, '--k', 'LineWidth', 0.5, 'HandleVisibility', 'off')
    end

    % Real data + analytical tail — faint, per-trial color and style
    plot(x_ref, y_interp, trial_styles{trial_num}, 'LineWidth', 1.5, ...
        'Color', [trial_rgb{trial_num}, line_alpha], ...
        'DisplayName', sprintf('Trial %d', trial_num))

    if trial_num == num_trials
        y_avg = mean(y_pos_all, 1, 'omitnan');
        plot(x_ref, y_avg, '-', 'Color', avg_color, 'LineWidth', 2.5, ...
            'DisplayName', 'Mean')

        % Clip y to data range so the background patch doesn't expand the axes
        data_max = max(abs(y_pos_all(:)));
        xlim([franken_x(1), franken_x(end)])
        ylim([-data_max * 1.15, data_max * 1.15])
        xlabel('$t$')
        ylabel('$\phi(t)$')
        legend('Location', 'northeast')
        title('Combined Signal Used in FFT')
    end

    hold off
end
