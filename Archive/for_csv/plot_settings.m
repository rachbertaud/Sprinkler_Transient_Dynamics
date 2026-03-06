function plot_settings(plot_switch)
    if(plot_switch == 1)
        fs = 20;  % set your preferred font size here

        % plotting settings
        set(groot, 'defaultTextInterpreter', 'latex');
        set(groot, 'defaultLegendInterpreter', 'latex');
        set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
        set(groot, 'defaultAxesFontSize', fs);
        set(groot, 'defaultTextFontSize', fs);
        set(groot, 'defaultLegendFontSize', 30);
        
        % Optional: line width and marker size defaults
        set(groot, 'defaultLineLineWidth', 3);
        set(groot, 'defaultLineMarkerSize', 15);
    else 
        plot_switch = 0;
    end 
end 
