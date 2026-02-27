function [spin_switch, Re, trial] = read_name(data_name)
    % Remove extension first
    [~, name_no_ext, ~] = fileparts(data_name);
    
    % Split on underscores
    parts = strsplit(name_no_ext, '_');
    
    % Extract each part
    if parts{1} == "forward"
        spin_switch = 1;
    elseif parts{1} == "reverse"
        spin_switch = 0;
    else
        error("Invalid direction '%s' in file name - please use 'forward' or 'reverse'", parts{1});
    end

    Re     = str2double(parts{2});  % 500
    trial     = parts{3};         % "trial2"
end
