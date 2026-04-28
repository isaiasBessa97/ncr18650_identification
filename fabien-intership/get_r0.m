function [r0_values, v_values, i_values] = get_r0(file_path, min_step)
% GET_R0 Calculates internal resistance (R0), Voltage and Current at pulse
%
% Inputs:
%   file_path : Full path to the test file (string)
%   min_step  : Minimum current step to trigger the calculation (e.g., 0.5 A)
%
% Outputs:
%   r0_values : Array containing the calculated R0 values (in Ohms)
%   v_values  : Array containing the voltage at the pulse
%   i_values  : Array containing the current at the pulse

    %% 1. Read Data
    data = readmatrix(file_path, 'Delimiter', ';', 'NumHeaderLines', 1);
    
    % Extract Voltage (Column 2) and Current (Column 3)
    V = data(:, 2);
    I = data(:, 3); 
    
    %% 2. Detect Pulses (Delta I)
    dI = diff(I);
    dV = diff(V);
    
    % Find the exact row indices where the current step exceeds the step
    pulse_indices = find(abs(dI) >= min_step);
    
    %% 3. Calculate R0, V and I
    num_pulses = length(pulse_indices);
    r0_values = zeros(num_pulses, 1);
    v_values = zeros(num_pulses, 1);
    i_values = zeros(num_pulses, 1);
    
    for j = 1:num_pulses
        idx = pulse_indices(j);
        
        % Calculate R0 (dV / dI)
        r0_values(j) = abs(dV(idx) / dI(idx));
        
        % Store the Voltage and Current exactly at the new pulse state
        v_values(j) = V(idx+1);
        i_values(j) = I(idx+1);
    end
    
    %% 4. Data Cleaning (Remove initialization artifacts)
    valid_indices = r0_values < 0.5;
    r0_values = r0_values(valid_indices);
    v_values = v_values(valid_indices);
    i_values = i_values(valid_indices);
end
