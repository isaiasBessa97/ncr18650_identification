function [soc_values, r0_values, v_values, i_values] = get_r0_rise_and_fall(file_path, Qn, initial_soc, min_step, delay_steps)
% GET_R0 Calculates R0 for ALL pulses (Rise and Fall).
%
% Inputs:
%   file_path   : File path
%   Qn          : Capacity (Ah)
%   initial_soc : Starting SoC (%)
%   min_step    : Detection threshold (e.g., 0.5 A)
%   delay_steps : Delay for the FALL only (e.g., 2)

    %% 1. Data Preparation
    soc_full = get_soc(file_path, Qn, initial_soc);
    data = readmatrix(file_path, 'Delimiter', ';', 'NumHeaderLines', 1);
    V = data(:, 2);
    I = data(:, 3); 
    
    %% 2. Transition Detection
    is_pulse_on = abs(I) > min_step;
    transitions = diff(is_pulse_on); 
    
    % Get all indices where the current changes (1=Rise, -1=Fall)
    all_trans_idx = find(transitions ~= 0);
    num_trans = length(all_trans_idx);
    
    %% 3. Initialization of output arrays
    r0_values = zeros(num_trans, 1);
    soc_values = zeros(num_trans, 1);
    v_values = zeros(num_trans, 1);
    i_values = zeros(num_trans, 1);
    
    valid_count = 0;

    %% 4. Hybrid calculation loop
    for j = 1:num_trans
        idx = all_trans_idx(j);
        trans_type = transitions(idx); % 1 for Rise, -1 for Fall
        
        if trans_type == 1 % --- RISE CASE (Pulse ON) ---
            % Immediate calculation (no delay needed)
            dV = abs(V(idx+1) - V(idx));
            dI = abs(I(idx+1) - I(idx));
            
        else % --- FALL CASE (Pulse OFF) ---
            % Apply the safety delay of X seconds
            target_idx = min(idx + delay_steps, length(V));
            dV = abs(V(target_idx) - V(idx));
            dI = abs(I(target_idx) - I(idx));
        end
        
        % Protection against division by zero and calculation of R0
        if dI > 0.1
            valid_count = valid_count + 1;
            r0_values(valid_count) = dV / dI;
            soc_values(valid_count) = soc_full(idx);
            v_values(valid_count) = V(idx);
            i_values(valid_count) = I(idx);
        end
    end
    
    % Final cleaning of arrays
    soc_values = soc_values(1:valid_count);
    r0_values = r0_values(1:valid_count);
    v_values = v_values(1:valid_count);
    i_values = i_values(1:valid_count);
    
    % Safety: Remove physically aberrant values (> 0.5 Ohm)
    % Safety: Remove physically aberrant values (> 0.5 Ohm)
    outlier_idx = r0_values > 0.5;
    soc_values(outlier_idx) = [];
    r0_values(outlier_idx) = [];
    v_values(outlier_idx) = []; 
    i_values(outlier_idx) = [];  
end