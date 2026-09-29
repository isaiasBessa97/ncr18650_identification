close all; clear all; clc;
%% Description
% Modified script to obtain a calculated pseudo-OCV based on the average of the
% charge and discharge of a cell at different temperatures (0.05 C-RATE).
% The preparation (start) and rest (end) phases are trimmed before calling get_ocv.

%% 1. Configuration
% Define the directories containing the charge and discharge files
charge_dir = 'C:\Users\PRH\Desktop\simulated_temperature_data\CCCV_tests';
discharge_dir = 'C:\Users\PRH\Desktop\simulated_temperature_data\CD_tests';

% Define the vector of temperatures to analyze (from 0 to 70 degrees in steps of 5)
temperatures = 0:5:70;

% Graph preparation
figure('Color', 'w'); hold on; grid on;
xlabel('State of Charge (SoC) [%]', 'FontWeight', 'bold');
ylabel('OCV (V)', 'FontWeight', 'bold');
title('Comparison of Pseudo-OCV curves by Temperature', 'FontSize', 14);

%% 2. Automated processing loop
for t = temperatures
    
    temp_str = sprintf('%03ddeg', t); 
    
    file_charge_name = sprintf('CCCV_005C_%s.csv', temp_str);
    file_discharge_name = sprintf('CDch_005C_%s.csv', temp_str);
    
    file_charge = fullfile(charge_dir, file_charge_name);
    file_discharge = fullfile(discharge_dir, file_discharge_name);
    
    if ~isfile(file_charge) || ~isfile(file_discharge)
        fprintf('Files for %d degrees not found. Moving to the next temperature.\n', t);
        continue; 
    end
    
    %% --- DATA READING AND TRIMMING (PREPARATION & REST) ---
    data_ch = readmatrix(file_charge, 'NumHeaderLines', 1);
    data_dis = readmatrix(file_discharge, 'NumHeaderLines', 1);
    
    % -- Charge Trimming (CCCV) --
    % 1. Find the start (when the current becomes negative)
    start_idx_ch = find(data_ch(:, 2) < -0.05, 1, 'first');
    if ~isempty(start_idx_ch)
        data_ch = data_ch(start_idx_ch:end, :);
        
        % 2. Find the end (when the test stops and the current falls back to 0)
        % We look for the first time the current rises above -0.001A
        end_idx_ch = find(data_ch(:, 2) > -0.001, 1, 'first');
        if ~isempty(end_idx_ch)
            data_ch = data_ch(1:end_idx_ch-1, :); % Trim the entire final rest period
        end
    end
    
    % -- Discharge Trimming (CDch) --
    % 1. Find the start (when the current becomes positive around 0.15A)
    start_idx_dis = find(data_dis(:, 2) > 0.05 & data_dis(:, 2) < 0.25, 1, 'first');
    if isempty(start_idx_dis)
        start_idx_dis = find(data_dis(:, 2) > 0.05, 1, 'first');
    end
    if ~isempty(start_idx_dis)
        data_dis = data_dis(start_idx_dis:end, :);
        
        % 2. Find the end (when the test stops and the current falls back to 0)
        % We look for the first time the current drops below 0.001A
        end_idx_dis = find(data_dis(:, 2) < 0.001, 1, 'first');
        if ~isempty(end_idx_dis)
            data_dis = data_dis(1:end_idx_dis-1, :); % Trim the entire final rest period
        end
    end
    
    %% --- CREATION OF TEMPORARY FILES FOR get_ocv.m ---
    temp_ch_file = fullfile(tempdir, 'temp_charge.txt');
    temp_dis_file = fullfile(tempdir, 'temp_discharge.txt');
    
    data_ch_formatted = [data_ch(:, 1), data_ch(:, 3), data_ch(:, 2), data_ch(:, 4)];
    data_dis_formatted = [data_dis(:, 1), data_dis(:, 3), data_dis(:, 2), data_dis(:, 4)];
    
    fid = fopen(temp_ch_file, 'w'); fprintf(fid, 'Time;Voltage;Current;Temp\n'); fclose(fid);
    writematrix(data_ch_formatted, temp_ch_file, 'Delimiter', ';', 'WriteMode', 'append');
    
    fid = fopen(temp_dis_file, 'w'); fprintf(fid, 'Time;Voltage;Current;Temp\n'); fclose(fid);
    writematrix(data_dis_formatted, temp_dis_file, 'Delimiter', ';', 'WriteMode', 'append');
    
    %% --- CALLING get_ocv.m AND PLOTTING ---
    [soc, ocv, qn] = get_ocv(temp_ch_file, temp_dis_file);
    
    delete(temp_ch_file); 
    delete(temp_dis_file);
    
    legend_str = sprintf('%d °C (Qn = %.3f Ah)', t, qn);
    plot(soc, ocv, 'LineWidth', 1.5, 'DisplayName', legend_str);
    
    fprintf('Temperature %d °C successfully processed!\n', t);
end

legend('Location', 'best');
hold off;