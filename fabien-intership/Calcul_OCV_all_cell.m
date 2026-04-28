close all; clear all; clc;
%% Description
% script developped to get a calculated pseudo OCV based on the average of the
% charge and discharge of the cell with a 0.05 C-RATE
% Input: The number of the cells observed with these previous conditions

%% 1. Configuration
% Prompt the user for the number of cells in the console
num_cells = input('Enter the number of cells to process: ');

% Define the root folder 
base_dir = 'C:\ncr18650_identification\dataset-thermal\';

% Prepare the plot
figure('Color', 'w'); hold on; grid on;
xlabel('State of Charge (SoC) [%]', 'FontWeight', 'bold');
ylabel('OCV (V)', 'FontWeight', 'bold');
title('Comparison of Pseudo-OCV Curves', 'FontSize', 14);

%% 2. Automated Processing Loop
for i = 1:num_cells
    
    % Create the folder name and identifier (i=1 -> 'BID001')
    % %03d means "integer with 3 digits, zero-padded"
    cell_id = sprintf('BID%03d', i); 
    folder_path = fullfile(base_dir, cell_id);
    
    % Build the filename pattern for 0.005 C-Rate
    % Ex: BID001_CCCV005*.txt
    charge_pattern = fullfile(folder_path, [cell_id, '_CCCV005*.txt']);
    discharge_pattern = fullfile(folder_path, [cell_id, '_CDch005*.txt']);
    
    % Use 'dir' to find the exact file matching the pattern
    files_ch = dir(charge_pattern);
    files_dis = dir(discharge_pattern);
    
    % Safety check: Verify if the files were actually found
    if isempty(files_ch) || isempty(files_dis)
        fprintf('Files for %s not found. Skipping to next cell.\n', cell_id);
        continue; % Skip directly to the next loop iteration
    end
    
    % Get the full path of the found files
    file_charge = fullfile(folder_path, files_ch(1).name);
    file_discharge = fullfile(folder_path, files_dis(1).name);
    
    % Call the created function calculating OCV
    [soc, ocv, qn] = get_ocv(file_charge, file_discharge);
    
    % Plot the curve with a dynamic legend including the measured capacity
    legend_str = sprintf('%s (Qn = %.3f Ah)', cell_id, qn);
    plot(soc, ocv, 'LineWidth', 1.5, 'DisplayName', legend_str);
    
    %console message to track progress
    fprintf('Cell %s processed successfully!\n', cell_id);
end

% Display the final legend
legend('Location', 'best');
hold off;