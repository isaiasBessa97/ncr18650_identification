close all; clear all; clc;

%% 1. Configuration
current_min_step = 0.5;
Qn = 3.05;         
initial_soc = 100; 
fall_delay_steps = 2; 

disp('Select the pulse test file...');
[file_name, path_name] = uigetfile({'*.txt;*.csv', 'Text/Data Files'; '*.*', 'All Files'}, 'Select the Pulse Test file');
if isequal(file_name,0); error('Cancelled by user'); end
file_path = fullfile(path_name, file_name);

%% 2. Extraction via the ALL-IN-ONE function
disp('Extracting R0 and SoC directly from the hybrid function...');
% We get soc_array directly from the function (no more '~')!
[soc_array, r0_array, v_array, i_array] = get_r0_rise_and_fall(file_path, Qn, initial_soc, current_min_step, fall_delay_steps);

%% 3. Rise (Pulse ON) and Fall (Pulse OFF) Classification
% i_array contains the current just BEFORE the transition.
% - If the current was close to 0, it's a current step up (Rise).
% - If the current was high, it's a current cut-off (Fall).
is_rise = abs(i_array) < 0.1; 
is_fall = ~is_rise;

% Cleanly separating the data
soc_rise = soc_array(is_rise);
r0_rise  = r0_array(is_rise);

soc_fall = soc_array(is_fall);
r0_fall  = r0_array(is_fall);

%% 4. Export to Excel
disp('Preparing data for export...');
[save_name, save_path] = uiputfile('*.xlsx', 'Save the Results as an Excel file');

if ~isequal(save_name, 0)
    full_save_path = fullfile(save_path, save_name);
    
    Table_Rise = table(soc_rise, r0_rise, 'VariableNames', {'SoC_Percentage', 'Internal_Resistance_R0_Ohms'});
    Table_Fall = table(soc_fall, r0_fall, 'VariableNames', {'SoC_Percentage', 'Internal_Resistance_R0_Ohms'});
    
    writetable(Table_Rise, full_save_path, 'Sheet', 'Pulse_ON_Rise');
    writetable(Table_Fall, full_save_path, 'Sheet', 'Pulse_OFF_Fall');
    
    fprintf('\nSUCCESS! Data saved to: %s\n', save_name);
else
    disp('Export cancelled by the user.');
end

%% 5. Plotting the Results
figure('Color', 'w'); hold on; grid on;

% Avoid warnings if one of the arrays is empty
if ~isempty(r0_rise), min_r0_rise = min(r0_rise); else, min_r0_rise = 0; end
if ~isempty(r0_fall), min_r0_fall = min(r0_fall); else, min_r0_fall = 0; end

% Plot Rise (Blue dots)
plot(soc_rise, r0_rise, 'ob', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'b', 'DisplayName', sprintf('Pulse ON (Rise: Immediate)| R0_min = %.4f', min_r0_rise));

% Plot Fall (Red squares)
plot(soc_fall, r0_fall, 'sr', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'r', 'DisplayName', sprintf('Pulse OFF (Fall: %d sec delay)| R0_min = %.4f', fall_delay_steps, min_r0_fall));

% Formatting
xlabel('State of Charge (SoC) [%]', 'FontWeight', 'bold');
ylabel('Internal Resistance R_0 (\Omega)', 'FontWeight', 'bold');
title('Internal Resistance vs SoC (Hybrid Method)', 'FontSize', 14);
legend('Location', 'best');

% Reverse X-axis to read from 100% down to 0%
set(gca, 'XDir', 'reverse'); 
hold off;
disp('Done!');