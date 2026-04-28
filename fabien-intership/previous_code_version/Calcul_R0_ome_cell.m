close all; clear all; clc;

%% 1. Configuration
% Define the minimum current step to be considered a pulse (in Amperes)
current_min_step = 0.5;

% Prompt user for the dynamic pulse test file
disp('Select the pulse test file...');
[file_name, path_name] = uigetfile({'*.txt;*.csv', 'Text/Data Files'; '*.*', 'All Files'}, 'Select the Pulse Test file');
if isequal(file_name,0); error('Cancelled by user'); end

file_path = fullfile(path_name, file_name);

%% 2. Calculate R0
% Call the extraction function
r0_array = get_r0(file_path, current_min_step);

% Check if any pulses were found
if isempty(r0_array)
    error('No current steps greater than %.2f A were found in this file.', current_min_step);
end

% Create the occurrence axis (1, 2, 3, 4...)
occurrence = 1:length(r0_array);

%% 3. Plotting
figure('Color', 'w'); hold on; grid on;

% Plot R0 vs Occurrence with lines and circular markers
plot(occurrence, r0_array, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'b', 'DisplayName', 'Internal Resistance R_0');

xlabel('Pulse Occurrence (#)', 'FontWeight', 'bold');
ylabel('Internal Resistance R_0 (\Omega)', 'FontWeight', 'bold');
title('Evolution of R_0 During Pulse Test', 'FontSize', 14);
legend('Location', 'best');

% Optional: Set axis limits to make the plot look cleaner
xlim([0, length(r0_array) + 1]);

hold off;