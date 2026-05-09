% close all; clear all; clc;

%% 1. Configuration
current_min_step = 0.5;
Qn = 3.05;         
initial_soc = 100; 
fall_delay_steps = 2; 

% File selection
disp('Select the pulse test file...');
[file_name, path_name] = uigetfile({'*.txt;*.csv', 'Text/Data Files'; '*.*', 'All Files'}, 'Select the Pulse Test file');
if isequal(file_name,0); error('Cancelled by user'); end
file_path = fullfile(path_name, file_name);

%% 2. Direct extraction using the new get_r0 function
disp('Extracting R0 and SoC directly from get_r0...');
% The function directly gives us the correct X (soc_fall) and Y (r0_fall)
[soc_fall, r0_fall, ~, ~] = get_r0(file_path, Qn, initial_soc, current_min_step, fall_delay_steps);

%% 3. Plotting the Results
figure('Color', 'w'); hold on; grid on;

min_r0_fall = min(r0_fall);
max_r0_fall = max(r0_fall);

% Plot Fall (Red squares)
plot(soc_fall, r0_fall, '-sr', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'r', 'DisplayName', sprintf('Pulse OFF (%d sec delay) | Min: %.4f | Max: %.4f', fall_delay_steps, min_r0_fall, max_r0_fall));

% Formatting
xlabel('State of Charge (SoC) [%]', 'FontWeight',