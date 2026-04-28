% close all; clear all; clc;

%% 1. Configuration
% Define the minimum current step to be considered a pulse (in Amperes)
current_min_step = 0.5;

% Prompt user for the FIRST dynamic pulse test file
disp('Select the FIRST pulse test file...');
[file_name1, path_name1] = uigetfile({'*.txt;*.csv', 'Text/Data Files'; '*.*', 'All Files'}, 'Select the FIRST Pulse Test file');
if isequal(file_name1,0); error('Cancelled by user'); end
file_path1 = fullfile(path_name1, file_name1);

% Prompt user for the SECOND dynamic pulse test file
disp('Select the SECOND pulse test file...');
[file_name2, path_name2] = uigetfile({'*.txt;*.csv', 'Text/Data Files'; '*.*', 'All Files'}, 'Select the SECOND Pulse Test file');
if isequal(file_name2,0); error('Cancelled by user'); end
file_path2 = fullfile(path_name2, file_name2);

%% 2. Calculate R0 and Extract V, I
% Call the extraction function for both files (now returns 3 variables)
[r0_array1, v_array1, i_array1] = get_r0(file_path1, current_min_step);
[r0_array2, v_array2, i_array2] = get_r0(file_path2, current_min_step);

% Check if any pulses were found
if isempty(r0_array1) || isempty(r0_array2)
    error('One or both files did not contain current steps greater than %.2f A.', current_min_step);
end

% Create the occurrence axis (1, 2, 3, 4...) for each test
occurrence1 = 1:length(r0_array1);
occurrence2 = 1:length(r0_array2);

%% 3. Calculate Averages, Differences & Export Data
% To subtract, arrays must be the exact same length. 
min_len = min(length(r0_array1), length(r0_array2));

% Sync all arrays to the minimum length
r0_1_sync = r0_array1(1:min_len);
r0_2_sync = r0_array2(1:min_len);
v1_sync   = v_array1(1:min_len);
i1_sync   = i_array1(1:min_len);
v2_sync   = v_array2(1:min_len);
i2_sync   = i_array2(1:min_len);

diff_r0 = r0_1_sync - r0_2_sync; % Point-by-point difference

% Calculate the global averages (Static R0)
mean_r0_1 = mean(r0_array1);
mean_r0_2 = mean(r0_array2);
static_diff = mean_r0_1 - mean_r0_2; 

% Prompt user to choose where to save the text file
disp('Choose where to save the differences file...');
[save_name, save_path] = uiputfile('*.txt', 'Save the Differences as');

if ~isequal(save_name, 0)
    % Open the file for writing ('w')
    fid = fopen(fullfile(save_path, save_name), 'w');
    
    % Write the NEW header row
    fprintf(fid, 'Pulse_Occurrence;Difference_R0;V_Test1;I_Test1;V_Test2;I_Test2\n');
    
    % Loop through ALL aligned points
    for k = 1:min_len
        fprintf(fid, '%d;%.6f;%.4f;%.4f;%.4f;%.4f\n', ...
            k, diff_r0(k), v1_sync(k), i1_sync(k), v2_sync(k), i2_sync(k));
    end
    
    % Write the static averages at the very bottom of the file (adapted for new columns)
    fprintf(fid, '\nSTATIC_DIFF_AVERAGE;%.6f;-;-;-;-\n', static_diff);
    
    % Always close the file
    fclose(fid);
    fprintf('Data successfully saved to: %s\n', save_name);
else
    disp('Export cancelled by user. Continuing with plotting...');
end

%% 4. Plotting
figure('Color', 'w'); hold on; grid on;

% Plot R0 vs Occurrence for Test 1 (Blue circles)
plot(occurrence1, r0_array1, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'b', 'Color', 'b', 'DisplayName', ['Test 1: ' file_name1]);

% Plot R0 vs Occurrence for Test 2 (Red squares)
plot(occurrence2, r0_array2, '-s', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'r', 'Color', 'r', 'DisplayName', ['Test 2: ' file_name2]);

% Plot the Average Lines (Horizontal dashed lines)
yline(mean_r0_1, '--b', 'LineWidth', 1.5, 'DisplayName', sprintf('Mean Test 1: %.4f \\Omega', mean_r0_1));
yline(mean_r0_2, '--r', 'LineWidth', 1.5, 'DisplayName', sprintf('Mean Test 2: %.4f \\Omega', mean_r0_2));

xlabel('Pulse Occurrence (#)', 'FontWeight', 'bold');
ylabel('Internal Resistance R_0 (\Omega)', 'FontWeight', 'bold');

% Add the static difference directly in the title
title_str = sprintf('Comparison of R_0 (Static Difference = %+.4f \\Omega)', static_diff);
title(title_str, 'FontSize', 14);

% Remove underscores in legend names to avoid subscript formatting
legend('Interpreter', 'none', 'Location', 'best');

% Set axis limits dynamically based on the longest test
max_occurrence = max(length(r0_array1), length(r0_array2));
xlim([0, max_occurrence + 1]);

hold off;