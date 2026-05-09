close all; clear all; clc;
%% Description
% This script has for purpose to see the internal resistance of 4 tests 
%% 1. Configuration
% Prompt the user for the number of cells in the console
num_cells = input('Enter the number of cells to process: ');

% Define the minimum current step to be considered a pulse (in Amperes)
current_min_step = 0.5;

% Define the root folder 
base_dir = 'C:\ncr18650_identification\dataset-thermal\';


% Define the patterns for the 4 tests here
% For example: {'_HPPC_10C*.txt', '_HPPC_25C*.txt', ...}
test_patterns = {'_MPD*.txt', '_RPD*.txt', '_DPD*.txt', '_RSD*.txt'};

num_tests = length(test_patterns);

% Initialize the result matrix with NaNs (Not a Number)
% NaNs are used so if a file is missing, it doesn't default to 0 and ruin your stats.
% Rows = Cells, Columns = Tests
R0_matrix = NaN(num_cells, num_tests);

%% 2. Automated Processing Loop
for i = 1:num_cells
    
    % Create the cell identifier (i=1 -> 'BID001')
    cell_id = sprintf('BID%03d', i); 
    folder_path = fullfile(base_dir, cell_id);
    
    fprintf('\n--- Processing %s ---\n', cell_id);
    
    % Inner loop: Process the 4 tests for this specific cell
    for j = 1:num_tests
        
        % Build the filename pattern for this specific test
        pattern = fullfile(folder_path, [cell_id, test_patterns{j}]);
        files = dir(pattern);
        
        % Safety check: Verify if the file was found
        if isempty(files)
            fprintf('  -> Warning: Test %d not found. Skipping.\n', j);
            continue; % Leave the matrix value as NaN and skip to next test
        end
        
        % Get the full path of the found file
        file_path = fullfile(folder_path, files(1).name);
        
        % Call the extraction function 
        % (Using ~ to ignore the v_array output since we don't need it here)
        [r0_array, ~, i_array] = get_r0(file_path, current_min_step);
        
        % --- FILTERING STEP: Remove occurrences where current is 0 ---
        valid_idx = (i_array ~= 0);
        r0_array = r0_array(valid_idx);
        
        % Check if any valid pulses remain
        if isempty(r0_array)
            fprintf('  -> Warning: No valid pulses (I ~= 0) found in Test %d.\n', j);
        else
            % Calculate the average R0 and store it in our matrix
            R0_matrix(i, j) = mean(r0_array);
            fprintf('  -> Test %d processed successfully (Mean R0 = %.4f Ohms).\n', j, R0_matrix(i, j));
        end
    end
end

%% 3. Display Final Results
fprintf('\n=========================================================\n');
disp('FINAL RESULTS: Average Internal Resistance R0 (Ohms)');
fprintf('=========================================================\n');

% Create a clean MATLAB Table for display
% Row names will be BID001, BID002...
cell_names = arrayfun(@(x) sprintf('BID%03d', x), (1:num_cells)', 'UniformOutput', false);

% Column names will be Test_1, Test_2...
test_names = {'MPD', 'RPD', 'DPD', 'RSD'};

% Build and display the table
results_table = array2table(R0_matrix, 'RowNames', cell_names, 'VariableNames', test_names);
disp(results_table);