close all; clear all; clc;

%% 1. Configuration and Test Selection
disp('--- Internal Resistance (R0) vs SoC Analysis (Pulse OFF only) ---');
disp('Available tests:');
disp('1 - MPD (e.g., MPDch_000deg.csv)');
disp('2 - DST (e.g., DST_000deg.csv)');
% disp('3 - New Test (e.g., NEW_000deg.csv)');

test_choice = input('Enter the test number: ');

% Use switch to easily add new tests
switch test_choice
    case 1
        test_prefix = 'MPDch';
        base_dir = 'C:\Users\PRH\Desktop\simulated_temperature_data\MPD_tests'; 
    case 2
        test_prefix = 'DST';
        base_dir = 'C:\Users\PRH\Desktop\simulated_temperature_data\DST_tests'; 
        
    % case 3
    %     test_prefix = 'NOUVEAU'; % Le préfixe du nom de vos fichiers
    %     base_dir = 'C:\chemin\vers\votre\dossier_NOUVEAU\'; 
    
    otherwise
        error('Ivalid choice. Please reload the script.');
end

% Function parameters
current_min_step = 0.5;
Qn = 3.05;         
initial_soc = 100; 
fall_delay_steps = 2; 

% Array of temperatures
temperatures = [0, 10, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65];

%% 2. Chart Preparation
figure('Color', 'w', 'Position', [100, 100, 1000, 600]); 
hold on; grid on;
xlabel('State of Charge (SoC) [%]', 'FontWeight', 'bold');
ylabel('Internal Resistance R_0 (\Omega)', 'FontWeight', 'bold');
title(sprintf('Internal Resistance vs SoC (Pulse OFF) - Test : %s', test_prefix), 'FontSize', 14);

% Generate a dynamic color palette
colors = lines(length(temperatures)); 

%% 3. Processing Loop for Temperatures
for k = 1:length(temperatures)
    t = temperatures(k);
    
    % Build the file name
    temp_str = sprintf('%03ddeg', t);
    file_name = sprintf('%s_%s.csv', test_prefix, temp_str);
    file_path = fullfile(base_dir, file_name);
    
    % Check if the file exists
    if ~isfile(file_path)
        fprintf('File not found for %d °C, moving to the next.\n', t);
        continue;
    end
    
    fprintf('Processing %d °C...\n', t);
    
    %% --- DATA ADAPTATION WITHOUT MODIFYING YOUR FUNCTIONS ---
    % 1. Robust CSV reading (, or ;)
    try
        data = readmatrix(file_path, 'NumHeaderLines', 1);
    catch
        opts = detectImportOptions(file_path);
        data = readmatrix(file_path, opts);
    end
    
    % 2. Column reorganization: [Time(1), Voltage(3), Current(2), Temp(4)]
    if size(data, 2) >= 4
        data_formatted = [data(:, 1), data(:, 3), data(:, 2), data(:, 4)];
    else
        data_formatted = [data(:, 1), data(:, 3), data(:, 2)]; % Fallback if no Temp column
    end
    
    % 3. Creation of the temporary .txt file (with the correct delimiter ;)
    temp_file = fullfile(tempdir, 'temp_pulse_fall_test.txt');
    fid = fopen(temp_file, 'w'); fprintf(fid, 'Time;Voltage;Current;Temp\n'); fclose(fid);
    writematrix(data_formatted, temp_file, 'Delimiter', ';', 'WriteMode', 'append');
    
    %% --- EXECUTION OF THE CALCULATION ---
    [soc_fall, r0_fall, ~, ~] = get_r0(temp_file, Qn, initial_soc, current_min_step, fall_delay_steps);
    
    delete(temp_file);
    
    %% --- PLOTTING ---
    if ~isempty(soc_fall)
        min_r0_fall = min(r0_fall);
        max_r0_fall = max(r0_fall);
        
        % Plot with squares (-s to connect the points as in your script)
        plot(soc_fall, r0_fall, '-s', 'LineWidth', 1.5, 'MarkerSize', 6, ...
            'Color', colors(k,:), 'MarkerFaceColor', colors(k,:), ...
            'DisplayName', sprintf('%d °C | Min: %.4f', t, min_r0_fall));
    end
end

%% 4. Final Formatting
xlim([0 100]); % Force display strictly between 0 and 100%
set(gca, 'XDir', 'reverse'); % Reverse the X-axis to read from 100% to 0%
legend('Location', 'bestoutside', 'NumColumns', 1);
hold off;
disp('Analysis complete!');