close all; clear all; clc;
%% Description
% This script allow to compute parameters
%% 1. Configuration
% Pulse detection threshold
current_min_step = 0.5;

% Battery specifications for SoC calculation
Qn = 3.05;         % Nominal capacity in Ah 
initial_soc = 100; % Starting SoC of the test (%)

% Prompt user for the dynamic pulse test file
disp('Select the pulse test file...');
[file_name, path_name] = uigetfile({'*.txt;*.csv', 'Text/Data Files'; '*.*', 'All Files'}, 'Select the Pulse Test file');
if isequal(file_name,0); error('Cancelled by user'); end
file_path = fullfile(path_name, file_name);

%% 2. Calculate Continuous SoC & Read Data
soc_full = get_soc(file_path, Qn, initial_soc);

% Need to extract Voltage and current for continued 
data = readmatrix(file_path, 'Delimiter', ';', 'NumHeaderLines', 1);
time = data(:, 1);
V = data(:, 2);
I = data(:, 3);

%% 3. Detect Pulses 
disp('Detecting pulses and relaxation periods...');

% State machine: 1 if current is ON (>0.1A), 0 if OFF
is_pulse_on = abs(I) > 0.1;

% Detect transitions: +1 means turning ON (Rise), -1 means turning OFF (Fall)
transitions = diff(is_pulse_on);

rise_idx = find(transitions == 1);
fall_idx = find(transitions == -1);

% Safety: Ensure the very first event is a RISE to keep arrays perfectly synchronized
if ~isempty(fall_idx) && ~isempty(rise_idx) && fall_idx(1) < rise_idx(1)
    fall_idx(1) = []; 
end

num_pulses = min(length(rise_idx), length(fall_idx));

% Initialize arrays to store our RC parameters
soc_results = zeros(num_pulses, 1);
R1_array = zeros(num_pulses, 1);
R2_array = zeros(num_pulses, 1);
C1_array = zeros(num_pulses, 1);
C2_array = zeros(num_pulses, 1);

%% 4. RC Parameters Identification (Strictly based on Article Equations)
disp('Identifying R1, C1, R2, C2 using explicit Article Equations...');

ft = fittype('y_e - x1_c*exp(-x/tau1) - x2_c*exp(-x/tau2)', ...
    'independent', 'x', ...
    'problem', 'y_e', ...
    'coefficients', {'tau1', 'tau2', 'x1_c', 'x2_c'});

% Configure fitting options with constraints
opts = fitoptions('Method', 'NonlinearLeastSquares', 'Display', 'Off');
opts.Lower = [0.1,  20,   0.01,   0.015];      % [tau1, tau2, x1_c, x2_c] lower bounds
opts.Upper = [1, 1000, 0.3, 0.3];      % [tau1, tau2, x1_c, x2_c] upper bounds
opts.StartPoint = [0.5, 100, 0.05, 0.05]; % Initial guess

valid_count = 0;

for k = 1:num_pulses
    ta_idx = rise_idx(k);   % Start of pulse (t_a)
    tc_idx = fall_idx(k);   % End of pulse / start of relaxation (t_c)
    
    % Find t_e (end of rest period just before the next pulse starts)
    if k < num_pulses
        te_idx = rise_idx(k+1) - 1;
    else
        te_idx = length(V);
    end
    
    if (te_idx - tc_idx) < 10
        continue;
    end
    
    % 4.1 Extract Known Variables for the specific pulse
    % ROBUST U: Get the true maximum current of this pulse
    u = max(abs(I(ta_idx+1 : tc_idx))); 
    t_on = time(tc_idx) - time(ta_idx); % Pulse duration (t_c - t_a)
    
    % Extract relaxation data
    t_relax = time(tc_idx+1 : te_idx) - time(tc_idx+1); 
    y_relax = V(tc_idx+1 : te_idx); 
    
    % y_e: Average of the last 10 points to avoid sensor noise 
    safe_points = min(10, length(y_relax));
    y_e = mean(y_relax(end-safe_points+1:end)); 
    
    % Force any noisy voltage spikes down to y_e so the fit doesn't face negative targets
    y_relax(y_relax > y_e) = y_e;
    
    % 4.2 Perform Curve Fitting 
    try
        [fitresult, ~] = fit(t_relax, y_relax, ft, opts, 'problem', y_e);
        
        tau1 = fitresult.tau1;
        tau2 = fitresult.tau2;
        x1_c = fitresult.x1_c;
        x2_c = fitresult.x2_c;
        
        % 4.3 Calculate RC Parameters
        R1 = x1_c / (u * (1 - exp(-t_on / tau1)));
        R2 = x2_c / (u * (1 - exp(-t_on / tau2)));
        C1 = tau1 / R1;
        C2 = tau2 / R2;
        
        % Save successful results
        valid_count = valid_count + 1;
        soc_results(valid_count) = soc_full(tc_idx);
        R1_array(valid_count) = R1;
        R2_array(valid_count) = R2;
        C1_array(valid_count) = C1;
        C2_array(valid_count) = C2;
        
    catch
        fprintf('Warning: Curve fitting failed for pulse %d. Skipping.\n', k);
    end
end

% Truncate unused pre-allocated space
soc_results = soc_results(1:valid_count);
R1_array = R1_array(1:valid_count);
R2_array = R2_array(1:valid_count);
C1_array = C1_array(1:valid_count);
C2_array = C2_array(1:valid_count);

disp('Identification complete!');

%% 5. Plotting the Parameters vs SoC
figure('Color', 'w', 'Position', [100, 100, 1000, 800]);

% R1
subplot(2, 2, 1); grid on; hold on;
plot(soc_results, R1_array, '-ob', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
xlabel('SoC [%]', 'FontWeight', 'bold'); ylabel('R_1 (\Omega)', 'FontWeight', 'bold');
title('Electrochemical Resistance (R_1)');
set(gca, 'XDir', 'reverse');

% R2
subplot(2, 2, 2); grid on; hold on;
plot(soc_results, R2_array, '-or', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
xlabel('SoC [%]', 'FontWeight', 'bold'); ylabel('R_2 (\Omega)', 'FontWeight', 'bold');
title('Concentration Resistance (R_2)');
set(gca, 'XDir', 'reverse');

% C1
subplot(2, 2, 3); grid on; hold on;
plot(soc_results, C1_array, '-og', 'LineWidth', 1.5, 'MarkerFaceColor', 'g');
xlabel('SoC [%]', 'FontWeight', 'bold'); ylabel('C_1 (F)', 'FontWeight', 'bold');
title('Electrochemical Capacitance (C_1)');
set(gca, 'XDir', 'reverse');

% C2
subplot(2, 2, 4); grid on; hold on;
plot(soc_results, C2_array, '-om', 'LineWidth', 1.5, 'MarkerFaceColor', 'm');
xlabel('SoC [%]', 'FontWeight', 'bold'); ylabel('C_2 (F)', 'FontWeight', 'bold');
title('Concentration Capacitance (C_2)');
set(gca, 'XDir', 'reverse');

sgtitle(sprintf('Equivalent Circuit RC Parameters vs SoC\nFile: %s', file_name), 'FontSize', 14, 'FontWeight', 'bold');