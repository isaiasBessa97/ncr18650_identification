% clear; clc; close all;

%% 1. Configuration and Loading
[f_ch, p_ch] = uigetfile('*.txt', 'Select the CHARGE file (0.05C)');
[f_dis, p_dis] = uigetfile('*.txt', 'Select the DISCHARGE file (0.05C)');
[f_test, p_test] = uigetfile('*.txt', 'Select the TEST file (MPD)');

file_charge = fullfile(p_ch, f_ch);
file_discharge = fullfile(p_dis, f_dis);
file_test = fullfile(p_test, f_test);

Qn = 3.05; % Nominal Capacity in Ah
initial_soc = 100; % Starting SoC in %

% Extract mean OCV curve
[soc_ocv, V_ocv_raw, ~] = get_ocv(file_charge, file_discharge);

% Clean and create OCV Polynomial (degree 9)
valid_idx = ~isnan(V_ocv_raw);
soc_valid = soc_ocv(valid_idx);
V_valid = V_ocv_raw(valid_idx);
p_coeffs_ocv = polyfit(soc_valid, V_valid, 9);

% True SoC from the test file 
soc_true = get_soc(file_test, Qn, initial_soc);
if max(soc_true) <= 1.05
    soc_true = soc_true * 100;
end

% Load test data 
data = readmatrix(file_test, 'Delimiter', ';', 'NumHeaderLines', 1);
time = data(:, 1);
V_meas = data(:, 2);
I_meas = data(:, 3); % I > 0 means discharge

N = length(time);
Ts = 1; % Sampling time in seconds
 

%% 2. RLS Initialization
lambda = 0.9999; % Forgetting factor
P = 1 * eye(5); % Covariance matrix
theta = [0.1; 0.1; 0.01; 0.01; 0.01]; % Initial weights [a1, a2, b0, b1, b2]

% BMS Internal States & Storage
soc_estimated = zeros(N, 1);
soc_estimated(1) = 100; % MODIFY HERE THE INITIAL SOC OF THE MODEL

theta_history = zeros(N, 5);
R0_hist = zeros(N, 1); R1_hist = zeros(N, 1); C1_hist = zeros(N, 1);
R2_hist = zeros(N, 1); C2_hist = zeros(N, 1);
V_model = zeros(N, 1);

% Past state variables (y(k-1), y(k-2), u(k-1), u(k-2))
y_past = [0; 0];
u_past = [0; 0];

%% 3. Real-Time Loop (The Core BMS Algorithm)
disp('Starting Real-Time RLS with Coulomb Counting...');

for k = 1:N
    % A. Read Sensors at step k (The only things the BMS actually "knows")
    u_k = I_meas(k); 
    
    % B. BMS SoC Estimation (Coulomb Counting)
    if k > 1
        % I > 0 is discharge, so we subtract from the estimated SoC
        % Ts is in seconds, so we divide by 3600 to get Ah
        soc_estimated(k) = soc_estimated(k-1) - (u_k * Ts / (Qn * 3600)) * 100;
        %SOC saturation to not go over 100 and under 0 to avoid divergence
        soc_estimated(k) = max(0, min(100, soc_estimated(k)));
    end
    
    % C. OCV Estimation based on Estimated SoC
    ocv_k = polyval(p_coeffs_ocv, soc_estimated(k));
    
    % D. Pure dynamics calculation (Observation variable y_k)
    y_k = ocv_k - V_meas(k); 
    
    % --- START OF RLS ---
    if k > 2
        % Observation vector Phi
        phi_k = [y_past(1); y_past(2); u_k; u_past(1); u_past(2)];
        
        % Predict model voltage BEFORE update
        y_pred = phi_k' * theta;
        V_model(k) = ocv_k - y_pred;
        
        % Update ONLY if in dynamic current phase or relaxing
        if abs(u_k) > 0.05 || abs(u_k - u_past(1)) > 0.05
            [theta, P] = rls_step(y_k, phi_k, theta, P, lambda);
        end
    else
        V_model(k) = V_meas(k); % Model follows measurement for the first 2 steps
    end
    % --- END OF RLS ---
    
    % E. Save mathematical parameters
    theta_history(k, :) = theta';
    
    % F. Convert to 2RC physical quantities
    [r0, r1, c1, r2, c2] = theta_to_2rc(theta, Ts);
    R0_hist(k) = r0; R1_hist(k) = r1; C1_hist(k) = c1;
    R2_hist(k) = r2; C2_hist(k) = c2;

    % G. Time shift for variables (k becomes k-1)
    y_past = [y_k; y_past(1)];
    u_past = [u_k; u_past(1)];
end



%% 4. Results & Validation
% Ignore the first few seconds (RLS convergence time)
valid_idx = time > 10; 

% Errors calculation (mV)
erreur_abs_V = V_meas(valid_idx) - V_model(valid_idx);
erreur_moyenne_V = mean(erreur_abs_V);
rmse_V = sqrt(mean((V_meas(valid_idx) - V_model(valid_idx)).^2));
moy_R0 = mean(R0_hist(valid_idx));

% --- Console Output ---
fprintf('\n====================================\n');
fprintf('   RLS 2RC MODEL PERFORMANCES       \n');
fprintf('====================================\n');
fprintf('Mean Absolute Error : %.4f V\n', erreur_moyenne_V);
fprintf('RMSE                : %.4f V\n', rmse_V);
fprintf('====================================\n\n');

% --- Plot 1: Voltage Validation ---
figure('Name', 'Real-Time RLS - Voltage Validation', 'Position', [100, 100, 1000, 500]);
plot(time, V_meas, 'k', 'LineWidth', 1.5); hold on;
plot(time, V_model, 'r--', 'LineWidth', 1.5);
ylabel('Voltage (V)', 'FontWeight', 'bold'); 
xlabel('Time (s)', 'FontWeight', 'bold');
title(sprintf('Measured vs Model 2RC (RMSE: %.4f V)', rmse_V), 'FontSize', 12);
legend('V measured (Experimental)', 'V model (RLS Estimation)', 'Location', 'best'); 
grid on;

% --- Plot 2: SoC True vs SoC Estimated ---
figure('Name', 'SoC Validation', 'Position', [150, 150, 1000, 400]);
plot(time, soc_true, 'k', 'LineWidth', 1.5); hold on;
plot(time, soc_estimated, 'b--', 'LineWidth', 1.5);
ylabel('State of Charge (%)', 'FontWeight', 'bold'); 
xlabel('Time (s)', 'FontWeight', 'bold');
title('Ground Truth SoC vs BMS Estimated SoC', 'FontSize', 12);
legend('True SoC', 'Estimated SoC (Coulomb Counting)', 'Location', 'best'); 
grid on;

% --- Plot 3: 2RC Parameters Evolution ---
figure('Name', '2RC Parameters vs Estimated SoC', 'Position', [200, 200, 800, 800]);

% Remove the massive starting peak (first 50s)
valid_p = time > 50;

% Top graph: Resistances
subplot(2,1,1); hold on; grid on;
scatter(soc_estimated(valid_p), R0_hist(valid_p), 15, 'k', 'filled', 'DisplayName', 'R_0 (Internal)');
scatter(soc_estimated(valid_p), R1_hist(valid_p), 15, 'b', 'filled', 'DisplayName', 'R_1 (Fast)');
scatter(soc_estimated(valid_p), R2_hist(valid_p), 15, 'r', 'filled', 'DisplayName', 'R_2 (Slow)');
xlabel('Estimated State of Charge (%)', 'FontWeight', 'bold'); 
ylabel('Resistance (\Omega)', 'FontWeight', 'bold');
title('Real-Time Estimated Resistances');
legend('Location', 'best');
ylim([0, 0.5]); 
set(gca, 'XDir', 'reverse'); 

% Bottom graph: Capacitances
subplot(2,1,2); hold on; grid on;
scatter(soc_estimated(valid_p), C1_hist(valid_p), 15, 'g', 'filled', 'DisplayName', 'C_1 (Fast)');
scatter(soc_estimated(valid_p), C2_hist(valid_p), 15, 'm', 'filled', 'DisplayName', 'C_2 (Slow)');
xlabel('Estimated State of Charge (%)', 'FontWeight', 'bold'); 
ylabel('Capacitance (F)', 'FontWeight', 'bold');
title('Real-Time Estimated Capacitances');
legend('Location', 'best');
ylim([0, 80000]); 
set(gca, 'XDir', 'reverse');