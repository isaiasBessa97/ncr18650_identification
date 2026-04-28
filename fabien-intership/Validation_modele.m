close all; clear all; clc;

%% 1. Cell Parameter Configuration
% These values define the physical properties of your Panasonic cell.
Qn_nominal = 3.05;      % Nominal capacity (Ah)
initial_soc = 100;      % Initial SoC for a test starting at full charge (%)
current_step = 0.5;     % Threshold to detect a current step (A)

%% 2. STEP 1: OCV (Open Circuit Voltage) Identification
% This step merges a slow charge and discharge test (0.05C) 
% to obtain the resting voltage curve "pseudo-OCV".
disp('--- STEP 1: OCV Modeling ---');
[f_ch, p_ch] = uigetfile('*.txt', 'Select the CHARGE file (0.05C)');
[f_dis, p_dis] = uigetfile('*.txt', 'Select the DISCHARGE file (0.05C)');

if isequal(f_ch,0) || isequal(f_dis,0); error('Missing OCV files'); end

% Call your get_ocv function
[soc_ocv, V_ocv_raw, ~] = get_ocv(fullfile(p_ch, f_ch), fullfile(p_dis, f_dis));

% OCV Polynomial Calculation (Least Squares)
% The article recommends a 9th degree to capture the complex chemistry.
poly_degree = 9; 
valid_idx = ~isnan(V_ocv_raw); % Filter to ignore empty edges (NaNs)
p_coeffs_ocv = polyfit(soc_ocv(valid_idx), V_ocv_raw(valid_idx), poly_degree);

fprintf('OCV Model (Degree %d) calculated.\n', poly_degree);

%% 3. STEP 2: Extraction of R0 and RC parameters (R1, C1, R2, C2)
% Here we use a dynamic dataset (like MPDch or RDch).
disp('--- STEP 2: Extraction of R0 and RC parameters ---');
[f_data, p_data] = uigetfile('*.txt', 'Select the test dataset (e.g., MPDch)');
if isequal(f_data,0); error('Missing dataset'); end
file_path = fullfile(p_data, f_data);

% A. R0 Extraction
[r0_vals, ~, ~] = get_r0(file_path, current_step);
R0_mean = mean(r0_vals); % Average for the constant model

% B. RC Parameters Extraction
% We use the simplified version without the 'min_current_fit' filter.
[soc_rc, R1_vec, C1_vec, R2_vec, C2_vec] = get_rc_parameters(file_path, Qn_nominal, initial_soc, current_step);

% Calculate averages (Algorithm 1, Step 9 of the article)
R1_m = mean(R1_vec);
C1_m = mean(C1_vec);
R2_m = mean(R2_vec);
C2_m = mean(C2_vec);

%% 4. Summary and Display of Results
fprintf('\n====================================================\n');
fprintf('FINAL PARAMETERS FOR YOUR MODEL\n');
fprintf('====================================================\n');
fprintf('R0 (Internal Resistance)  : %.4f Ohm\n', R0_mean);
fprintf('R1 (Fast Polarization)    : %.4f Ohm | C1 : %.1f F\n', R1_m, C1_m);
fprintf('R2 (Slow Polarization)    : %.4f Ohm | C2 : %.1f F\n', R2_m, C2_m);
fprintf('====================================================\n');

% Graph of RC parameters vs SoC
figure('Name', 'Evolution of Identified Parameters');
subplot(2,1,1); hold on; grid on;
plot(soc_rc, R1_vec, 'o-b', 'DisplayName', 'R1 (Electrochemical)');
plot(soc_rc, R2_vec, 's-r', 'DisplayName', 'R2 (Concentration)');
ylabel('Resistance [\Omega]'); legend; set(gca, 'XDir', 'reverse');

subplot(2,1,2); hold on; grid on;
plot(soc_rc, C1_vec, 'o-g', 'DisplayName', 'C1');
plot(soc_rc, C2_vec, 's-m', 'DisplayName', 'C2');
ylabel('Capacitance [F]'); xlabel('SoC [%]'); legend; set(gca, 'XDir', 'reverse');

sgtitle(['RC Parameters vs SoC - File: ', f_data]);


%% --- MODEL SIMULATION BLOCK ---
% Reload raw data for comparison
data = readmatrix(file_path, 'Delimiter', ';', 'NumHeaderLines', 1);
t_vec = data(:, 1); V_mes = data(:, 2); I_vec = data(:, 3);
soc_sim = get_soc(file_path, Qn_nominal, initial_soc);

% Variable preparation
N = length(I_vec);
x1 = zeros(N, 1); x2 = zeros(N, 1); V_sim = zeros(N, 1);
tau1 = R1_m * C1_m; tau2 = R2_m * C2_m;

% Simulation loop (Logic from your photos)
for k = 2:N
    dt = t_vec(k) - t_vec(k-1);
    
    % Update polarization voltages (x1 and x2)
    x1(k) = exp(-dt/tau1)*x1(k-1) + R1_m*(1 - exp(-dt/tau1))*I_vec(k-1);
    x2(k) = exp(-dt/tau2)*x2(k-1) + R2_m*(1 - exp(-dt/tau2))*I_vec(k-1);
    
    % Calculate simulated voltage: OCV - R0*I - x1 - x2
    V_sim(k) = polyval(p_coeffs_ocv, soc_sim(k)) - R0_mean * I_vec(k) - x1(k) - x2(k);
end

%% --- VALIDATION GRAPH (VOLTAGE COMPARISON) ---
figure('Color', 'w', 'Name', 'Validation: Measurement vs Simulation');

% Upper subplot: Both voltage curves
subplot(2,1,1); hold on; grid on;
plot(t_vec, V_mes, 'k', 'LineWidth', 1.5, 'DisplayName', 'Measured Voltage');
plot(t_vec, V_sim, '--r', 'LineWidth', 1.5, 'DisplayName', 'Simulated Voltage');
ylabel('Voltage [V]', 'FontWeight', 'bold');
title('Model Validation: Voltage Comparison');
legend('Location', 'best');

% Lower subplot: The error (V_mes - V_sim)
subplot(2,1,2); hold on; grid on;
error_vec = V_mes - V_sim;
plot(t_vec, error_vec, 'b', 'DisplayName', 'Error (V_{mes} - V_{sim})');
yline(0, '--k', 'HandleVisibility', 'off'); % Zero line
ylabel('Error [V]', 'FontWeight', 'bold');
xlabel('Time [s]', 'FontWeight', 'bold');
title(sprintf('Prediction Error (RMSE = %.4f V)', sqrt(mean(error_vec.^2, 'omitnan'))));