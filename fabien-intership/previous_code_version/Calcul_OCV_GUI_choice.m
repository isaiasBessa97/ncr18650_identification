close all; clear all; clc;

%% Description
% Script developed to get the OCV curve based on CC-CV charge and CD
% discharge

%% 1. Parameters and Loading
% Prompt user for the CHARGE file
disp('Select charge file...');
[file_charge, path_charge] = uigetfile('*.*', 'Select the CHARGE CSV file');
if isequal(file_charge,0); error('Cancelled by user'); end

% Prompt user for the DISCHARGE file
disp('Select discharge file...');
[file_discharge, path_discharge] = uigetfile('*.*', 'Select the DISCHARGE CSV file');
if isequal(file_discharge,0); error('Cancelled by user'); end

% Read files using the full path
data_ch = readmatrix(fullfile(path_charge, file_charge), 'Delimiter', ';', 'NumHeaderLines', 1);
data_dis = readmatrix(fullfile(path_discharge, file_discharge), 'Delimiter', ';', 'NumHeaderLines', 1);

% Extract columns (Column 2 = Voltage, Column 3 = Current)
V_ch = data_ch(:, 2);  I_ch = abs(data_ch(:, 3));
V_dis = data_dis(:, 2); I_dis = abs(data_dis(:, 3));

% Calculate nominal capacity in Ah (sum of current because current sampled
% each second)
Qn = sum(I_dis) / 3600; 

%% 2. SoC Calculation via Iterative Method (Coulomb Counting)
% Since dt = 1s, it does not appear in the equation (multiplication by 1)

% --- CHARGE Processing ---
soc_ch = zeros(length(V_ch), 1);
soc_ch(1) = 0; % Assuming we start from 0%
for ii = 2:length(V_ch)
    % Adding current since it is a charge phase
    soc_ch(ii) = soc_ch(ii-1) + (1 / (3600 * Qn)) * I_ch(ii);
end

% --- DISCHARGE Processing ---
soc_dis = zeros(length(V_dis), 1);
soc_dis(1) = 1; % Assuming we start from 100% (1.0)
for ii = 2:length(V_dis)
    % Subtracting current for the discharge phase
    soc_dis(ii) = soc_dis(ii-1) - (1 / (3600 * Qn)) * I_dis(ii);
end

% Convert to percentage for display
soc_ch_pct = soc_ch * 100;
soc_dis_pct = soc_dis * 100;

%% 3. Average Curve Calculation (OCV)
% Problem: the two files do not have the same number of rows.
% To calculate the average (V_ch + V_dis)/2, they must be on the same SoC axis.

% Create a common SoC axis from 0 to 100% with 1000 points
soc_axe = linspace(0, 100, 1000)';

% Align the voltages on this new axis (Linear interpolation)
V_ch_aligned = interp1(soc_ch_pct, V_ch, soc_axe, 'linear', 'extrap');
V_dis_aligned = interp1(soc_dis_pct, V_dis, soc_axe, 'linear', 'extrap');

% Calculate the average
V_average = (V_ch_aligned + V_dis_aligned) / 2;

%% 4. Plotting
figure; hold on; grid on;
plot(soc_ch_pct, V_ch, 'b', 'DisplayName', 'Charge (0.05C)');
plot(soc_dis_pct, V_dis, 'r', 'DisplayName', 'Discharge (0.05C)');
plot(soc_axe, V_average, 'k', 'LineWidth', 2, 'DisplayName', 'Average (OCV)');

xlabel('State of Charge (%)');
ylabel('Voltage (V)');
title('Voltage vs SoC Curves');
legend('Location', 'best');