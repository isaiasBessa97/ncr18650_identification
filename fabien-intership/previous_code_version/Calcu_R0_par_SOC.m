%close all; clear all; clc;

%% 1. Configuration
% Pulse detection threshold
current_min_step = 0.5;

% Battery specifications for SoC calculation
Qn = 3.05;         % Nominal capacity in Ah (Adjust this for your cell)
initial_soc = 100; % Starting SoC of the test (%)

% Prompt user for the dynamic pulse test file
disp('Select the pulse test file...');
[file_name, path_name] = uigetfile({'*.txt;*.csv', 'Text/Data Files'; '*.*', 'All Files'}, 'Select the Pulse Test file');
if isequal(file_name,0); error('Cancelled by user'); end
file_path = fullfile(path_name, file_name);

%% 2. Calculate Continuous SoC & Read Data
disp('Calculating Continuous SoC...');
% Call the external function (make sure get_soc.m is in the same folder)
soc_full = get_soc(file_path, Qn, initial_soc);

% IMPORTANT: We must read V and I here in the main script 
% so we can use them for the synchronization in Step 4.
data = readmatrix(file_path, 'Delimiter', ';', 'NumHeaderLines', 1);
V = data(:, 2);
I = data(:, 3);

%% 3. Retrieve R0 using your UNMODIFIED get_r0 function
disp('Extracting R0 from get_r0...');
[r0_array, v_array, i_array] = get_r0(file_path, current_min_step);

%% 4. Match SoC to R0 (Recreating the indices externally)
% To find out EXACTLY what SoC corresponds to what R0, we mirror the 
% detection logic of get_r0 here to grab the correct SoC points.
dI = diff(I);
dV = diff(V);
raw_indices = find(abs(dI) >= current_min_step);

% Apply the exact same filter as get_r0 (R0 < 0.5)
r0_raw = abs(dV(raw_indices) ./ dI(raw_indices));
valid_indices = raw_indices(r0_raw < 0.5);

% Now we map the SoC exactly to the valid pulses!
soc_at_pulses = soc_full(valid_indices);

% Safety check: Ensure arrays matched perfectly
if length(soc_at_pulses) ~= length(r0_array)
    error('Mismatch in array lengths. Data alignment failed.');
end

%% 5. Classify into Rise (a->b) and Fall (c->d)
% i_array contains the current AFTER the pulse step.
% If |I| > 0.1 A, the pulse is ON (Rise). If |I| < 0.1 A, it returned to 0 (Fall).
is_rise = abs(i_array) > 0.1; 

% --- RISE Data (Pulse ON) ---
soc_rise = soc_at_pulses(is_rise);
r0_rise  = r0_array(is_rise);

% --- FALL Data (Pulse OFF) ---
soc_fall = soc_at_pulses(~is_rise);
r0_fall  = r0_array(~is_rise);

%% 6. Export Results to Excel File
disp('Preparing data for export...');

% Ask the user where to save the new Excel file
[save_name, save_path] = uiputfile('*.xlsx', 'Save the Results as an Excel file');

if ~isequal(save_name, 0)
    full_save_path = fullfile(save_path, save_name);
    
    % Create MATLAB Tables for clean exporting
    % We group SoC and R0 into columns with clear headers
    Table_Rise = table(soc_rise, r0_rise, 'VariableNames', {'SoC_Percentage', 'Internal_Resistance_R0_Ohms'});
    Table_Fall = table(soc_fall, r0_fall, 'VariableNames', {'SoC_Percentage', 'Internal_Resistance_R0_Ohms'});
    
    % Write the tables to the Excel file in separate sheets (tabs)
    writetable(Table_Rise, full_save_path, 'Sheet', 'Pulse_ON_Rise');
    writetable(Table_Fall, full_save_path, 'Sheet', 'Pulse_OFF_Fall');
    
    fprintf('\nSUCCESS! Data saved to: %s\n', save_name);
    disp('You can now open the Excel file to see the two separate tabs.');
else
    disp('Export cancelled by the user.');

    %% 7. Plotting the Results
figure('Color', 'w'); hold on; grid on;

% Plot Rise (Points bleus pour la montée du courant)
plot(soc_rise, r0_rise, 'ob', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'b', 'DisplayName', 'Pulse ON (Rise: a \rightarrow b)');

% Plot Fall (Carrés rouges pour la descente du courant)
plot(soc_fall, r0_fall, 'sr', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'r', 'DisplayName', 'Pulse OFF (Fall: c \rightarrow d)');

% Mise en forme du graphique
xlabel('State of Charge (SoC) [%]', 'FontWeight', 'bold');
ylabel('Internal Resistance R_0 (\Omega)', 'FontWeight', 'bold');
title('Internal Resistance vs SoC', 'FontSize', 14);
legend('Location', 'best');

% Inverser l'axe X pour aller de 100% vers 0% (standard pour la décharge)
set(gca, 'XDir', 'reverse'); 

hold off;
end