close all; clear all; clc;
%% Description 
% This script find the coefficient of a polynomial curve modeling the
% pseudo OCV curve we calculated with the average of a charge and discharge
% with a C-Rate of 0.05 to be more accurate
%% 1. File Selection
disp('Select the CHARGE file for OCV...');
[file_ch_name, path_ch] = uigetfile({'*.txt;*.csv', 'Data Files'; '*.*', 'All Files'}, 'Select Charge File');
if isequal(file_ch_name,0); error('Cancelled by user'); end
file_charge = fullfile(path_ch, file_ch_name);

disp('Select the DISCHARGE file for OCV...');
[file_dis_name, path_dis] = uigetfile({'*.txt;*.csv', 'Data Files'; '*.*', 'All Files'}, 'Select Discharge File');
if isequal(file_dis_name,0); error('Cancelled by user'); end
file_discharge = fullfile(path_dis, file_dis_name);

%% 2. Calculate Pseudo-OCV
disp('Calculating Pseudo-OCV...');
% Make sure get_ocv.m is in the same folder!
[soc_axe, V_average, Qn] = get_ocv(file_charge, file_discharge);

%% 3. Polynomial Fitting (Least Squares)
fprintf('\n--- OCV Polynomial Modeling ---\n');
% Ask the user for the polynomial degree in the command window
poly_degree = input('Enter the desired polynomial degree (Examples 1,2,3,4,5,6,7,8,9): ');

disp('Running Least Squares fitting...');


% On trouve les indices où V_average n'est pas un NaN (i.o.w that we search
% the points where we have a voltage for a SOC 
valid_idx = ~isnan(V_average);

% On crée des tableaux "propres" pour le calcul mathématique
soc_valid = soc_axe(valid_idx);
V_valid = V_average(valid_idx);

% On utilise polyfit uniquement sur les données valides
p_coeffs = polyfit(soc_valid, V_valid, poly_degree);

% On évalue le modèle polynomial pour tout l'axe SoC
V_model = polyval(p_coeffs, soc_axe);

%% 4. Evaluate the Model (RMSE)
% On calcule les valeurs du modèle uniquement pour les points valides
V_model_valid = polyval(p_coeffs, soc_valid);

% 1. RMSE absolu (en Volts)
rmse = sqrt(mean((V_valid - V_model_valid).^2));

% 2. RMSE en Pourcentage (rapport entre l'erreur et la vraie tension)
rmse_pct = 100 * sqrt(mean(((V_valid - V_model_valid) ./ V_valid).^2));

% Affichage dans la console (on utilise %% pour afficher le symbole %)
fprintf('Model RMSE : %.4f V (Erreur relative : %.4f %%)\n', rmse, rmse_pct);

%% 5. Plotting the Comparison
figure('Color', 'w'); hold on; grid on;

% Plot the Original Pseudo-OCV (Solid thick blue line)
plot(soc_axe, V_average, 'b', 'LineWidth', 2.5, 'DisplayName', 'Actual Pseudo-OCV (Data)');

% Plot the Fitted Polynomial Model (Dashed yellow/red line)
plot(soc_axe, V_model, '--r', 'LineWidth', 2, 'DisplayName', sprintf('Polynomial Model (Degree %d)', poly_degree));

% Formatting
xlabel('State of Charge (SoC) [%]', 'FontWeight', 'bold');
ylabel('Open Circuit Voltage (OCV) [V]', 'FontWeight', 'bold');

% On ajoute le pourcentage au titre du graphique !
title_str = sprintf('Pseudo-OCV vs Model (Degree %d | RMSE = %.4f V [%.4f %%])', poly_degree, rmse, rmse_pct);
title(title_str, 'FontSize', 14);
legend('Location', 'best');

hold off;
%% 6. Display Coefficients for External Use
fprintf('\n====================================================\n');
fprintf('Polynomial Coefficients (Highest degree to constant):\n');
fprintf('====================================================\n');
disp(p_coeffs');