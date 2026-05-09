% close all; clc;

%% 1. Configuration des Paramètres de la Cellule
% Propriétés physiques de la cellule Panasonic
Qn_nominal = 3.05;      % Capacité nominale (Ah)
initial_soc = 100;      % SoC initial du test (%)
Ts = 1;                 % Temps d'échantillonnage de simulation (s)
current_step = 0.5;

%% 2. ÉTAPE 1 : Modélisation de l'OCV (Polynôme)
disp('--- ÉTAPE 1: Modélisation de l''OCV ---');
[f_ch, p_ch] = uigetfile('*.txt;*.csv', 'Sélectionnez le fichier CHARGE (0.05C)');
[f_dis, p_dis] = uigetfile('*.txt;*.csv', 'Sélectionnez le fichier DISCHARGE (0.05C)');

if isequal(f_ch,0) || isequal(f_dis,0); error('Fichiers OCV manquants'); end

% Appel de ta fonction get_ocv
[soc_ocv, V_ocv_raw, ~] = get_ocv(fullfile(p_ch, f_ch), fullfile(p_dis, f_dis));

% Calcul du polynôme 
poly_degree = 9; 
valid_idx = ~isnan(V_ocv_raw); 
soc_valid = soc_ocv(valid_idx);
V_valid = V_ocv_raw(valid_idx);
p_coeffs_ocv = polyfit(soc_valid, V_valid, poly_degree);

disp('Modèle OCV généré avec succès.');

%% 3. STEP 2: Extraction of R0 and RC parameters
disp('--- STEP 2: Extraction of R0 and RC parameters ---');
[f_data, p_data] = uigetfile('*.txt', 'Select the test dataset (e.g., MPDch)');
if isequal(f_data,0); error('Missing dataset'); end
file_path = fullfile(p_data, f_data);

% A. R0 Extraction
[soc_r0, r0_vals, ~, ~] = get_r0(file_path, Qn_nominal, initial_soc, current_step, 2);

% B. RC Parameters Extraction
[soc_rc, R1_vec, C1_vec, R2_vec, C2_vec] = get_rc_parameters(file_path, Qn_nominal, initial_soc, current_step);

% C. NETTOYAGE POUR LUT (Obligatoire pour interp1)
[soc_rc_lut, idx_rc] = unique(soc_rc);
R1_lut = R1_vec(idx_rc);
C1_lut = C1_vec(idx_rc);
R2_lut = R2_vec(idx_rc);
C2_lut = C2_vec(idx_rc);

[soc_r0_lut, idx_r0] = unique(soc_r0);
R0_lut = r0_vals(idx_r0);

%% 4. ÉTAPE 3 : Chargement du Profil de Validation

% Lecture des données réelles du cycleur
data_val = readmatrix(file_path, 'Delimiter', ';', 'NumHeaderLines', 1);
t_vec = data_val(:, 1);
V_mes = data_val(:, 2);
I_vec = data_val(:, 3); % Vecteur courant u(k)
N = length(I_vec);

%% 5. ÉTAPE 4 : Boucle de Simulation Dynamique (Modèle LTV)
disp('--- ÉTAPE 4: Simulation du modèle LUT dynamique en cours... ---');

x = zeros(3, N); 
V_sim = zeros(N, 1);
x(:, 1) = [0; 0; initial_soc / 100]; 

for k = 1:N
    soc_percent = x(3, k) * 100; 
    
    % INTERPOLATION (sur les LUT nettoyées)
    % Le 'max(..., 1e-4)' empêche physiquement d'avoir une résistance ou capacité à 0 ou négative
    R1_k =  max(1e-4, interp1(soc_rc_lut, R1_lut, soc_percent, 'linear', 'extrap'));
    R2_k =  max(1e-4, interp1(soc_rc_lut, R2_lut, soc_percent, 'linear', 'extrap'));
    R0_k = max(1e-4, interp1(soc_r0_lut, R0_lut, soc_percent, 'linear', 'extrap'));
    
    % Pour les condensateurs, on les empêche juste d'être trop petits pour éviter la division par zéro
    C1_k = min(10000,max(100,interp1(soc_rc_lut, C1_lut, soc_percent, 'linear', 'extrap')));
    C2_k = min(10000,max(100,interp1(soc_rc_lut, C2_lut, soc_percent, 'linear', 'extrap')));
    
    phi_k = polyval(p_coeffs_ocv, soc_percent);
    
    % SÉCURITÉ EULER : On empêche le terme (1 - Ts/RC) de devenir négatif ou nul
    % C'est ce qui faisait exploser ton graphique à 10^307
    A11 = max(0.01, 1 - Ts/(R1_k*C1_k));
    A22 = max(0.01, 1 - Ts/(R2_k*C2_k));
    
    A_k = [ A11, 0,   0 ;
            0,   A22, 0 ;
            0,   0,   1 ];

    B_k = [ Ts/C1_k ;
            Ts/C2_k ;
           -Ts/(3600*Qn_nominal) ]; 
           
    C_k = [-1, -1, 0];
    D_k = -R0_k;
    
    V_sim(k) = C_k * x(:, k) + D_k * I_vec(k) + phi_k;
    
    if k < N
        x(:, k+1) = A_k * x(:, k) + B_k * I_vec(k);
    end
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

%% --- VISUALISATION DES PARAMÈTRES EXTRAITS ---
figure('Color', 'w', 'Name', 'Paramètres RC et R0 vs SoC', 'Position', [100, 100, 900, 700]);

% Graph du haut : Les Résistances (R0, R1, R2)
subplot(2,1,1); hold on; grid on;
plot(soc_r0, R0_lut, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'DisplayName', 'R_0 (Interne)');
plot(soc_rc, R1_lut, '-ob', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'DisplayName', 'R_1 (Rapide)');
plot(soc_rc, R2_lut, '-or', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'DisplayName', 'R_2 (Lente)');
xlabel('SoC [%]', 'FontWeight', 'bold');
ylabel('Résistance [\Omega]', 'FontWeight', 'bold');
title('Évolution des Résistances identifiées en fonction du SoC');
legend('Location', 'best');
set(gca, 'XDir', 'reverse'); % Pour afficher de 100% à 0%

% Graph du bas : Les Capacités (C1, C2)
subplot(2,1,2); hold on; grid on;
plot(soc_rc, C1_lut, '-og', 'LineWidth', 1.5, 'MarkerFaceColor', 'g', 'DisplayName', 'C_1 (Rapide)');
plot(soc_rc, C2_lut, '-om', 'LineWidth', 1.5, 'MarkerFaceColor', 'm', 'DisplayName', 'C_2 (Lente)');
xlabel('SoC [%]', 'FontWeight', 'bold');
ylabel('Capacité [F]', 'FontWeight', 'bold');
title('Évolution des Capacités identifiées en fonction du SoC');
legend('Location', 'best');
set(gca, 'XDir', 'reverse');