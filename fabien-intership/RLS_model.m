% clear; clc; close all;

%% 1. Configuration et Chargement
[f_ch, p_ch] = uigetfile('*.txt', 'Select the CHARGE file (0.05C)');
[f_dis, p_dis] = uigetfile('*.txt', 'Select the DISCHARGE file (0.05C)');
[f_test, p_test] = uigetfile('*.txt', 'Select the TEST file (MPD)');

file_charge = fullfile(p_ch, f_ch);
file_discharge = fullfile(p_dis, f_dis);
file_test = fullfile(p_test, f_test);

Qn = 3.05; 
initial_soc = 50; % SoC de départ

% On récupère la courbe OCV moyenne
[soc_ocv, V_ocv_raw, ~] = get_ocv(file_charge, file_discharge);

% Nettoyage et création du Polynôme OCV (degré 9)
valid_idx = ~isnan(V_ocv_raw);
soc_valid = soc_ocv(valid_idx);
V_valid = V_ocv_raw(valid_idx);
p_coeffs_ocv = polyfit(soc_valid, V_valid, 9);

% On récupère le SoC continu de l'essai
soc_array = get_soc(file_test, Qn, initial_soc);
soc_array_model = get_soc(file_test, Qn, 100);
%verifie si on est entre 0 et 100 et non entre 0 et 1
if max(soc_array) <= 1.05
    soc_array = soc_array * 100;
end

% Chargement des données du test
data = readmatrix(file_test, 'Delimiter', ';', 'NumHeaderLines', 1);
time = data(:, 1);
V_meas = data(:, 2);
I_meas = data(:, 3); % I > 0 en décharge

N = length(time);
Ts = 1; % Temps d'échantillonnage
 

%% 2. Initialisation du RLS
lambda = 0.9999; % Facteur d'oubli (Mémoire : proche de 1)
P = 1 * eye(5); % Matrice de confiance (1000 = très peu confiant au début)
theta = [0.1; 0.1; 0.01; 0.01; 0.01]; % Poids initiaux [a1, a2, b0, b1, b2]

% Tableaux pour stocker les résultats "temps réel"
theta_history = zeros(N, 5);
R0_hist = zeros(N, 1); R1_hist = zeros(N, 1); C1_hist = zeros(N, 1);
R2_hist = zeros(N, 1); C2_hist = zeros(N, 1);
V_model = zeros(N, 1);

% Variables d'état passées (y(k-1), y(k-2), u(k-1), u(k-2))
y_past = [0; 0];
u_past = [0; 0];

%% 3. Boucle Temps Réel (Le cœur de la simulation BMS)
disp('Lancement du RLS Temps Réel...');

for k = 1:N
    % A. Lecture des capteurs à l'instant k
    u_k = I_meas(k); 
    soc_k = soc_array(k);
    
    % B. Calcul de l'OCV actuel (Interpolation)
    ocv_k = polyval(p_coeffs_ocv, soc_k);
    
    % C. Calcul de la dynamique pure (Variable d'observation y_k)
    y_k = ocv_k - V_meas(k); 
    
    % --- DÉBUT DU RLS ---
    if k > 2
        % Construction du vecteur d'observation Phi
        phi_k = [y_past(1); y_past(2); u_k; u_past(1); u_past(2)];
        
        % Prédiction de la tension du modèle AVANT mise à jour
        y_pred = phi_k' * theta;
        V_model(k) = ocv_k - y_pred;
        
  
        % On met à jour SEULEMENT si on est dans une phase de courant dynamique
        % ou si on vient juste de couper le courant (relaxation)
        if abs(u_k) > 0.05 || abs(u_k - u_past(1)) > 0.05
            [theta, P] = rls_step(y_k, phi_k, theta, P, lambda);
        end
    else
        V_model(k) = V_meas(k); % Les 2 premiers pas, le modèle suit la mesure
    end
    % --- FIN DU RLS ---
    
    % D. Sauvegarde des paramètres mathématiques
    theta_history(k, :) = theta';
    
    % E. Conversion en grandeurs physiques 2RC
    [r0, r1, c1, r2, c2] = theta_to_2rc(theta, Ts);
    R0_hist(k) = r0; R1_hist(k) = r1; C1_hist(k) = c1;
    R2_hist(k) = r2; C2_hist(k) = c2;

    % F. Model validation

    
    % H. Glissement temporels des variables (k devient k-1)
    y_past = [y_k; y_past(1)];
    u_past = [u_k; u_past(1)];
end

disp('Terminé !');

%% 4. Affichage des Résultats
% --- Calcul des Erreurs ---
% On ignore les toutes premières secondes (le temps que le RLS converge) pour ne pas fausser la moyenne
valid_idx = time > 10; 

% Calcul de l'Erreur Moyenne Absolue (MAE) en millivolts
erreur_abs_V = V_meas(valid_idx) - V_model(valid_idx);
erreur_moyenne_V = mean(erreur_abs_V);

% Calcul de la Root Mean Square Error (RMSE) en millivolts
rmse_V = sqrt(mean((V_meas(valid_idx) - V_model(valid_idx)).^2));

%Calcul 
moy_R0 = mean(R0_hist);

% --- Affichage dans la Command Window (Console) ---
fprintf('\n====================================\n');
fprintf('   PERFORMANCES DU MODÈLE RLS 2RC   \n');
fprintf('====================================\n');
fprintf('Erreur Moyenne Absolue : %.2f V\n', erreur_moyenne_V);
fprintf('Erreur RMSE            : %.2f V\n', rmse_V);
fprintf('====================================\n\n');

% --- Création du Graphique Unique ---
figure('Name', 'RLS Temps Réel - Validation Tension', 'Position', [100, 100, 1000, 500]);

% On trace uniquement les tensions sur un seul grand graphique
plot(time, V_meas, 'k', 'LineWidth', 1.5); hold on;
plot(time, V_model, 'r--', 'LineWidth', 1.5);

ylabel('Tension (V)', 'FontWeight', 'bold'); 
xlabel('Temps (s)', 'FontWeight', 'bold');

% On intègre l'erreur moyenne directement dans le titre du graphique
title(sprintf('Measure vs Model 2RC (Error RMSE : %.2f V)', rmse_V), 'FontSize', 12);
legend('V measured (Expérimental)', 'V model (Model RLS)', 'Location', 'best'); 
grid on;

% Affichage de l'évolution de R0 en fonction du SoC (Inchangé)
figure('Name', 'Évolution de R0');
valid_r0_idx = time > 50 & R0_hist > 0 & R0_hist < 0.1; 
scatter(soc_array(valid_r0_idx), R0_hist(valid_r0_idx), 10, time(valid_r0_idx), 'filled');
title(sprintf('R0 function of SOC mean = %2f', moy_R0));
xlabel('State of Charge (%)'); ylabel('R0 (\Omega)');
set(gca, 'XDir', 'reverse'); 
grid on;
% --- Affichage de l'évolution des Paramètres (R0, R1, R2, C1, C2) ---
figure('Name', 'Évolution des Paramètres 2RC vs SoC', 'Position', [150, 150, 800, 800]);

% On enlève juste le gros pic de démarrage (les 50 premières secondes)
valid_p = time > 50;

% --- Graphique du haut : Les Résistances ---
subplot(2,1,1); hold on; grid on;
scatter(soc_array(valid_p), R0_hist(valid_p), 15, 'k', 'filled', 'DisplayName', 'R_0 (Interne)');
scatter(soc_array(valid_p), R1_hist(valid_p), 15, 'b', 'filled', 'DisplayName', 'R_1 (Rapide)');
scatter(soc_array(valid_p), R2_hist(valid_p), 15, 'r', 'filled', 'DisplayName', 'R_2 (Lente)');

xlabel('State of Charge (%)', 'FontWeight', 'bold'); 
ylabel('Résistance (\Omega)', 'FontWeight', 'bold');
title('Évolution des Résistances estimées en temps réel');
legend('Location', 'best');
ylim([0, 0.5]); % Limite fixe (zoome sur les valeurs intéressantes)
set(gca, 'XDir', 'reverse'); 

% --- Graphique du bas : Les Capacités ---
subplot(2,1,2); hold on; grid on;
scatter(soc_array(valid_p), C1_hist(valid_p), 15, 'g', 'filled', 'DisplayName', 'C_1 (Rapide)');
scatter(soc_array(valid_p), C2_hist(valid_p), 15, 'm', 'filled', 'DisplayName', 'C_2 (Lente)');

xlabel('State of Charge (%)', 'FontWeight', 'bold'); 
ylabel('Capacité (F)', 'FontWeight', 'bold');
title('Évolution des Capacités estimées en temps réel');
legend('Location', 'best');
ylim([0, 80000]); % Limite fixe (ajuste à 50000 ou 100000 si besoin selon ta batterie)
set(gca, 'XDir', 'reverse');