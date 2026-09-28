close all; clear all; clc;
%% Description
% Script modifié pour obtenir un pseudo OCV calculé basé sur la moyenne de la
% charge et décharge d'une cellule à différentes températures (0.05 C-RATE).

%% 1. Configuration
% Définir les dossiers contenant les fichiers de charge et de décharge 
% (À modifier avec les chemins réels sur votre ordinateur)
charge_dir = 'C:\ncr18650_identification\fabien-intership\simulated_temperature_data\CCCV_tests';
discharge_dir = 'C:\ncr18650_identification\fabien-intership\simulated_temperature_data\CD_tests';

% Définir le vecteur des températures à analyser (de 0 à 70 degrés par pas de 5)
temperatures = 0:5:70;

% Préparation du graphique
figure('Color', 'w'); hold on; grid on;
xlabel('State of Charge (SoC) [%]', 'FontWeight', 'bold');
ylabel('OCV (V)', 'FontWeight', 'bold');
title('Comparaison des courbes Pseudo-OCV par Température', 'FontSize', 14);

%% 2. Boucle de traitement automatisé
for t = temperatures
    
    % Création de l'identifiant de température (ex: t=5 -> '005deg')
    % %03d signifie "entier avec 3 chiffres, rempli de zéros"
    temp_str = sprintf('%03ddeg', t); 
    
    % Construction du nom exact des fichiers CSV basés sur vos images
    % Ex: CCCV_005C_005deg.csv et CDch_005C_005deg.csv
    file_charge_name = sprintf('CCCV_005C_%s.csv', temp_str);
    file_discharge_name = sprintf('CDch_005C_%s.csv', temp_str);
    
    % Création des chemins complets
    file_charge = fullfile(charge_dir, file_charge_name);
    file_discharge = fullfile(discharge_dir, file_discharge_name);
    
    % Vérification de sécurité : vérifier si les fichiers existent bien
    if ~isfile(file_charge) || ~isfile(file_discharge)
        fprintf('Fichiers pour %d degrés non trouvés. Passage à la température suivante.\n', t);
        continue; % Passe directement à la prochaine itération de la boucle
    end
    
    % Appel de la fonction pour calculer l'OCV
    [soc, ocv, qn] = get_ocv(file_charge, file_discharge);
    
    % Tracer la courbe avec une légende dynamique incluant la température et la capacité
    legend_str = sprintf('%d °C (Qn = %.3f Ah)', t, qn);
    plot(soc, ocv, 'LineWidth', 1.5, 'DisplayName', legend_str);
    
    % Message console pour suivre la progression
    fprintf('Température %d °C traitée avec succès !\n', t);
end

% Affichage de la légende finale
legend('Location', 'best');
hold off;