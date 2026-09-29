close all; clear all; clc;

%% 1. Configuration et Choix du Test
disp('--- Analyse de la Résistance Interne (R0) vs SoC ---');
disp('Tests disponibles :');
disp('1 - MPD (ex: MPDch_000deg.csv)');
disp('2 - DST (ex: DST_000deg.csv)');
% Ajoutez simplement une ligne ici pour afficher votre nouveau test dans le menu
% disp('3 - Nouveau Test (ex: NOUVEAU_000deg.csv)');

test_choice = input('Entrez le numéro du test : ');

% Utilisation de switch pour ajouter facilement de nouveaux tests
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
        error('Choix invalide. Veuillez relancer le script.');
end

% Paramètres de la fonction
current_min_step = 0.5;
Qn = 3.05;         
initial_soc = 100; 
fall_delay_steps = 2; 

% Tableau des températures (vous pouvez y ajouter toutes les températures possibles)
temperatures = [0, 10, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65];

%% 2. Préparation du Graphique
figure('Color', 'w', 'Position', [100, 100, 1000, 600]); 
hold on; grid on;
xlabel('State of Charge (SoC) [%]', 'FontWeight', 'bold');
ylabel('Internal Resistance R_0 (\Omega)', 'FontWeight', 'bold');
title(sprintf('Internal Resistance vs SoC - Test : %s', test_prefix), 'FontSize', 14);

% Générer une palette de couleurs dynamiques
colors = lines(length(temperatures)); 

%% 3. Boucle de Traitement sur les Températures
for k = 1:length(temperatures)
    t = temperatures(k);
    
    % Construction du nom du fichier
    temp_str = sprintf('%03ddeg', t);
    file_name = sprintf('%s_%s.csv', test_prefix, temp_str);
    file_path = fullfile(base_dir, file_name);
    
    % Vérification de l'existence du fichier
    if ~isfile(file_path)
        fprintf('Fichier introuvable pour %d °C, passage au suivant.\n', t);
        continue;
    end
    
    fprintf('Traitement en cours pour %d °C...\n', t);
    
    %% --- ADAPTATION DES DONNÉES SANS MODIFIER VOS FONCTIONS ---
    % 1. Lecture robuste du CSV (, ou ;)
    try
        data = readmatrix(file_path, 'NumHeaderLines', 1);
    catch
        opts = detectImportOptions(file_path);
        data = readmatrix(file_path, opts);
    end
    
    % 2. Réorganisation des colonnes : [Time(1), Voltage(3), Current(2), Temp(4)]
    if size(data, 2) >= 4
        data_formatted = [data(:, 1), data(:, 3), data(:, 2), data(:, 4)];
    else
        data_formatted = [data(:, 1), data(:, 3), data(:, 2)]; % Sécurité si pas de Temp
    end
    
    % 3. Création du fichier temporaire .txt (avec le bon délimiteur ;)
    temp_file = fullfile(tempdir, 'temp_pulse_test.txt');
    fid = fopen(temp_file, 'w'); fprintf(fid, 'Time;Voltage;Current;Temp\n'); fclose(fid);
    writematrix(data_formatted, temp_file, 'Delimiter', ';', 'WriteMode', 'append');
    
    %% --- EXÉCUTION DU CALCUL ---
    [soc_array, r0_array, v_array, i_array] = get_r0_rise_and_fall(temp_file, Qn, initial_soc, current_min_step, fall_delay_steps);
    
    delete(temp_file);
    
    %% --- SÉPARATION ET TRACÉ ---
    is_rise = abs(i_array) < 0.1; 
    is_fall = ~is_rise;

    soc_rise = soc_array(is_rise);
    r0_rise  = r0_array(is_rise);

    soc_fall = soc_array(is_fall);
    r0_fall  = r0_array(is_fall);
    
    % Tracé Rise (Cercles pleins)
    if ~isempty(soc_rise)
        plot(soc_rise, r0_rise, 'o', 'LineWidth', 1, 'MarkerSize', 6, ...
            'Color', colors(k,:), 'MarkerFaceColor', colors(k,:), ...
            'DisplayName', sprintf('%d °C - Pulse ON', t));
    end
    
    % Tracé Fall (Carrés vides)
    if ~isempty(soc_fall)
        plot(soc_fall, r0_fall, 's', 'LineWidth', 1.5, 'MarkerSize', 6, ...
            'Color', colors(k,:), ... 
            'DisplayName', sprintf('%d °C - Pulse OFF', t));
    end
end

%% 4. Formatage Final
xlim([0 100]); % Force l'affichage strictement entre 0 et 100%
set(gca, 'XDir', 'reverse'); % Inversion de l'axe X pour lire de 100% à 0%
legend('Location', 'bestoutside', 'NumColumns', 1);
hold off;
disp('Analyse terminée !');