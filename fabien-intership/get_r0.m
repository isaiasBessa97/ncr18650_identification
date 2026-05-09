function [soc_values, r0_values, v_values, i_values] = get_r0(file_path, Qn, initial_soc, min_step, delay_steps)
% GET_R0 Calcule la résistance interne (R0) et le SoC correspondant lors des coupures de courant.
%
% Inputs:
%   file_path   : Chemin complet vers le fichier de test (string)
%   Qn          : Capacité nominale de la cellule en Ah
%   initial_soc : SoC de départ du test (%)
%   min_step    : Seuil de détection d'une impulsion (ex: 0.5 A)
%   delay_steps : Nombre de secondes pour laisser la tension chuter (ex: 2)
%
% Outputs:
%   soc_values  : Tableau contenant le SoC à l'instant de chaque impulsion
%   r0_values   : Tableau contenant les valeurs de R0 calculées (Ohms)
%   v_values    : Tableau contenant la tension avant la coupure
%   i_values    : Tableau contenant le courant avant la coupure

    %% 1. Calcul du SoC Continu
    % Appel de la fonction externe pour générer l'axe SoC
    soc_full = get_soc(file_path, Qn, initial_soc);

    %% 2. Lecture des Données
    data = readmatrix(file_path, 'Delimiter', ';', 'NumHeaderLines', 1);
    
    % Extraction Tension (Colonne 2) et Courant (Colonne 3)
    V = data(:, 2);
    I = data(:, 3); 
    
    %% 3. Détection des Impulsions (Coupures / Falls)
    % On cherche les moments où le courant passe de ON à OFF
    is_pulse_on = abs(I) > min_step;
    transitions = diff(is_pulse_on);
    
    % -1 signifie que le courant s'est coupé (fin de l'impulsion)
    fall_idx = find(transitions == -1); 
    
    %% 4. Calcul de R0 et extraction du SoC
    num_pulses = length(fall_idx);
    r0_values = zeros(num_pulses, 1);
    soc_values = zeros(num_pulses, 1);
    v_values = zeros(num_pulses, 1);
    i_values = zeros(num_pulses, 1);
    
    valid_count = 0;
    
    for j = 1:num_pulses
        idx = fall_idx(j); % L'instant exact avant la coupure
        
        % Sécurité : on vérifie qu'on ne déborde pas de la fin du fichier avec le délai
        if idx + delay_steps <= length(V)
            valid_count = valid_count + 1;
            
            % dI : Courant juste avant la coupure (I est à 0 juste après)
            dI = abs(I(idx)); 
            
            % dV : Chute de tension mesurée "delay_steps" secondes après la coupure
            dV = abs(V(idx + delay_steps) - V(idx));
            
            % Calcul de R0 = dV / dI
            r0_values(valid_count) = dV / dI;
            
            % Enregistrement des variables au moment exact de la coupure
            soc_values(valid_count) = soc_full(idx);
            v_values(valid_count) = V(idx);
            i_values(valid_count) = I(idx);
        end
    end
    
    % On nettoie les tableaux des zéros inutiles
    soc_values = soc_values(1:valid_count);
    r0_values = r0_values(1:valid_count);
    v_values = v_values(1:valid_count);
    i_values = i_values(1:valid_count);
end