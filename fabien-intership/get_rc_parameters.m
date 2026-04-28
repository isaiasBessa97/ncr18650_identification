function [soc_results, R1_array, C1_array, R2_array, C2_array] = get_rc_parameters(file_path, Qn, initial_soc, current_min_step)
% GET_RC_PARAMETERS Extrait les paramètres R1, C1, R2, C2 d'un essai.
% (Filtre de courant minimum supprimé selon ta demande)

    % 1. Chargement des données et calcul du SoC
    soc_full = get_soc(file_path, Qn, initial_soc);
    
    data = readmatrix(file_path, 'Delimiter', ';', 'NumHeaderLines', 1);
    time = data(:, 1); V = data(:, 2); I = data(:, 3);
    
    % 2. Détection des impulsions (Méthode robuste par état)
    is_pulse_on = abs(I) > current_min_step;
    transitions = diff(is_pulse_on);
    rise_idx = find(transitions == 1);
    fall_idx = find(transitions == -1);
    
    if ~isempty(fall_idx) && ~isempty(rise_idx) && fall_idx(1) < rise_idx(1)
        fall_idx(1) = []; 
    end
    num_pulses = min(length(rise_idx), length(fall_idx));
    
    % Initialisation
    soc_results = zeros(num_pulses, 1);
    R1_array = zeros(num_pulses, 1); R2_array = zeros(num_pulses, 1);
    C1_array = zeros(num_pulses, 1); C2_array = zeros(num_pulses, 1);
    
    % 3. Configuration du Curve Fitting (Équation 42 de l'article)
    ft = fittype('y_e - x1_c*exp(-x/tau1) - x2_c*exp(-x/tau2)', ...
        'independent', 'x', 'problem', 'y_e', ...
        'coefficients', {'tau1', 'tau2', 'x1_c', 'x2_c'});
    
    opts = fitoptions('Method', 'NonlinearLeastSquares', 'Display', 'Off');
    opts.Lower = [0.1,  20,   0,   0]; opts.Upper = [20, 2000, 0.1, 0.1];
    opts.StartPoint = [5, 100, 0.05, 0.05]; 
    
    valid_count = 0;
    
    % 4. Traitement de chaque impulsion
    for k = 1:num_pulses
        t0_idx = rise_idx(k); tc_idx = fall_idx(k);
        
        if k < num_pulses
            te_idx = rise_idx(k+1) - 1;
        else
            te_idx = length(V);
        end
        
        if (te_idx - tc_idx) < 10; continue; end
        
        u = max(abs(I(t0_idx+1 : tc_idx)));
        t_on = time(tc_idx) - time(t0_idx); 
        t_relax = time(tc_idx+1 : te_idx) - time(tc_idx+1); 
        y_relax = V(tc_idx+1 : te_idx); 
        
        safe_points = min(10, length(y_relax));
        y_e = mean(y_relax(end-safe_points+1:end)); 
        y_relax(y_relax > y_e) = y_e;
        
        try
            [fitresult, ~] = fit(t_relax, y_relax, ft, opts, 'problem', y_e);
            
            % Calcul des paramètres R1, C1, R2, C2 (Équations 44 et 45)
            R1 = fitresult.x1_c / (u * (1 - exp(-t_on / fitresult.tau1)));
            R2 = fitresult.x2_c / (u * (1 - exp(-t_on / fitresult.tau2)));
            
            valid_count = valid_count + 1;
            soc_results(valid_count) = soc_full(tc_idx);
            R1_array(valid_count) = R1; C1_array(valid_count) = fitresult.tau1 / R1;
            R2_array(valid_count) = R2; C2_array(valid_count) = fitresult.tau2 / R2;
        catch
        end
    end
    
    soc_results = soc_results(1:valid_count); R1_array = R1_array(1:valid_count);
    C1_array = C1_array(1:valid_count); R2_array = R2_array(1:valid_count);
    C2_array = C2_array(1:valid_count);
end