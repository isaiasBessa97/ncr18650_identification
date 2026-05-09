function [R0, R1, C1, R2, C2] = theta_to_2rc(theta, Ts)
% THETA_TO_2RC Convertit les paramètres discrets ARX en paramètres physiques 2RC
% Utilise la transformation bilinéaire (Tustin).
% theta = [a1; a2; b0; b1; b2]
% Ts = Temps d'échantillonnage (secondes)

    a1 = theta(1); a2 = theta(2);
    b0 = theta(3); b1 = theta(4); b2 = theta(5);
    
    % Variables intermédiaires pour simplifier les équations de Tustin
    E = 1 + a1 - a2;
    F = 1 - a1 - a2;
    G = 1 + a2;
    
    % Protection contre la division par zéro en début de RLS
    if abs(E) < 1e-4 || abs(F) < 1e-4
        R0=0; R1=0; C1=0; R2=0; C2=0; return;
    end

    % Calcul de R0
    R0 = (b0 - b1 + b2) / E;
    
    % Calcul des constantes de temps (tau1, tau2)
    % Basé sur les racines du polynôme caractéristique
    term1 = Ts * G / F;
    term2 = (Ts^2 * E) / (4 * F);
    
    % Résolution de l'équation du second degré pour les Taus
    delta = term1^2 - 4 * term2;
    if delta < 0; delta = 0; end % Sécurité physique
    
    tau1 = (term1 + sqrt(delta)) / 2;
    tau2 = (term1 - sqrt(delta)) / 2;
    
    % Somme des résistances (R_total = R0 + R1 + R2)
    R_tot = (b0 + b1 + b2) / F;
    
    % R1 et R2 (Approximation classique issue de l'identification)
    % La répartition exacte demande la résolution d'un système non linéaire
    % On utilise ici la méthode simplifiée d'extraction :
    R1 = (R_tot - R0) * (tau1 / (tau1 + tau2)); 
    R2 = (R_tot - R0) * (tau2 / (tau1 + tau2));
    
    % Calcul des capacités
    if R1 > 0; C1 = tau1 / R1; else; C1 = 0; end
    if R2 > 0; C2 = tau2 / R2; else; C2 = 0; end
end