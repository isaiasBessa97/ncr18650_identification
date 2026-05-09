function [theta_new, P_new] = rls_step(y_k, phi_k, theta_prev, P_prev, lambda)
% RLS_STEP Exécute une itération de l'algorithme Recursive Least Squares
%
% Inputs:
%   y_k        : Valeur mesurée à l'instant k (scalaire)
%   phi_k      : Vecteur des observations [y(k-1); y(k-2); u(k); u(k-1); u(k-2)]
%   theta_prev : Vecteur des paramètres estimés à l'instant k-1 
%   P_prev     : Matrice de covariance à l'instant k-1 
%   lambda     : Facteur d'oubli 
%
% Outputs:
%   theta_new  : Nouveaux paramètres mis à jour
%   P_new      : Nouvelle matrice de covariance

    % 1. Calcul du Gain de Kalman (K)
    % k(n) = P(n-1)*x(n) / (lambda + x^T(n)*P(n-1)*x(n))
    K = (P_prev * phi_k) / (lambda + phi_k' * P_prev * phi_k);
    
    % 2. Calcul de l'erreur d'estimation a priori (e)
    % e(n) = d(n) - y_pred(n)
    y_pred = phi_k' * theta_prev;
    e = y_k - y_pred;
    
    % 3. Mise à jour du vecteur de paramètres (theta)
    theta_new = theta_prev + K * e;
    
    % 4. Mise à jour de la matrice de covariance (P)
    P_new = (1 / lambda) * (P_prev - K * phi_k' * P_prev);
    
end