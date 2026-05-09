import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# 1. FONCTIONS UTILITAIRES (DÉJÀ REMPLIES)
# =============================================================================

def rls_step(y_k, phi_k, theta_prev, P_prev, lmbda):
    
    phi_k = phi_k.reshape(-1, 1)
    den = lmbda + phi_k.T @ P_prev @ phi_k
    K = (P_prev @ phi_k) / den
    
    y_pred = (phi_k.T @ theta_prev).item()
    e = y_k - y_pred
    
    theta_new = theta_prev + K * e
    P_new = (1 / lmbda) * (P_prev - K @ phi_k.T @ P_prev)
    
    return theta_new, P_new

def theta_to_2rc(theta, Ts):
    
    a1, a2, b0, b1, b2 = theta.flatten()
    
    E, F, G = 1 + a1 - a2, 1 - a1 - a2, 1 + a2
    if abs(E) < 1e-4 or abs(F) < 1e-4:
        return 0.05, 0.01, 5000, 0.01, 40000 # Valeurs refuges si division par zéro

    R0 = (b0 - b1 + b2) / E
    
    term1 = Ts * G / F
    term2 = (Ts**2 * E) / (4 * F)
    delta = max(0, term1**2 - 4 * term2)
    
    tau1 = (term1 + np.sqrt(delta)) / 2
    tau2 = (term1 - np.sqrt(delta)) / 2
    
    R_tot = (b0 + b1 + b2) / F
    sum_tau = tau1 + tau2
    
    R1 = (R_tot - R0) * (tau1 / sum_tau) if sum_tau > 0 else 0
    R2 = (R_tot - R0) * (tau2 / sum_tau) if sum_tau > 0 else 0
    
    C1 = tau1 / R1 if R1 > 0 else 0
    C2 = tau2 / R2 if R2 > 0 else 0
    
    return R0, R1, C1, R2, C2
# =============================================================================
# 1. DÉFINITION DES PARAMÈTRES ET VARIABLES (INITIALISATION)
# =============================================================================

# --- Paramètres de la Batterie ---
Qn_Ah = 3.05                    # Capacité nominale en Ah
Qn_As = Qn_Ah * 3600            # Capacité en Ampères-secondes (Coulombs)
Ts = 1.0                        # Temps d'échantillonnage (secondes)

# Paramètres du modèle 2RC (Supposés connus/constants pour ce Kalman standard)
# (Si vous voulez les estimer en même temps, il faudra un EKF ou un Dual Kalman)
R0 = 0.05                       # Résistance interne (Ohm)
R1 = 0.01                       # Résistance branche 1 (Ohm)
C1 = 5000.0                     # Capacité branche 1 (Farad)
R2 = 0.02                       # Résistance branche 2 (Ohm)
C2 = 40000.0                    # Capacité branche 2 (Farad)

# Modèle OCV Linéaire (Vocv = alpha * SoC + beta)
# Ici, SoC est exprimé en % (de 0 à 100)
alpha = 0.008                   # Pente (V / %)
beta = 3.1                      # Ordonnée à l'origine (V à 0% SoC)

# --- Matrices d'État du système Continu -> Discret ---
# Vecteur d'état X = [SoC, V_c1, V_c2].T

# Matrice A : Dynamique interne du système (Transition d'état)
A = np.array([
    [1, 0, 0],
    [0, np.exp(-Ts / (R1 * C1)), 0],
    [0, 0, np.exp(-Ts / (R2 * C2))]
])

# Matrice B : Impact de l'entrée (Courant u_k) sur les états
# Attention au signe : on suppose que I > 0 décharge la batterie
B = np.array([
    [-(Ts * 100) / Qn_As],               # Variation du SoC en %
    [R1 * (1 - np.exp(-Ts / (R1 * C1)))], # Variation de Vc1
    [R2 * (1 - np.exp(-Ts / (R2 * C2)))]  # Variation de Vc2
]).reshape(-1, 1) # Assure que c'est un vecteur colonne (3x1)

# Matrice C : Comment les états sont liés à la mesure (Observation)
# V_meas = alpha*SoC + beta - Vc1 - Vc2 - R0*I
C = np.array([[alpha, -1, -1]]) # La constante 'beta' sera gérée à part

# Matrice D : Lien direct entre l'entrée (I) et la mesure (V)
D = np.array([[-R0]])

# --- Variables du Filtre de Kalman ---
# État initial X_k (On démarre à 100% de SoC, et condensateurs déchargés)
x_k = np.array([[100.0], [0.0], [0.0]]) 

# Matrice de covariance de l'erreur initiale P_k
# Représente notre confiance dans x_k (Petit = très confiant)
P_k = np.eye(3) * 0.1 

# Matrice de bruit de processus Q (Ce que vous avez noté Q_k sur votre cahier)
# Modélise l'incertitude du modèle mathématique
Q = np.diag([1e-5, 1e-4, 1e-4]) 

# Bruit de mesure R (scalaire car on ne mesure que la tension)
# Bruit des capteurs (ex: 10 mV d'incertitude au carré)
R = np.array([[0.01**2]])

# =============================================================================
# 2. SIMULATION DES DONNÉES (À remplacer par votre lecture de fichier)
# =============================================================================
N = 1000
time = np.arange(0, N * Ts, Ts)
I_meas = np.ones(N) * 1.5 # Courant de décharge constant de 1.5A pour l'exemple
V_meas = np.zeros(N)      # On va stocker ici la tension mesurée simulée (ou lue)

# On simule une tension mesurée fictive pour faire tourner le code
# (Dans votre vrai code, chargez V_meas depuis le fichier texte)
soc_vrai = 100.0
for k in range(N):
    soc_vrai -= (I_meas[k] * Ts * 100) / Qn_As
    V_meas[k] = (alpha * soc_vrai + beta) - (0.01 * I_meas[k]) + np.random.randn()*0.01

# =============================================================================
# 3. BOUCLE TEMPS RÉEL : LE FILTRE DE KALMAN
# =============================================================================

# Historique pour l'affichage
soc_est_hist = np.zeros(N)
V_est_hist = np.zeros(N)

print("Lancement du Filtre de Kalman...")

for k in range(N):
    u_k = I_meas[k]
    y_k = V_meas[k] # La vraie mesure de tension
    
    # -------------------------------------------------------------
    # ÉTAPE 1 : PRÉDICTION (Correspond à vos notes manuscrites)
    # -------------------------------------------------------------
    
    # Prédiction de l'état (X_k|k-1)
    x_pred = A @ x_k + B * u_k
    
    # Prédiction de la covariance (P_k|k-1)
    P_pred = A @ P_k @ A.T + Q
    
    # Prédiction de la mesure de tension (Y_pred)
    # y = Cx + Du + beta (on ajoute beta car notre OCV est affine, pas strictement linéaire)
    y_pred = C @ x_pred + D * u_k + beta 
    
    # -------------------------------------------------------------
    # ÉTAPE 2 : CORRECTION (Mise à jour avec la mesure)
    # -------------------------------------------------------------
    
    # Innovation (Erreur de prédiction de la mesure)
    innovation = y_k - y_pred
    
    # Calcul du Gain de Kalman (K)
    # K = P_pred * C.T * inv(C * P_pred * C.T + R)
    S = C @ P_pred @ C.T + R
    K = P_pred @ C.T @ np.linalg.inv(S)
    
    # Mise à jour de l'état estimé (X_k|k)
    x_k = x_pred + K @ innovation
    
    # Mise à jour de la matrice de covariance (P_k|k)
    # P = (I - K*C) * P_pred
    I_mat = np.eye(3)
    P_k = (I_mat - K @ C) @ P_pred
    
    # -------------------------------------------------------------
    # SAUVEGARDE
    # -------------------------------------------------------------
    soc_est_hist[k] = x_k[0, 0]
    
    # La tension estimée finale est calculée avec le nouvel état mis à jour
    V_est_hist[k] = (C @ x_k + D * u_k + beta).item()

print("Terminé !")

# =============================================================================
# 4. AFFICHAGE
# =============================================================================
plt.figure(figsize=(12, 6))

plt.subplot(2, 1, 1)
plt.plot(time, V_meas, label='Tension Mesurée', color='black')
plt.plot(time, V_est_hist, '--', label='Tension Estimée (Kalman)', color='red')
plt.ylabel('Tension (V)')
plt.legend()
plt.grid()
plt.title('Suivi de Tension par Filtre de Kalman')

plt.subplot(2, 1, 2)
plt.plot(time, soc_est_hist, label='SoC Estimé (%)', color='blue')
plt.ylabel('SoC (%)')
plt.xlabel('Temps (s)')
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()