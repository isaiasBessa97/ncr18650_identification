import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import cholesky

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


def get_ocv(file_charge, file_discharge):
    """
    Calcule la courbe OCV moyenne à partir des fichiers de charge et décharge.
    
    Arguments:
      file_charge    : Chemin complet vers le fichier de charge (str)
      file_discharge : Chemin complet vers le fichier de décharge (str)
      
    Retours:
      soc_axe        : Axe SoC commun (0 à 100%)
      V_average      : Vecteur de tension Pseudo-OCV résultant
      Qn             : Capacité nominale calculée (Ah)
    """
    

    data_ch = np.loadtxt(file_charge, delimiter=';', skiprows=1)
    data_dis = np.loadtxt(file_discharge, delimiter=';', skiprows=1)

    # Extraction des colonnes (Attention: Python commence à l'index 0 !)
    # MATLAB (:, 2) -> Python [:, 1] pour la tension
    # MATLAB (:, 3) -> Python [:, 2] pour le courant
    V_ch = data_ch[:, 1]
    I_ch = np.abs(data_ch[:, 2])
    
    V_dis = data_dis[:, 1]
    I_dis = np.abs(data_dis[:, 2])

    # Calcul de la capacité nominale (Qn) en Ah basée sur la décharge
    Qn = np.sum(I_dis) / 3600.0

    # On utilise np.cumsum() qui fait la somme cumulée instantanément.
    
    # --- CHARGE ---
    soc_ch = np.zeros(len(V_ch))
    soc_ch[1:] = np.cumsum(I_ch[1:]) / (3600.0 * Qn)

    # --- DÉCHARGE ---
    soc_dis = np.ones(len(V_dis)) # On initialise tout à 1
    soc_dis[1:] = 1.0 - (np.cumsum(I_dis[1:]) / (3600.0 * Qn))

    soc_ch_pct = soc_ch * 100.0
    soc_dis_pct = soc_dis * 100.0

    soc_axe = np.linspace(0, 100, 1000)

    # Interpolation linéaire avec np.interp
    V_ch_aligned = np.interp(soc_axe, soc_ch_pct, V_ch)
    
    # ATTENTION : np.interp exige que l'axe X (soc_dis_pct) soit strictement croissant.
    # Comme le SoC de décharge descend de 100 à 0, on doit inverser les tableaux 
    # avec la syntaxe [::-1] avant d'interpoler.
    V_dis_aligned = np.interp(soc_axe, soc_dis_pct[::-1], V_dis[::-1])

    # Moyenne
    V_average = (V_ch_aligned + V_dis_aligned) / 2.0

    return soc_axe, V_average, Qn


def get_soc(file_path, Qn, initial_soc):
    """
    Calcule le State of Charge (SoC) continu pour un essai donné.
    
    Arguments:
      file_path   : Chemin complet vers le fichier de test (str)
      Qn          : Capacité nominale de la cellule en Ah (float)
      initial_soc : SoC de départ de l'essai (ex: 100 pour une cellule chargée) (float)
      
    Retourne:
      soc_array   : Tableau numpy contenant la valeur du SoC pour chaque pas de temps
    """
    
    # Lecture des données
    # delimiter=';' et skiprows=1 remplacent les options readmatrix de MATLAB
    data = np.loadtxt(file_path, delimiter=';', skiprows=1)
    
    # En Python, l'indexation commence à 0. 
    # La colonne 1 de MATLAB devient 0, la colonne 3 devient 2.
    time = data[:, 0]
    I = data[:, 2]
    
    # Initialisation du tableau SoC
    num_points = len(I)
    soc_array = np.zeros(num_points)
    soc_array[0] = initial_soc
    
    # Boucle de comptage coulombien (Coulomb Counting)
    # L'indexation commence à 1 (qui correspond au 2ème élément)
    for ii in range(1, num_points):
        # Calcul du dt (pas de temps en secondes)
        dt = time[ii] - time[ii-1]
        
        if dt == 0:
            dt = 1  # Sécurité (Fallback)
            
        # Calcul du SoC en pourcentage
        soc_array[ii] = soc_array[ii-1] - (100 * dt / (3600 * Qn)) * I[ii]
        
    return soc_array



def generate_sigma_points(x_k, P, kappa):
    """
    Génère les 2N + 1 points sigma sous forme de vecteurs colonnes.
    
    Arguments:
    x_k   -- Vecteur d'état actuel (DOIT être de forme (N, 1))
    P     -- Matrice de covariance actuelle (forme (N, N))
    kappa -- Paramètre d'étalement
    """
    # 1. Déterminer la dimension N à partir des lignes du vecteur colonne
    N = x_k.shape[0] 
    nombre_de_points = 2 * N + 1
    
    # 2. Créer la matrice pour stocker les points (N lignes, 2N+1 colonnes)
    sigma_points = np.zeros((N, nombre_de_points))
    
    # 3. Placer x_k à l'index 0 (On extrait les valeurs de la colonne)
    sigma_points[:, 0] = x_k[:, 0]
    
    # 4. Calcul de la racine carrée de la matrice (Cholesky)
    L = cholesky((N + kappa) * P, lower=True)
    
    # 5. La Boucle : Création des voisins
    for i in range(N):
        # L[:, i] est une ligne/colonne de dispersion. 
        # On l'ajoute directement à la colonne de sigma_points
        
        # Voisins "positifs" (Index 1 à N)
        sigma_points[:, i + 1] = x_k[:, 0] + L[:, i]
        
        # Voisins "négatifs" (Index N+1 à 2N)
        sigma_points[:, i + 1 + N] = x_k[:, 0] - L[:, i]
        
    return sigma_points
# --- 1. Configuration et Chargement ---
Ts = 1.0  #Time sample period
Qn = 3.05 
initial_soc = 100

# Attention aux chemins Windows : utiliser r"..." pour "raw string"
file_charge = r"C:\ncr18650_identification\dataset-thermal\BID003\BID003_CCCV005.0_02022026.txt"
file_discharge = r"C:\ncr18650_identification\dataset-thermal\BID003\BID003_CDch005.0_02022026.txt"
file_test = r"C:\ncr18650_identification\dataset-thermal\BID003\BID003_MPDch_24022026.txt"

soc_ocv, V_ocv_raw, _ = get_ocv(file_charge, file_discharge)

valid_idx = ~np.isnan(V_ocv_raw)
soc_valid = soc_ocv[valid_idx]
V_valid = V_ocv_raw[valid_idx]

# Création du polynôme (degré 1 pour que ce soit lineaire)
p_coeffs_ocv = np.polyfit(soc_valid, V_valid, 9)

#Calcul de la derivee des coefficients 
dp_coeffs_ocv = np.polyder(p_coeffs_ocv)

# Chargement du SoC réel
soc_true = get_soc(file_test, Qn, initial_soc)
if np.max(soc_true) <= 1.05:
    soc_true = soc_true * 100

# Chargement des données de test
data = np.loadtxt(file_test, delimiter=';', skiprows=1)
time = data[:, 0]
V_meas = data[:, 1]
I_meas = data[:, 2] # I > 0 en décharge
N = len(time)

# --- 2. Initialisation RLS ---
lmbda = 0.9999
P = 1*np.eye(5)
theta = np.array([[0.1], [0.1], [0.01], [0.01], [0.01]]) # Vecteur colonne 5x1

soc_estimated = np.zeros(N)
soc_estimated[0] = 70 # SOC initial du modèle

theta_history = np.zeros((N, 5))
# R0_hist = np.zeros(N); R1_hist = np.zeros(N); C1_hist = np.zeros(N)
# R2_hist = np.zeros(N); C2_hist = np.zeros(N)
V_model = np.zeros(N)

y_past = np.zeros(2)
u_past = np.zeros(2)
# --- Initialisation Kalman --- 
b3 = -Ts/(3600*Qn)*100
#Vecteur d'etat
x_k = np.array([[0],[0],[soc_estimated[0]]])
A = np.array([[0,0,0],[0,0,0],[0,0,0]])
B = np.array([[0],[0],[b3]])
#C = np.array([[-1,-1,p_coeffs_ocv[0]]])
D = np.array([[-0.1]]) # -R0
#P_KF = np.diag([9, 9, 12]) # Plus on augmente plus le gain sera fort au depart
P_KF = np.diag([1, 1, 0.5]) # Plus on augmente plus le gain sera fort au depart 
P_zn = np.array([[0.]])
P_xz = np.array([[0.],[0.],[0.]])
wk = 0 # Process noise
Q = np.diag([1e-4, 1e-4, 1e-3]) # Process noise covariance 
R_kf = np.array([[0.0001]]) # Covariance of the measurement noise
kappa = 1
weight = np.array([[0.4],[0.1],[0.1],[0.1],[0.1],[0.1],[0.1]])  #Reminder: the index 0 is the measured value and sum of coefficients must = 1
soc_estimated[0] = x_k[2, 0]
soc_estimated[0] = max(0, min(100, soc_estimated[0]))
Kn = np.array([[0],[0],[0]])
        
    # C. OCV Estimation
ocv_k = np.polyval(p_coeffs_ocv, soc_estimated[0])
# --- 3. Boucle Temps Réel (BMS) ---
print('Starting Real-Time RLS with Coulomb Counting...')

for k in range(N):

   

    u_k = I_meas[k]
    
    # D. Pure dynamics
    y_rls = ocv_k - V_meas[k]
    
    # --- START OF RLS ---
    if k > 1: # Équivalent à k > 2 en MATLAB
        phi_k = np.array([y_past[0], y_past[1], u_k, u_past[0], u_past[1]])
        
        if abs(u_k) > 0.05 or abs(u_k - u_past[0]) > 0.05:
            theta, P = rls_step(y_rls, phi_k, theta, P, lmbda)
 
    # E. Save parameters
    theta_history[k, :] = theta.flatten()
    
    # F. Convert to 2RC
    r0, r1, c1, r2, c2 = theta_to_2rc(theta, Ts)
    # R0_hist[k] = r0; R1_hist[k] = r1; C1_hist[k] = c1
    # R2_hist[k] = r2; C2_hist[k] = c2

    # Securite division par zero
    r0 = max(r0, 1e-4)
    r1 = max(r1, 1e-4)
    c1 = max(c1, 1.0)   
    r2 = max(r2, 1e-4)
    c2 = max(c2, 1.0)   

    a1 = np.exp(-Ts / (r1 * c1))
    a2 = np.exp(-Ts / (r2 * c2))
    b1 = r1 * (1 - a1)
    b2 = r2 * (1 - a2)
    
    #a1 = 1-(Ts/(r1*c1))
    #a2 = 1-(Ts/(r2*c2))
    #b1 = Ts/c1
    #b2 = Ts/c2

    #---- UPDATE MATRIX ----

    A = np.array([[a1,0,0],[0,a2,0],[0,0,1]])
    B = np.array([[b1],[b2],[b3]])
    D = np.array([[-r0]])
    # C never change
   # --- Utilisation EKF (Méthode Matricielle) ---

    # 1. PRÉDICTION 
    sigma_points_pred = np.zeros((3,7))

    #ATTENTION ICI RISQUE DE PROBLEMES DE COMPATIBILITE DES TABLEAUX 
    sigma_points = generate_sigma_points(x_k,P_KF, kappa)
    
    for i in range(7):
          
          x_sigma_pred = A @ sigma_points[:,i].reshape(3,1) + B * u_k
          sigma_points_pred[:,i] =  x_sigma_pred.flatten()  
#RAPPEL: sigma_points_pred[:,i] renvoie une matrice (1,3) donc c'est pour cela qu'on utilise flatten qui vient mettre notre (3,1) en (1,3)

    
    x_pred = sigma_points_pred @ weight

    #Calcul of covariance matrix
    P_pred = np.copy(Q)


    #Matrice de covarian Pn 
    for i in range(7):  # ICI 7 CAR ON UTILISE UN VECTEUR D'ETAT DE 3 DIMENSIONS
        ecart = sigma_points_pred[:,i].reshape(3,1) - x_pred
        P_pred += weight[i,0] * (ecart @ ecart.T) 



    output_sigma_points_pred = np.zeros((1,7))
    
    for i in range(7):
        ocv_point = np.polyval(p_coeffs_ocv, sigma_points_pred[2, i])

        y_pred_ukf = -sigma_points_pred[0, i] - sigma_points_pred[1, i] + (D * u_k).item() + ocv_point
        output_sigma_points_pred[:,i] = y_pred_ukf

    y_pred = (output_sigma_points_pred @ weight).item()
    #Matrice de covarian Pz 
    #Matrice de covarian P_xz 
    P_zn = R_kf.item() 
    P_xz = np.zeros((3, 1))

    for i in range(7):  # ICI 7 CAR ON UTILISE UN VECTEUR D'ETAT DE 3 DIMENSIONS
        ecart_x = sigma_points_pred[:,i].reshape(3,1) - x_pred
        
        ecart_z = output_sigma_points_pred[0,i] - y_pred
        
        P_zn += weight[i, 0].item() * (ecart_z ** 2)
        P_xz += weight[i, 0].item() * ecart_x * ecart_z
    

    # Gain de Kalman (on garde bien la variable C)
    Kn = P_xz / P_zn

    innovation = V_meas[k] - y_pred
    # Mise à jour des états et de la covariance
    x_k = x_pred + Kn * innovation
    # Traduction littérale : P_nn = P_pred - K @ S @ K.T
    P_KF = P_pred - P_zn * (Kn @ Kn.T)
    
    ocv_k = np.polyval(p_coeffs_ocv, x_k[2,0])
    # Sauvegarde de la tension
    V_model[k] = y_pred

    soc_estimated[k] = x_k[2, 0]
    soc_estimated[k] = max(0, min(100, soc_estimated[k]))
    # G. Time shift (k devient k-1)
    y_past = np.array([y_rls, y_past[0]])
    u_past = np.array([u_k, u_past[0]])

print('Simulation Finished!')


# =============================================================================
# 4. CALCUL DES ERREURS (RMSE) ET AFFICHAGE DES RÉSULTATS
# =============================================================================

# --- Calcul du RMSE ---
# On ignore les 10 premières secondes (le temps que le Kalman converge)
valid_idx = time > 10 

# Formule du RMSE : Racine carrée de la moyenne des erreurs au carré
rmse_V = np.sqrt(np.mean((V_meas[valid_idx] - V_model[valid_idx])**2))
rmse_soc = np.sqrt(np.mean((soc_true[valid_idx] - soc_estimated[valid_idx])**2))

# Affichage dans la console
print(f"\nPerformances du Filtre de Kalman (après 10s) :")
print(f" -> RMSE Tension : {rmse_V:.4f} V")
print(f" -> RMSE SoC     : {rmse_soc:.2f} %")

# --- Graphique 1 : Comparaison des Tensions ---
plt.figure(figsize=(12, 5))
plt.plot(time, V_meas, label='Tension Measured (Expérimentale)', color='black', linewidth=1.5)
plt.plot(time, V_model, label='Tension modele(Kalman)', color='red', linestyle='--')
# On intègre le RMSE directement dans le titre
plt.title(f'Voltage Comparaison: Measured vs Model Kalman (RMSE = {rmse_V:.4f} V)', fontweight='bold')
plt.xlabel('Temps (s)', fontweight='bold')
plt.ylabel('Tension (V)', fontweight='bold')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.7)
plt.tight_layout()

# --- Graphique 2 : Comparaison du SoC ---
plt.figure(figsize=(12, 5))
plt.plot(time, soc_true, label='SoC REAL (Intégration Théorique)', color='black', linewidth=1.5)
plt.plot(time, soc_estimated, label='SoC Model (Filtre Kalman)', color='blue', linestyle='--')
# On intègre le RMSE directement dans le titre
plt.title(f'State of charge(SoC) : Real vs Estimated (RMSE = {rmse_soc:.2f} %)', fontweight='bold')
plt.xlabel('Temps (s)', fontweight='bold')
plt.ylabel('State of Charge (%)', fontweight='bold')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.7)
plt.tight_layout()

# Afficher les graphiques
plt.show()