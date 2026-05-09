import numpy as np
import matplotlib.pyplot as plt

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

import numpy as np

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

def get_soc(file_test, Qn, initial_soc):
    # Remplacer par ta vraie fonction qui calcule le SoC de référence
    data = np.loadtxt(file_test, delimiter=';', skiprows=1)
    time = data[:, 0]
    I_meas = data[:, 2]
    # Simple intégration
    soc = np.zeros(len(time))
    soc[0] = initial_soc
    for i in range(1, len(time)):
        soc[i] = soc[i-1] - (I_meas[i] * 1.0 / (Qn * 3600)) * 100
    return soc

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

# Création du polynôme (degré 9 comme dans le MATLAB)
p_coeffs_ocv = np.polyfit(soc_valid, V_valid, 9)

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
P = np.eye(5)
theta = np.array([[0.1], [0.1], [0.01], [0.01], [0.01]]) # Vecteur colonne 5x1

soc_estimated = np.zeros(N)
soc_estimated[0] = 60 # SOC initial du modèle

theta_history = np.zeros((N, 5))
R0_hist = np.zeros(N); R1_hist = np.zeros(N); C1_hist = np.zeros(N)
R2_hist = np.zeros(N); C2_hist = np.zeros(N)
V_model = np.zeros(N)

y_past = np.zeros(2)
u_past = np.zeros(2)

# --- 3. Boucle Temps Réel (BMS) ---
print('Starting Real-Time RLS with Coulomb Counting...')

for k in range(N):
    u_k = I_meas[k]
    
    # B. BMS SoC Estimation
    if k > 0:
        soc_estimated[k] = soc_estimated[k-1] - (u_k * Ts / (Qn * 3600)) * 100
        soc_estimated[k] = max(0, min(100, soc_estimated[k]))
        
    # C. OCV Estimation
    ocv_k = np.polyval(p_coeffs_ocv, soc_estimated[k])
    
    # D. Pure dynamics
    y_k = ocv_k - V_meas[k]
    
    # --- START OF RLS ---
    if k > 1: # Équivalent à k > 2 en MATLAB
        phi_k = np.array([y_past[0], y_past[1], u_k, u_past[0], u_past[1]])
        
        y_pred = (phi_k.T @ theta).item()
        V_model[k] = ocv_k - y_pred
        
        if abs(u_k) > 0.05 or abs(u_k - u_past[0]) > 0.05:
            theta, P = rls_step(y_k, phi_k, theta, P, lmbda)
    else:
        V_model[k] = V_meas[k]
        
    # E. Save parameters
    theta_history[k, :] = theta.flatten()
    
    # F. Convert to 2RC
    r0, r1, c1, r2, c2 = theta_to_2rc(theta, Ts)
    R0_hist[k] = r0; R1_hist[k] = r1; C1_hist[k] = c1
    R2_hist[k] = r2; C2_hist[k] = c2
    
    # G. Time shift (k devient k-1)
    y_past = np.array([y_k, y_past[0]])
    u_past = np.array([u_k, u_past[0]])

print('Simulation Finished!')

# --- 4. Validation et Affichage ---
valid_idx = time > 10

erreur_abs_V = V_meas[valid_idx] - V_model[valid_idx]
erreur_moyenne_V = np.mean(erreur_abs_V)
rmse_V = np.sqrt(np.mean(erreur_abs_V**2))

print('\n====================================')
print('   RLS 2RC MODEL PERFORMANCES       ')
print('====================================')
print(f'Mean Absolute Error : {erreur_moyenne_V:.4f} V')
print(f'RMSE                : {rmse_V:.4f} V')
print('====================================\n')

# Style des graphiques (facultatif, pour rendre plus joli)
plt.style.use('default')

# Plot 1: Voltage
plt.figure('Real-Time RLS - Voltage Validation', figsize=(10, 5))
plt.plot(time, V_meas, 'k', linewidth=1.5, label='V measured (Experimental)')
plt.plot(time, V_model, 'r--', linewidth=1.5, label='V model (RLS Estimation)')
plt.ylabel('Voltage (V)', fontweight='bold')
plt.xlabel('Time (s)', fontweight='bold')
plt.title(f'Measured vs Model 2RC (RMSE: {rmse_V:.4f} V)')
plt.legend(loc='best')
plt.grid(True)

# Plot 2: SoC
plt.figure('SoC Validation', figsize=(10, 4))
plt.plot(time, soc_true, 'k', linewidth=1.5, label='True SoC')
plt.plot(time, soc_estimated, 'b--', linewidth=1.5, label='Estimated SoC (Coulomb Counting)')
plt.ylabel('State of Charge (%)', fontweight='bold')
plt.xlabel('Time (s)', fontweight='bold')
plt.title('Ground Truth SoC vs BMS Estimated SoC')
plt.legend(loc='best')
plt.grid(True)

# Plot 3: 2RC Parameters
valid_p = time > 50

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), num='2RC Parameters vs Estimated SoC')

# Top graph: Resistances
ax1.scatter(soc_estimated[valid_p], R0_hist[valid_p], s=15, c='k', label='R_0 (Internal)')
ax1.scatter(soc_estimated[valid_p], R1_hist[valid_p], s=15, c='b', label='R_1 (Fast)')
ax1.scatter(soc_estimated[valid_p], R2_hist[valid_p], s=15, c='r', label='R_2 (Slow)')
ax1.set_xlabel('Estimated State of Charge (%)', fontweight='bold')
ax1.set_ylabel('Resistance (Ohms)', fontweight='bold')
ax1.set_title('Real-Time Estimated Resistances')
ax1.legend(loc='best')
ax1.set_ylim([0, 0.5])
ax1.invert_xaxis() # set(gca, 'XDir', 'reverse')
ax1.grid(True)

# Bottom graph: Capacitances
ax2.scatter(soc_estimated[valid_p], C1_hist[valid_p], s=15, c='g', label='C_1 (Fast)')
ax2.scatter(soc_estimated[valid_p], C2_hist[valid_p], s=15, c='m', label='C_2 (Slow)')
ax2.set_xlabel('Estimated State of Charge (%)', fontweight='bold')
ax2.set_ylabel('Capacitance (F)', fontweight='bold')
ax2.set_title('Real-Time Estimated Capacitances')
ax2.legend(loc='best')
ax2.set_ylim([0, 80000])
ax2.invert_xaxis()
ax2.grid(True)

plt.tight_layout()
plt.show()