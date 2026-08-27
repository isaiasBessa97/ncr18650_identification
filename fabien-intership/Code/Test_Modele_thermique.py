import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.linalg import cholesky

# =============================================================================
# --- 1. FONCTIONS DE BASE ---
# =============================================================================
def load_battery_data(file_path):
    df = pd.read_csv(file_path, sep=None, engine='python')
    if 'Tempo_s' in df.columns:
        time = df['Tempo_s'].values
        I_meas = df['Corrente_A'].values
        V_meas = df['Tensao_V'].values
        Ts_meas = df['Temperatura_C'].values
        Ta_meas = np.full(len(time), Ts_meas[0])
    else:
        data = df.values
        time = data[:, 0]
        V_meas = data[:, 1]
        I_meas = data[:, 2]
        Ta_meas = data[:, 3]
        Ts_meas = (data[:, 4] + data[:, 5]) / 2.0
    return time, V_meas, I_meas, Ts_meas, Ta_meas

def get_ocv(file_charge, file_discharge):
    _, V_ch, I_ch, _, _ = load_battery_data(file_charge)
    _, V_dis, I_dis, _, _ = load_battery_data(file_discharge)
    
    I_ch, I_dis = np.abs(I_ch), np.abs(I_dis)
    Qn = np.sum(I_dis) / 3600.0

    soc_ch = np.zeros(len(V_ch))
    soc_ch[1:] = np.cumsum(I_ch[1:]) / (3600.0 * Qn)
    soc_dis = np.ones(len(V_dis)) 
    soc_dis[1:] = 1.0 - (np.cumsum(I_dis[1:]) / (3600.0 * Qn))

    soc_axe = np.linspace(0, 100, 1000)
    V_ch_aligned = np.interp(soc_axe, soc_ch * 100.0, V_ch)
    V_dis_aligned = np.interp(soc_axe, (soc_dis * 100.0)[::-1], V_dis[::-1])

    return soc_axe, (V_ch_aligned + V_dis_aligned) / 2.0, Qn

def get_soc(file_path, Qn, initial_soc):
    time, _, I_meas, _, _ = load_battery_data(file_path)
    soc_array = np.zeros(len(I_meas))
    soc_array[0] = initial_soc
    for ii in range(1, len(I_meas)):
        dt = time[ii] - time[ii-1] if time[ii] - time[ii-1] != 0 else 1  
        soc_array[ii] = soc_array[ii-1] - (100 * dt / (3600 * Qn)) * I_meas[ii]
    return soc_array

def generate_sigma_points(x_k, P, kappa):
    N = x_k.shape[0] 
    sigma_points = np.zeros((N, 2 * N + 1))
    sigma_points[:, 0] = x_k[:, 0]
    L = cholesky((N + kappa) * P, lower=True)
    for i in range(N):
        sigma_points[:, i + 1] = x_k[:, 0] + L[:, i]
        sigma_points[:, i + 1 + N] = x_k[:, 0] - L[:, i]
    return sigma_points

# =============================================================================
# --- 2. CONFIGURATION ET ÉQUATIONS DU JUMEAU NUMÉRIQUE ---
# =============================================================================
Ts = 1.0  
initial_soc = 100


# Loi d'Arrhenius : R0 = A * exp(B * Ts)
A_R0 = 6.9191e-02      
B_R0 = -0.0193      

# Polynômes d'ordre 2 : a*Ts^2 + b*Ts + c
coeffs_R1 = [-8.824978e-05,5.466936e-03, 5.7490e-02]     # [a, b, c]
coeffs_C1 = [-1.4373, 1.5117e+02, 6.9571e+03]   # [a, b, c]
coeffs_R2 = [9.132070e-07, -5.766334e-05, 1.0217e-02]    # [a, b, c]
coeffs_C2 = [4.8260e-01, 8.0579e+01, 9.0909e+03]   # [a, b, c]
#Coeffice=ients trouves grace au fichier Model_recup_temp_UKF.py
# Fichiers OCV
file_charge = r"C:\ncr18650_identification\dataset-thermal\BID003\BID003_CCCV005.0_02022026.txt"
file_discharge = r"C:\ncr18650_identification\dataset-thermal\BID003\BID003_CDch005.0_02022026.txt"

# Le nouveau fichier de test à valider (ex: profil dynamique dynamique)
file_test = r"C:\Users\PRH\Downloads\MPDch_045deg.csv"

soc_ocv, V_ocv_raw, Qn = get_ocv(file_charge, file_discharge)
p_coeffs_ocv = np.polyfit(soc_ocv[~np.isnan(V_ocv_raw)], V_ocv_raw[~np.isnan(V_ocv_raw)], 9)
soc_true = get_soc(file_test, Qn, initial_soc)
if np.max(soc_true) <= 1.05: soc_true *= 100

time, V_meas, I_meas, Ts_meas, _ = load_battery_data(file_test)
N = len(time)

# =============================================================================
# --- 3. INITIALISATION DU FILTRE UKF ---
# =============================================================================

soc_estimated = np.zeros(N)
soc_estimated[0] = 100 
theta_history = np.zeros((N, 5))
V_model = np.zeros(N)

y_past = np.zeros(2)
u_past = np.zeros(2)

b3 = -Ts/(3600*Qn)*100
x_k = np.array([[0],[0],[soc_estimated[0]]])
A_mat = np.array([[0,0,0],[0,0,0],[0,0,0]])
B_mat = np.array([[0],[0],[b3]])
D_mat = np.array([[-0.1]]) 
P_KF = np.diag([1, 1, 0.5]) 
P_zn = np.array([[0.]])
P_xz = np.array([[0.],[0.],[0.]])
Q = np.diag([1e-4, 1e-4, 1e-3]) 
R_kf = np.array([[0.0001]]) 
kappa = 1
weight = np.array([[0.4],[0.1],[0.1],[0.1],[0.1],[0.1],[0.1]])  
soc_estimated[0] = max(0, min(100, x_k[2, 0]))
Kn = np.array([[0],[0],[0]])

print('Démarrage de la Simulation de Validation (Sans RLS)...')

# =============================================================================
# --- 4. BOUCLE DE VALIDATION RAPIDE ---
# =============================================================================
for k in range(N):
    u_k = I_meas[k]
    T_cell = Ts_meas[k]
    
    # 1. Calcul instantané des composants grâce aux équations thermiques
    r0 = max(1e-4, A_R0 * np.exp(B_R0 * T_cell))
    r1 = max(1e-4, np.polyval(coeffs_R1, T_cell))
    c1 = max(1.0, np.polyval(coeffs_C1, T_cell))
    r2 = max(1e-4, np.polyval(coeffs_R2, T_cell))
    c2 = max(1.0, np.polyval(coeffs_C2, T_cell))

    # 2. Mise à jour des matrices d'état
    a1, a2 = np.exp(-Ts / (r1 * c1)), np.exp(-Ts / (r2 * c2))
    A_mat = np.array([[a1,0,0],[0,a2,0],[0,0,1]])
    B_mat = np.array([[r1*(1-a1)], [r2*(1-a2)], [b3]])
    D_mat = np.array([[-r0]])

    # 3. UKF Prediction
    sigma_points = generate_sigma_points(x_k, P_KF, kappa)
    sigma_points_pred = np.zeros((3,7))
    for i in range(7):
        sigma_points_pred[:,i] = (A_mat @ sigma_points[:,i].reshape(3,1) + B_mat * u_k).flatten()  

    x_pred = sigma_points_pred @ weight
    P_pred = np.copy(Q)
    for i in range(7):  
        ecart = sigma_points_pred[:,i].reshape(3,1) - x_pred
        P_pred += weight[i,0] * (ecart @ ecart.T) 

    output_sigma_pred = np.zeros((1,7))
    for i in range(7):
        y_pred_ukf = -sigma_points_pred[0, i] - sigma_points_pred[1, i] + (D_mat * u_k).item() + np.polyval(p_coeffs_ocv, sigma_points_pred[2, i])
        output_sigma_pred[:,i] = y_pred_ukf

    y_pred = (output_sigma_pred @ weight).item()
    
    # 4. UKF Correction
    P_zn, P_xz = R_kf.item(), np.zeros((3, 1))
    for i in range(7): 
        ecart_x = sigma_points_pred[:,i].reshape(3,1) - x_pred
        ecart_z = output_sigma_pred[0,i] - y_pred
        P_zn += weight[i, 0].item() * (ecart_z ** 2)
        P_xz += weight[i, 0].item() * ecart_x * ecart_z
    
    Kn = P_xz / P_zn
    x_k = x_pred + Kn * (V_meas[k] - y_pred)
    P_KF = P_pred - P_zn * (Kn @ Kn.T)
    
    V_model[k] = y_pred
    soc_estimated[k] = max(0, min(100, x_k[2, 0]))

print('Simulation Terminée !')

# =============================================================================
# --- 5. RÉSULTATS ET AFFICHAGE ---
# =============================================================================
valid_idx = time > 10 
rmse_V = np.sqrt(np.mean((V_meas[valid_idx] - V_model[valid_idx])**2))
rmse_soc = np.sqrt(np.mean((soc_true[valid_idx] - soc_estimated[valid_idx])**2))

print(f"\n--- Validation du Modèle Thermo-Électrique ---")
print(f"RMSE Tension : {rmse_V:.4f} V")
print(f"RMSE SoC     : {rmse_soc:.2f} %")

plt.figure(figsize=(10, 8))

plt.subplot(2, 1, 1)
plt.plot(time, V_meas, label='Tension Mesurée', color='black')
plt.plot(time, V_model, label='Tension Prédite (Modèle Fixe)', color='red', linestyle='--')
plt.title(f'Précision de la Tension (RMSE = {rmse_V:.4f} V)', fontweight='bold')
plt.ylabel('Tension (V)')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.7)

plt.subplot(2, 1, 2)
plt.plot(time, soc_true, label='SoC Réel (Ampli Coulométrique)', color='black')
plt.plot(time, soc_estimated, label='SoC Estimé (UKF)', color='blue', linestyle='--')
plt.title(f'Suivi de l\'État de Charge (RMSE = {rmse_soc:.2f} %)', fontweight='bold')
plt.xlabel('Temps (s)')
plt.ylabel('SoC (%)')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.7)

plt.tight_layout()
plt.show()