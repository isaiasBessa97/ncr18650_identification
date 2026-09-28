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

# --- RÉINTRODUCTION DES FONCTIONS RLS ---
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
        return 0.05, 0.01, 5000, 0.01, 40000 
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
# ----------------------------------------

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
coeffs_R1 = [-8.824978e-05,5.466936e-03, 5.7490e-02]
coeffs_C1 = [-1.4373, 1.5117e+02, 6.9571e+03]
coeffs_R2 = [9.132070e-07, -5.766334e-05, 1.0217e-02]
coeffs_C2 = [4.8260e-01, 8.0579e+01, 9.0909e+03]

file_charge = r"C:\ncr18650_identification\dataset-thermal\BID003\BID003_CCCV005.0_02022026.txt"
file_discharge = r"C:\ncr18650_identification\dataset-thermal\BID003\BID003_CDch005.0_02022026.txt"
file_test = r"C:\Users\PRH\Downloads\MPDch_045deg.csv"

soc_ocv, V_ocv_raw, Qn = get_ocv(file_charge, file_discharge)
p_coeffs_ocv = np.polyfit(soc_ocv[~np.isnan(V_ocv_raw)], V_ocv_raw[~np.isnan(V_ocv_raw)], 9)
soc_true = get_soc(file_test, Qn, initial_soc)
if np.max(soc_true) <= 1.05: soc_true *= 100

time, V_meas, I_meas, Ts_meas, _ = load_battery_data(file_test)
N = len(time)

# =============================================================================
# --- 3. INITIALISATION (RLS + UKF) ---
# =============================================================================
# Variables RLS
lmbda = 0.9999
P = 1 * np.eye(5)
theta = np.array([[0.1], [0.1], [0.01], [0.01], [0.01]]) 
y_past = np.zeros(2)
u_past = np.zeros(2)

soc_estimated = np.zeros(N)
soc_estimated[0] = 80 # Gardé à 80 pour observer la convergence rapide du RLS
V_model = np.zeros(N)

b3 = -Ts/(3600*Qn)*100
x_k = np.array([[0],[0],[soc_estimated[0]]])
A_mat = np.array([[0,0,0],[0,0,0],[0,0,0]])
B_mat = np.array([[0],[0],[b3]])
D_mat = np.array([[-0.1]]) 
P_KF = np.diag([1, 1, 0.5]) 
Q = np.diag([1e-4, 1e-4, 1e-3]) 
R_kf = np.array([[0.0001]]) 
kappa = 1
weight = np.array([[0.4],[0.1],[0.1],[0.1],[0.1],[0.1],[0.1]])  
ocv_k = np.polyval(p_coeffs_ocv, soc_estimated[0])

print('Démarrage de la Simulation Hybride (RLS bridé par le modèle thermique)...')

# =============================================================================
# --- 4. BOUCLE HYBRIDE ---
# =============================================================================
for k in range(N):
    u_k = I_meas[k]
    T_cell = Ts_meas[k]
    
    # ----------------------------------------------------
    # A. La Base Thermique (La Théorie Rigide)
    # ----------------------------------------------------
    r0_th = max(1e-4, A_R0 * np.exp(B_R0 * T_cell))
    r1_th = max(1e-4, np.polyval(coeffs_R1, T_cell))
    c1_th = max(1.0, np.polyval(coeffs_C1, T_cell))
    r2_th = max(1e-4, np.polyval(coeffs_R2, T_cell))
    c2_th = max(1.0, np.polyval(coeffs_C2, T_cell))

    # ----------------------------------------------------
    # B. L'Observateur RLS (L'Adaptation Rapide)
    # ----------------------------------------------------
    y_rls = ocv_k - V_meas[k]
    if k > 1: 
        phi_k = np.array([y_past[0], y_past[1], u_k, u_past[0], u_past[1]])
        if abs(u_k) > 0.05 or abs(u_k - u_past[0]) > 0.05:
            theta, P = rls_step(y_rls, phi_k, theta, P, lmbda)
            
    r0_rls, r1_rls, c1_rls, r2_rls, c2_rls = theta_to_2rc(theta, Ts)

    # ----------------------------------------------------
    # C. Le Bridage (Hybridation)
    # ----------------------------------------------------
    # On autorise le RLS à varier de ± X% autour du modèle thermique
    marge = 1
    
    r0 = np.clip(r0_rls, r0_th * (1 - marge), r0_th * (1 + marge))
    r1 = np.clip(r1_rls, r1_th * (1 - marge), r1_th * (1 + marge))
    c1 = np.clip(c1_rls, c1_th * (1 - marge), c1_th * (1 + marge))
    r2 = np.clip(r2_rls, r2_th * (1 - marge), r2_th * (1 + marge))
    c2 = np.clip(c2_rls, c2_th * (1 - marge), c2_th * (1 + marge))

    # ----------------------------------------------------
    # D. Suite du UKF Normal
    # ----------------------------------------------------
    a1, a2 = np.exp(-Ts / (r1 * c1)), np.exp(-Ts / (r2 * c2))
    A_mat = np.array([[a1,0,0],[0,a2,0],[0,0,1]])
    B_mat = np.array([[r1*(1-a1)], [r2*(1-a2)], [b3]])
    D_mat = np.array([[-r0]])

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
    
    P_zn, P_xz = R_kf.item(), np.zeros((3, 1))
    for i in range(7): 
        ecart_x = sigma_points_pred[:,i].reshape(3,1) - x_pred
        ecart_z = output_sigma_pred[0,i] - y_pred
        P_zn += weight[i, 0].item() * (ecart_z ** 2)
        P_xz += weight[i, 0].item() * ecart_x * ecart_z
    
    Kn = P_xz / P_zn
    x_k = x_pred + Kn * (V_meas[k] - y_pred)
    P_KF = P_pred - P_zn * (Kn @ Kn.T)
    
    ocv_k = np.polyval(p_coeffs_ocv, x_k[2,0])
    V_model[k] = y_pred
    soc_estimated[k] = max(0, min(100, x_k[2, 0]))
    
    y_past = np.array([y_rls, y_past[0]])
    u_past = np.array([u_k, u_past[0]])

print('Simulation Terminée !')

# =============================================================================
# --- 5. RÉSULTATS ET AFFICHAGE ---
# =============================================================================
valid_idx = time > 10 
rmse_V = np.sqrt(np.mean((V_meas[valid_idx] - V_model[valid_idx])**2))
rmse_soc = np.sqrt(np.mean((soc_true[valid_idx] - soc_estimated[valid_idx])**2))

print(f"\n--- Validation Hybride (RLS + Thermique) ---")
print(f"RMSE Tension : {rmse_V:.4f} V")
print(f"RMSE SoC     : {rmse_soc:.2f} %")

plt.figure(figsize=(10, 8))

plt.subplot(2, 1, 1)
plt.plot(time, V_meas, label='Tension Mesurée', color='black')
plt.plot(time, V_model, label='Tension Prédite (Hybride)', color='red', linestyle='--')
plt.title(f'Précision de la Tension (RMSE = {rmse_V:.4f} V)', fontweight='bold')
plt.ylabel('Tension (V)')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.7)

plt.subplot(2, 1, 2)
plt.plot(time, soc_true, label='SoC Réel', color='black')
plt.plot(time, soc_estimated, label='SoC Estimé (UKF Hybride)', color='blue', linestyle='--')
plt.title(f'Suivi de l\'État de Charge (RMSE = {rmse_soc:.2f} %)', fontweight='bold')
plt.xlabel('Temps (s)')
plt.ylabel('SoC (%)')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.7)

plt.tight_layout()
plt.show()