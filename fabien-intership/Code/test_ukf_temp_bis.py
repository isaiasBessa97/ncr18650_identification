import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import cholesky

# =============================================================================
# --- FONCTIONS DE BASE ---
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

def theta_to_gamma(theta_th, Ts):
    a1, a2, a3 = theta_th.flatten()
    gamma_1 = a2 / Ts
    gamma_2 = a3 / Ts
    return gamma_1, gamma_2

def get_ocv(file_charge, file_discharge):
    data_ch = np.loadtxt(file_charge, delimiter=';', skiprows=1)
    data_dis = np.loadtxt(file_discharge, delimiter=';', skiprows=1)

    V_ch = data_ch[:, 1]
    I_ch = np.abs(data_ch[:, 2])
    V_dis = data_dis[:, 1]
    I_dis = np.abs(data_dis[:, 2])

    Qn = np.sum(I_dis) / 3600.0

    soc_ch = np.zeros(len(V_ch))
    soc_ch[1:] = np.cumsum(I_ch[1:]) / (3600.0 * Qn)

    soc_dis = np.ones(len(V_dis)) 
    soc_dis[1:] = 1.0 - (np.cumsum(I_dis[1:]) / (3600.0 * Qn))

    soc_ch_pct = soc_ch * 100.0
    soc_dis_pct = soc_dis * 100.0

    soc_axe = np.linspace(0, 100, 1000)

    V_ch_aligned = np.interp(soc_axe, soc_ch_pct, V_ch)
    V_dis_aligned = np.interp(soc_axe, soc_dis_pct[::-1], V_dis[::-1])

    V_average = (V_ch_aligned + V_dis_aligned) / 2.0

    return soc_axe, V_average, Qn

def get_soc(file_path, Qn, initial_soc):
    data = np.loadtxt(file_path, delimiter=';', skiprows=1)
    time = data[:, 0]
    I = data[:, 2]
    
    num_points = len(I)
    soc_array = np.zeros(num_points)
    soc_array[0] = initial_soc
    
    for ii in range(1, num_points):
        dt = time[ii] - time[ii-1]
        if dt == 0:
            dt = 1  
        soc_array[ii] = soc_array[ii-1] - (100 * dt / (3600 * Qn)) * I[ii]
        
    return soc_array

def generate_sigma_points(x_k, P, kappa):
    N = x_k.shape[0] 
    nombre_de_points = 2 * N + 1
    sigma_points = np.zeros((N, nombre_de_points))
    sigma_points[:, 0] = x_k[:, 0]
    L = cholesky((N + kappa) * P, lower=True)
    
    for i in range(N):
        sigma_points[:, i + 1] = x_k[:, 0] + L[:, i]
        sigma_points[:, i + 1 + N] = x_k[:, 0] - L[:, i]
        
    return sigma_points

# =============================================================================
# --- 1. CONFIGURATION ET CHARGEMENT ---
# =============================================================================
Ts = 1.0  
Qn = 3.05 
initial_soc = 100

file_charge = r"C:\ncr18650_identification\dataset-thermal\BID003\BID003_CCCV005.0_02022026.txt"
file_discharge = r"C:\ncr18650_identification\dataset-thermal\BID003\BID003_CDch005.0_02022026.txt"
file_test = r"C:\ncr18650_identification\dataset-thermal\BID003\BID003_RSDch_24022026.txt"

soc_ocv, V_ocv_raw, _ = get_ocv(file_charge, file_discharge)
valid_idx = ~np.isnan(V_ocv_raw)
soc_valid = soc_ocv[valid_idx]
V_valid = V_ocv_raw[valid_idx]
p_coeffs_ocv = np.polyfit(soc_valid, V_valid, 9)

soc_true = get_soc(file_test, Qn, initial_soc)
if np.max(soc_true) <= 1.05:
    soc_true = soc_true * 100

data = np.loadtxt(file_test, delimiter=';', skiprows=1)
time = data[:, 0]
V_meas = data[:, 1]
I_meas = data[:, 2] 
N = len(time)

Ta_meas = data[:, 3]
Ts_plus = data[:, 4]
Ts_minus = data[:, 5]
Ts_meas = (Ts_plus + Ts_minus) / 2.0

# =============================================================================
# --- 2. INITIALISATION DES FILTRES ---
# =============================================================================
# A. RLS Électrique
lmbda = 0.9999
P = 1*np.eye(5)
theta = np.array([[0.1], [0.1], [0.01], [0.01], [0.01]]) 
theta_history = np.zeros((N, 5))
y_past = np.zeros(2)
u_past = np.zeros(2)

# B. RLS Thermique
lmbda_th = 0.999
P_th = 1 * np.eye(3)
theta_th = np.array([[0.99], [0.01], [0.001]]) 
T_s_est = np.zeros(N)
T_s_est[0] = Ts_meas[0]
T_s_past = Ts_meas[0]
T_a_past = Ta_meas[0]
H_past = 0.0

# C. Filtre UKF
soc_estimated = np.zeros(N)
soc_estimated[0] = 80 
V_model = np.zeros(N)
b3 = -Ts/(3600*Qn)*100
x_k = np.array([[0],[0],[soc_estimated[0]]])
P_KF = np.diag([1, 1, 0.5]) 
Q = np.diag([1e-4, 1e-4, 1e-3]) 
R_kf = np.array([[0.0001]]) 
kappa = 1
weight = np.array([[0.4],[0.1],[0.1],[0.1],[0.1],[0.1],[0.1]])  
ocv_k = np.polyval(p_coeffs_ocv, soc_estimated[0])

# D. RLS en Cascade pour le Couplage Thermo-Électrique (_bis)
lmbda_bis = 0.999
theta_R0_fit = np.array([[-3.0], [0.0]]) # Loi d'Arrhenius
P_R0_fit = 10 * np.eye(2)
theta_R1_fit = np.array([[0.0], [0.0], [0.01]])
P_R1_fit = 10 * np.eye(3)
theta_C1_fit = np.array([[0.0], [0.0], [1000.0]])
P_C1_fit = 10 * np.eye(3)
theta_R2_fit = np.array([[0.0], [0.0], [0.01]])
P_R2_fit = 10 * np.eye(3)
theta_C2_fit = np.array([[0.0], [0.0], [5000.0]])
P_C2_fit = 10 * np.eye(3)

# Historiques
R0_hist, R1_hist, C1_hist, R2_hist, C2_hist = (np.zeros(N) for _ in range(5))
gamma1_hist, gamma2_hist = np.zeros(N), np.zeros(N)
R0_bis, R1_bis, C1_bis, R2_bis, C2_bis = (np.zeros(N) for _ in range(5))

# =============================================================================
# --- 3. BOUCLE TEMPS RÉEL ---
# =============================================================================
print('Démarrage de la Simulation Temps Réel...')

for k in range(N):
    u_k = I_meas[k]
    
    # --- A. RLS ÉLECTRIQUE ---
    y_rls = ocv_k - V_meas[k]
    if k > 1: 
        phi_k = np.array([y_past[0], y_past[1], u_k, u_past[0], u_past[1]])
        if abs(u_k) > 0.05 or abs(u_k - u_past[0]) > 0.05:
            theta, P = rls_step(y_rls, phi_k, theta, P, lmbda)
 
    theta_history[k, :] = theta.flatten()
    r0, r1, c1, r2, c2 = theta_to_2rc(theta, Ts)
    r0 = max(r0, 1e-4); r1 = max(r1, 1e-4); c1 = max(c1, 1.0)   
    r2 = max(r2, 1e-4); c2 = max(c2, 1.0)   
    R0_hist[k], R1_hist[k], C1_hist[k], R2_hist[k], C2_hist[k] = r0, r1, c1, r2, c2

    a1 = np.exp(-Ts / (r1 * c1))
    a2 = np.exp(-Ts / (r2 * c2))
    b1 = r1 * (1 - a1)
    b2 = r2 * (1 - a2)

    A_mat = np.array([[a1,0,0],[0,a2,0],[0,0,1]])
    B_mat = np.array([[b1],[b2],[b3]])
    D_mat = np.array([[-r0]])

    # --- B. RLS THERMIQUE ---
    H_k = abs(u_k) * abs(ocv_k - V_meas[k])
    if k > 1:
        phi_th = np.array([T_s_past, -T_a_past, H_past])
        theta_th, P_th = rls_step(Ts_meas[k], phi_th, theta_th, P_th, lmbda_th)
    
    phi_th_current = np.array([T_s_past, -T_a_past, H_past])
    T_s_est[k] = (phi_th_current.T @ theta_th).item()

    g1, g2 = theta_to_gamma(theta_th, Ts)
    gamma1_hist[k], gamma2_hist[k] = g1, g2
    T_s_past, T_a_past, H_past = Ts_meas[k], Ta_meas[k], H_k

    # --- C. FILTRE UKF ---
    sigma_points = generate_sigma_points(x_k, P_KF, kappa)
    sigma_points_pred = np.zeros((3,7))
    for i in range(7):
        x_sigma_pred = A_mat @ sigma_points[:,i].reshape(3,1) + B_mat * u_k
        sigma_points_pred[:,i] =  x_sigma_pred.flatten()  

    x_pred = sigma_points_pred @ weight
    P_pred = np.copy(Q)
    for i in range(7):  
        ecart = sigma_points_pred[:,i].reshape(3,1) - x_pred
        P_pred += weight[i,0] * (ecart @ ecart.T) 

    output_sigma_points_pred = np.zeros((1,7))
    for i in range(7):
        ocv_point = np.polyval(p_coeffs_ocv, sigma_points_pred[2, i])
        y_pred_ukf = -sigma_points_pred[0, i] - sigma_points_pred[1, i] + (D_mat * u_k).item() + ocv_point
        output_sigma_points_pred[:,i] = y_pred_ukf

    y_pred = (output_sigma_points_pred @ weight).item()
    P_zn = R_kf.item() 
    P_xz = np.zeros((3, 1))

    for i in range(7): 
        ecart_x = sigma_points_pred[:,i].reshape(3,1) - x_pred
        ecart_z = output_sigma_points_pred[0,i] - y_pred
        P_zn += weight[i, 0].item() * (ecart_z ** 2)
        P_xz += weight[i, 0].item() * ecart_x * ecart_z
    
    Kn = P_xz / P_zn
    innovation = V_meas[k] - y_pred
    x_k = x_pred + Kn * innovation
    P_KF = P_pred - P_zn * (Kn @ Kn.T)
    
    ocv_k = np.polyval(p_coeffs_ocv, x_k[2,0])
    V_model[k] = y_pred
    soc_estimated[k] = max(0, min(100, x_k[2, 0]))
    
    y_past = np.array([y_rls, y_past[0]])
    u_past = np.array([u_k, u_past[0]])

    # --- D. RLS EN CASCADE (IDENTIFICATION THERMO-ÉLECTRIQUE _BIS) ---
    # On laisse passer 50 s pour éviter d'apprendre des données transitoires fausses
    if k > 50:
        # Arrhenius pour R0
        y_R0 = np.log(r0)
        phi_R0 = np.array([1.0, Ts_meas[k]])
        theta_R0_fit, P_R0_fit = rls_step(y_R0, phi_R0, theta_R0_fit, P_R0_fit, lmbda_bis)
        
        # Polynômes d'ordre 2 pour le reste
        phi_poly = np.array([Ts_meas[k]**2, Ts_meas[k], 1.0])
        theta_R1_fit, P_R1_fit = rls_step(r1, phi_poly, theta_R1_fit, P_R1_fit, lmbda_bis)
        theta_C1_fit, P_C1_fit = rls_step(c1, phi_poly, theta_C1_fit, P_C1_fit, lmbda_bis)
        theta_R2_fit, P_R2_fit = rls_step(r2, phi_poly, theta_R2_fit, P_R2_fit, lmbda_bis)
        theta_C2_fit, P_C2_fit = rls_step(c2, phi_poly, theta_C2_fit, P_C2_fit, lmbda_bis)

    # Calcul des _bis avec les équations du moment
    A_rt = np.exp(theta_R0_fit[0, 0])
    B_rt = theta_R0_fit[1, 0]
    R0_bis[k] = A_rt * np.exp(B_rt * Ts_meas[k])
    
    phi_poly_current = np.array([Ts_meas[k]**2, Ts_meas[k], 1.0])
    R1_bis[k] = (phi_poly_current @ theta_R1_fit).item()
    C1_bis[k] = (phi_poly_current @ theta_C1_fit).item()
    R2_bis[k] = (phi_poly_current @ theta_R2_fit).item()
    C2_bis[k] = (phi_poly_current @ theta_C2_fit).item()

print('Simulation Terminée !')

# =============================================================================
# --- 4. AFFICHAGES ---
# =============================================================================
valid_idx = time > 50 
rmse_V = np.sqrt(np.mean((V_meas[valid_idx] - V_model[valid_idx])**2))
rmse_soc = np.sqrt(np.mean((soc_true[valid_idx] - soc_estimated[valid_idx])**2))
rmse_temp = np.sqrt(np.mean((Ts_meas[valid_idx] - T_s_est[valid_idx])**2))

# --- FIGURE 1 : VALIDATION GLOBALE ---
plt.figure(figsize=(12, 10))
plt.subplot(3, 1, 1)
plt.plot(time, V_meas, label='Tension Mesurée', color='black')
plt.plot(time, V_model, label='Tension Modèle (UKF)', color='red', linestyle='--')
plt.title(f'Comparaison Tension (RMSE = {rmse_V:.4f} V)', fontweight='bold')
plt.ylabel('Tension (V)')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.7)

plt.subplot(3, 1, 2)
plt.plot(time, soc_true, label='SoC Réel', color='black')
plt.plot(time, soc_estimated, label='SoC Estimé (UKF)', color='blue', linestyle='--')
plt.title(f'État de Charge (RMSE = {rmse_soc:.2f} %)', fontweight='bold')
plt.ylabel('SoC (%)')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.7)

plt.subplot(3, 1, 3)
plt.plot(time, Ts_meas, label='Température Mesurée', color='black')
plt.plot(time, T_s_est, label='Température Estimée (RLS)', color='orange', linestyle='--')
plt.title(f'Température (RMSE = {rmse_temp:.2f} °C)', fontweight='bold')
plt.xlabel('Temps (s)')
plt.ylabel('Température (°C)')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.7)
plt.tight_layout()

# --- FIGURE 2 : RLS ÉLECTRIQUE VS MODÈLES _BIS (Cascade) ---
fig, axs = plt.subplots(3, 2, figsize=(14, 10))
fig.suptitle("Évolution des Paramètres 2RC et Modèles Cascade Temps Réel (_bis)", fontweight='bold', fontsize=14)

axs[0, 0].scatter(Ts_meas[valid_idx], R0_hist[valid_idx], s=2, alpha=0.5, color='blue', label='Composant RLS Brut')
axs[0, 0].scatter(Ts_meas[valid_idx], R0_bis[valid_idx], s=2, alpha=0.5, color='red', label='Modèle _bis (Arrhenius)')
axs[0, 0].set_title('R0 = f(Ts)')
axs[0, 0].legend()
axs[0, 0].grid(True)

axs[1, 0].scatter(Ts_meas[valid_idx], R1_hist[valid_idx], s=2, alpha=0.5, color='orange', label='Brut')
axs[1, 0].scatter(Ts_meas[valid_idx], R1_bis[valid_idx], s=2, alpha=0.5, color='black', label='_bis')
axs[1, 0].set_title('R1 = f(Ts)')
axs[1, 0].grid(True)

axs[2, 0].scatter(Ts_meas[valid_idx], C1_hist[valid_idx], s=2, alpha=0.5, color='green', label='Brut')
axs[2, 0].scatter(Ts_meas[valid_idx], C1_bis[valid_idx], s=2, alpha=0.5, color='black', label='_bis')
axs[2, 0].set_title('C1 = f(Ts)')
axs[2, 0].grid(True)

axs[0, 1].scatter(Ts_meas[valid_idx], R2_hist[valid_idx], s=2, alpha=0.5, color='red', label='Brut')
axs[0, 1].scatter(Ts_meas[valid_idx], R2_bis[valid_idx], s=2, alpha=0.5, color='black', label='_bis')
axs[0, 1].set_title('R2 = f(Ts)')
axs[0, 1].grid(True)

axs[1, 1].scatter(Ts_meas[valid_idx], C2_hist[valid_idx], s=2, alpha=0.5, color='purple', label='Brut')
axs[1, 1].scatter(Ts_meas[valid_idx], C2_bis[valid_idx], s=2, alpha=0.5, color='black', label='_bis')
axs[1, 1].set_title('C2 = f(Ts)')
axs[1, 1].grid(True)

axs[2, 1].axis('off')
plt.tight_layout()

# --- FIGURE 3 : ÉVOLUTION DES PARAMÈTRES GAMMA THERMIQUES ---
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(time[valid_idx], gamma1_hist[valid_idx], color='firebrick')
plt.title('Identification : Gamma 1 (Refroidissement)')
plt.xlabel('Temps (s)')
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot(time[valid_idx], gamma2_hist[valid_idx], color='teal')
plt.title('Identification : Gamma 2 (Échauffement)')
plt.xlabel('Temps (s)')
plt.grid(True)
plt.tight_layout()

plt.show()