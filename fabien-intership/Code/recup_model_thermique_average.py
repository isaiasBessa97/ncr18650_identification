import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.linalg import cholesky
import os

# =============================================================================
# --- 1. FONCTIONS OUTILS ET CHARGEMENT ---
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

def get_ocv(file_charge, file_discharge):
    _, V_ch, I_ch, _, _ = load_battery_data(file_charge)
    _, V_dis, I_dis, _, _ = load_battery_data(file_discharge)
    I_ch = np.abs(I_ch)
    I_dis = np.abs(I_dis)
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
    num_points = len(I_meas)
    soc_array = np.zeros(num_points)
    soc_array[0] = initial_soc
    for ii in range(1, num_points):
        dt = time[ii] - time[ii-1]
        if dt == 0:
            dt = 1  
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
# --- 2. CONFIGURATION GLOBALE ---
# =============================================================================
Ts = 1.0  
initial_soc = 100

file_charge = r"C:\Users\PRH\Downloads\CCCV_005C_025deg.csv"
file_discharge = r"C:\Users\PRH\Downloads\CDch_005C_025deg.csv"

test_files = [
    r"C:\Users\PRH\Downloads\DST_065deg.csv",
    r"C:\Users\PRH\Downloads\DST_055deg.csv",
    r"C:\Users\PRH\Downloads\DST_045deg.csv",
    r"C:\Users\PRH\Downloads\DST_035deg.csv",
    r"C:\Users\PRH\Downloads\DST_025deg.csv",
    r"C:\Users\PRH\Downloads\DST_060deg.csv",
    r"C:\Users\PRH\Downloads\DST_050deg.csv",
    r"C:\Users\PRH\Downloads\DST_040deg.csv",
    r"C:\Users\PRH\Downloads\DST_030deg.csv",
    r"C:\Users\PRH\Downloads\DST_020deg.csv",
    r"C:\Users\PRH\Downloads\DST_010deg.csv",
    r"C:\Users\PRH\Downloads\DST_000deg.csv"
]

soc_ocv, V_ocv_raw, Qn = get_ocv(file_charge, file_discharge)
valid_idx = ~np.isnan(V_ocv_raw)
p_coeffs_ocv = np.polyfit(soc_ocv[valid_idx], V_ocv_raw[valid_idx], 9)

all_Ts = []
all_R0 = []
all_R1 = []
all_C1 = []
all_R2 = []
all_C2 = []

# =============================================================================
# --- 3. TRAITEMENT EN BOUCLE DE CHAQUE FICHIER ---
# =============================================================================
for file_path in test_files:
    print(f"\n--- Traitement du fichier : {file_path} ---")
    
    time, V_meas, I_meas, Ts_meas, Ta_meas = load_battery_data(file_path)
    soc_true = get_soc(file_path, Qn, initial_soc)
    if np.max(soc_true) <= 1.05:
        soc_true = soc_true * 100
        
    N = len(time)

    lmbda = 0.9999
    P = 1*np.eye(5)
    theta = np.array([[0.1], [0.1], [0.01], [0.01], [0.01]]) 
    y_past, u_past = np.zeros(2), np.zeros(2)

    soc_estimated = np.zeros(N)
    soc_estimated[0] = 80 
    V_model = np.zeros(N)
    T_s_est = np.zeros(N)
    T_s_est[0] = Ts_meas[0]

    b3 = -Ts/(3600*Qn)*100
    x_k = np.array([[0],[0],[soc_estimated[0]]])
    P_KF = np.diag([1, 1, 0.5]) 
    Q = np.diag([1e-4, 1e-4, 1e-3]) 
    R_kf = np.array([[0.0001]]) 
    kappa = 1
    weight = np.array([[0.4],[0.1],[0.1],[0.1],[0.1],[0.1],[0.1]])  
    
    ocv_k = np.polyval(p_coeffs_ocv, soc_estimated[0])

    R0_temp = np.zeros(N); R1_temp = np.zeros(N); C1_temp = np.zeros(N)
    R2_temp = np.zeros(N); C2_temp = np.zeros(N)

    for k in range(N):
        u_k = I_meas[k]
        
        y_rls = ocv_k - V_meas[k]
        if k > 1: 
            phi_k = np.array([y_past[0], y_past[1], u_k, u_past[0], u_past[1]])
            if abs(u_k) > 0.05 or abs(u_k - u_past[0]) > 0.05:
                theta, P = rls_step(y_rls, phi_k, theta, P, lmbda)
     
        r0, r1, c1, r2, c2 = theta_to_2rc(theta, Ts)
        
        # =====================================================================
        # --- SATURATION DYNAMIQUE EN TEMPS RÉEL (REMPLACE LE FILTRE FINAL) ---
        # =====================================================================
        r0 = np.clip(r0, 0.005, 0.15)
        r1 = np.clip(r1, 0.001, 0.1)
        c1 = np.clip(c1, 100.0, 7000.0)
        r2 = np.clip(r2, 0.001, 0.2)
        c2 = np.clip(c2, 500.0, 40000.0)
        
        R0_temp[k], R1_temp[k], C1_temp[k], R2_temp[k], C2_temp[k] = r0, r1, c1, r2, c2

        a1, a2 = np.exp(-Ts / (r1 * c1)), np.exp(-Ts / (r2 * c2))
        A = np.array([[a1,0,0],[0,a2,0],[0,0,1]])
        B = np.array([[r1*(1-a1)], [r2*(1-a2)], [b3]])
        D = np.array([[-r0]])

        sigma_points = generate_sigma_points(x_k, P_KF, kappa)
        sigma_points_pred = np.zeros((3,7))
        for i in range(7):
            sigma_points_pred[:,i] = (A @ sigma_points[:,i].reshape(3,1) + B * u_k).flatten()  

        x_pred = sigma_points_pred @ weight
        P_pred = np.copy(Q)
        for i in range(7):  
            ecart = sigma_points_pred[:,i].reshape(3,1) - x_pred
            P_pred += weight[i,0] * (ecart @ ecart.T) 

        output_sigma_pred = np.zeros((1,7))
        for i in range(7):
            y_pred_ukf = -sigma_points_pred[0, i] - sigma_points_pred[1, i] + (D * u_k).item() + np.polyval(p_coeffs_ocv, sigma_points_pred[2, i])
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

    # Les données tracées sont maintenant propres pour chaque dataset
    filename = os.path.basename(file_path)
    fig_ds, axs_ds = plt.subplots(3, 2, figsize=(12, 8))
    fig_ds.suptitle(f"Évolution des paramètres 2RC - {filename}", fontweight='bold', fontsize=14)

    axs_ds[0, 0].plot(time, R0_temp, color='blue', alpha=0.7)
    axs_ds[0, 0].set_title('R0 en fonction du temps')
    axs_ds[0, 0].set_ylabel('R0 (Ohms)')
    axs_ds[0, 0].grid(True, linestyle=':')

    axs_ds[1, 0].plot(time, R1_temp, color='orange', alpha=0.7)
    axs_ds[1, 0].set_title('R1 en fonction du temps')
    axs_ds[1, 0].set_ylabel('R1 (Ohms)')
    axs_ds[1, 0].grid(True, linestyle=':')

    axs_ds[2, 0].plot(time, C1_temp, color='green', alpha=0.7)
    axs_ds[2, 0].set_title('C1 en fonction du temps')
    axs_ds[2, 0].set_xlabel('Temps (s)')
    axs_ds[2, 0].set_ylabel('C1 (Farads)')
    axs_ds[2, 0].grid(True, linestyle=':')

    axs_ds[0, 1].plot(time, R2_temp, color='red', alpha=0.7)
    axs_ds[0, 1].set_title('R2 en fonction du temps')
    axs_ds[0, 1].set_ylabel('R2 (Ohms)')
    axs_ds[0, 1].grid(True, linestyle=':')

    axs_ds[1, 1].plot(time, C2_temp, color='purple', alpha=0.7)
    axs_ds[1, 1].set_title('C2 en fonction du temps')
    axs_ds[1, 1].set_xlabel('Temps (s)')
    axs_ds[1, 1].set_ylabel('C2 (Farads)')
    axs_ds[1, 1].grid(True, linestyle=':')

    axs_ds[2, 1].axis('off')
    
    plt.tight_layout()
    plt.show(block=False)

    # Calcul simple des médianes à partir de t > 50s (plus de masques complexes nécessaires)
    valid_time = time > 50
    
    all_Ts.append(np.median(Ts_meas[valid_time]))
    all_R0.append(np.median(R0_temp[valid_time]))
    all_R1.append(np.median(R1_temp[valid_time]))
    all_C1.append(np.median(C1_temp[valid_time]))
    all_R2.append(np.median(R2_temp[valid_time]))
    all_C2.append(np.median(C2_temp[valid_time]))

print("\nTous les fichiers ont été traités !")

all_Ts = np.array(all_Ts)
all_R0 = np.array(all_R0)
all_R1 = np.array(all_R1)
all_C1 = np.array(all_C1)
all_R2 = np.array(all_R2)
all_C2 = np.array(all_C2)

# =============================================================================
# --- 4. IDENTIFICATION DES MODÈLES (AUCUNE RETOUCHE GLOBALE NÉCESSAIRE) ---
# =============================================================================
coeffs_R0 = np.polyfit(all_Ts, np.log(all_R0), 1)
B_R0 = coeffs_R0[0]
A_R0 = np.exp(coeffs_R0[1])

coeffs_R1 = np.polyfit(all_Ts, all_R1, 2)
coeffs_C1 = np.polyfit(all_Ts, all_C1, 3)
coeffs_R2 = np.polyfit(all_Ts, all_R2, 2)
coeffs_C2 = np.polyfit(all_Ts, all_C2, 3)

print("\n" + "="*50)
print("ÉQUATIONS DU MODÈLE THERMO-ÉLECTRIQUE")
print("="*50)
print(f"R0(Ts) = {A_R0:.4e} * exp({B_R0:.4f} * Ts)")
print(f"R1(Ts) = {coeffs_R1[0]:.6e}*Ts^2 + {coeffs_R1[1]:.6e}*Ts + {coeffs_R1[2]:.4e}")
print(f"C1(Ts) = {coeffs_C1[0]:.4e}*Ts^3 + {coeffs_C1[1]:.4e}*Ts^2 + {coeffs_C1[2]:.4e}*Ts + {coeffs_C1[3]:.4e}")
print(f"R2(Ts) = {coeffs_R2[0]:.6e}*Ts^2 + {coeffs_R2[1]:.6e}*Ts + {coeffs_R2[2]:.4e}")
print(f"C2(Ts) = {coeffs_C2[0]:.4e}*Ts^3 + {coeffs_C2[1]:.4e}*Ts^2 + {coeffs_C2[2]:.4e}*Ts + {coeffs_C2[3]:.4e}")
print("="*50)

# =============================================================================
# --- 5. AFFICHAGE DES RÉSULTATS GLOBAUX ---
# =============================================================================
fig, axs = plt.subplots(3, 2, figsize=(14, 10))
fig.suptitle("Composants 2RC (Médianes par Dataset & Polynômes Stables)", fontweight='bold', fontsize=14)

Ts_axis = np.linspace(np.min(all_Ts), np.max(all_Ts), 200)

axs[0, 0].scatter(all_Ts, all_R0, s=50, color='blue', label='Médiane mesurée')
axs[0, 0].plot(Ts_axis, A_R0 * np.exp(B_R0 * Ts_axis), color='red', linewidth=2, label='Fit Arrhenius')
axs[0, 0].set_title('R0 = f(Ts)')
axs[0, 0].set_ylabel('R0 (Ohms)')
axs[0, 0].legend()
axs[0, 0].grid(True, linestyle=':')

axs[1, 0].scatter(all_Ts, all_R1, s=50, color='orange', label='Médiane mesurée')
axs[1, 0].plot(Ts_axis, np.polyval(coeffs_R1, Ts_axis), color='black', linewidth=2)
axs[1, 0].set_title('R1 = f(Ts)')
axs[1, 0].grid(True, linestyle=':')
axs[1, 0].legend()

axs[2, 0].scatter(all_Ts, all_C1, s=50, color='green', label='Médiane mesurée')
axs[2, 0].plot(Ts_axis, np.polyval(coeffs_C1, Ts_axis), color='black', linewidth=2)
axs[2, 0].set_title('C1 = f(Ts)')
axs[2, 0].set_xlabel('Température (°C)')
axs[2, 0].grid(True, linestyle=':')
axs[2, 0].legend()

axs[0, 1].scatter(all_Ts, all_R2, s=50, color='red', label='Médiane mesurée')
axs[0, 1].plot(Ts_axis, np.polyval(coeffs_R2, Ts_axis), color='black', linewidth=2)
axs[0, 1].set_title('R2 = f(Ts)')
axs[0, 1].grid(True, linestyle=':')
axs[0, 1].legend()

axs[1, 1].scatter(all_Ts, all_C2, s=50, color='purple', label='Médiane mesurée')
axs[1, 1].plot(Ts_axis, np.polyval(coeffs_C2, Ts_axis), color='black', linewidth=2)
axs[1, 1].set_title('C2 = f(Ts)')
axs[1, 1].set_xlabel('Température (°C)')
axs[1, 1].grid(True, linestyle=':')
axs[1, 1].legend()

axs[2, 1].axis('off') 
plt.tight_layout()

plt.show()