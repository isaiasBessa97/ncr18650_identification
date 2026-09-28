import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.linalg import cholesky

# =============================================================================
# --- FONCTIONS DE LECTURE UNIVERSELLES ---
# =============================================================================
def load_battery_data(file_path):
    """
    Fonction universelle pour charger les fichiers .txt ou .csv.
    Détecte automatiquement le séparateur et les colonnes.
    """
    df = pd.read_csv(file_path, sep=None, engine='python')
    
    # Détection du nouveau format CSV 
    if 'Tempo_s' in df.columns:
        time = df['Tempo_s'].values
        I_meas = df['Corrente_A'].values
        V_meas = df['Tensao_V'].values
        Ts_meas = df['Temperatura_C'].values
        # Température ambiante supposée égale à la température initiale de la cellule
        Ta_meas = np.full(len(time), Ts_meas[0])
        
    # Ancien format TXT
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
        return 0.05, 0.01, 5000, 0.01, 40000 # Valeurs refuges

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
    
    num_points = len(I_meas)
    soc_array = np.zeros(num_points)
    soc_array[0] = initial_soc
    
    for ii in range(1, num_points):
        dt = time[ii] - time[ii-1] if time[ii] - time[ii-1] != 0 else 1
        soc_array[ii] = soc_array[ii-1] - (100 * dt / (3600 * Qn)) * I_meas[ii]
        
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
# --- 1. Configuration et Chargement ---
# =============================================================================
Ts = 1.0  
Qn = 3.05 
initial_soc = 100

file_charge = r"C:\ncr18650_identification\dataset-thermal\BID003\BID003_CCCV005.0_02022026.txt"
file_discharge = r"C:\ncr18650_identification\dataset-thermal\BID003\BID003_CDch005.0_02022026.txt"
file_test = r"C:\Users\PRH\Downloads\MPDch_045deg.csv"

soc_ocv, V_ocv_raw, _ = get_ocv(file_charge, file_discharge)
valid_idx = ~np.isnan(V_ocv_raw)
p_coeffs_ocv = np.polyfit(soc_ocv[valid_idx], V_ocv_raw[valid_idx], 9)
dp_coeffs_ocv = np.polyder(p_coeffs_ocv)

soc_true = get_soc(file_test, Qn, initial_soc)
if np.max(soc_true) <= 1.05:
    soc_true = soc_true * 100

# Chargement intelligent (CSV ou TXT)
time, V_meas, I_meas, Ts_meas, Ta_meas = load_battery_data(file_test)
N = len(time)

# =============================================================================
# --- 2. Initialisation RLS & Kalman ---
# =============================================================================
lmbda = 0.9999
P = 1*np.eye(5)
theta = np.array([[0.1], [0.1], [0.01], [0.01], [0.01]]) 

soc_estimated = np.zeros(N)
soc_estimated[0] = 100 
theta_history = np.zeros((N, 5))
V_model = np.zeros(N)

y_past = np.zeros(2)
u_past = np.zeros(2)

b3 = -Ts/(3600*Qn)*100
x_k = np.array([[0],[0],[soc_estimated[0]]])
A = np.array([[0,0,0],[0,0,0],[0,0,0]])
B = np.array([[0],[0],[b3]])
D = np.array([[-0.1]]) 
P_KF = np.diag([1, 1, 0.5]) 
P_zn = np.array([[0.]])
P_xz = np.array([[0.],[0.],[0.]])
Q = np.diag([1e-4, 1e-4, 1e-3]) 
R_kf = np.array([[0.0001]]) 
kappa = 1
weight = np.array([[0.4],[0.1],[0.1],[0.1],[0.1],[0.1],[0.1]])  
soc_estimated[0] = max(0, min(100, x_k[2, 0]))
Kn = np.array([[0],[0],[0]])
        
ocv_k = np.polyval(p_coeffs_ocv, soc_estimated[0])

# =============================================================================
# --- 3. Boucle Temps Réel ---
# =============================================================================
print('Starting Real-Time RLS with Coulomb Counting...')

for k in range(N):
    u_k = I_meas[k]
    y_rls = ocv_k - V_meas[k]
    
    if k > 1: 
        phi_k = np.array([y_past[0], y_past[1], u_k, u_past[0], u_past[1]])
        if abs(u_k) > 0.05 or abs(u_k - u_past[0]) > 0.05:
            theta, P = rls_step(y_rls, phi_k, theta, P, lmbda)
 
    theta_history[k, :] = theta.flatten()
    r0, r1, c1, r2, c2 = theta_to_2rc(theta, Ts)

    r0 = max(r0, 1e-4)
    r1 = max(r1, 1e-4)
    c1 = max(c1, 1.0)   
    r2 = max(r2, 1e-4)
    c2 = max(c2, 1.0)   

    a1 = np.exp(-Ts / (r1 * c1))
    a2 = np.exp(-Ts / (r2 * c2))
    b1 = r1 * (1 - a1)
    b2 = r2 * (1 - a2)

    A = np.array([[a1,0,0],[0,a2,0],[0,0,1]])
    B = np.array([[b1],[b2],[b3]])
    D = np.array([[-r0]])

    sigma_points_pred = np.zeros((3,7))
    sigma_points = generate_sigma_points(x_k,P_KF, kappa)
    
    for i in range(7):
        x_sigma_pred = A @ sigma_points[:,i].reshape(3,1) + B * u_k
        sigma_points_pred[:,i] =  x_sigma_pred.flatten()  
    
    x_pred = sigma_points_pred @ weight
    P_pred = np.copy(Q)

    for i in range(7):  
        ecart = sigma_points_pred[:,i].reshape(3,1) - x_pred
        P_pred += weight[i,0] * (ecart @ ecart.T) 

    output_sigma_points_pred = np.zeros((1,7))
    for i in range(7):
        ocv_point = np.polyval(p_coeffs_ocv, sigma_points_pred[2, i])
        y_pred_ukf = -sigma_points_pred[0, i] - sigma_points_pred[1, i] + (D * u_k).item() + ocv_point
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

print('Simulation Finished!')

# =============================================================================
# 4. CALCUL DES ERREURS (RMSE) ET AFFICHAGE DES RÉSULTATS
# =============================================================================
valid_idx = time > 10 
rmse_V = np.sqrt(np.mean((V_meas[valid_idx] - V_model[valid_idx])**2))
rmse_soc = np.sqrt(np.mean((soc_true[valid_idx] - soc_estimated[valid_idx])**2))

print(f"\nPerformances du Filtre de Kalman (après 10s) :")
print(f" -> RMSE Tension : {rmse_V:.4f} V")
print(f" -> RMSE SoC     : {rmse_soc:.2f} %")

# --- Configuration globale de Matplotlib pour imiter le style MATLAB ---
plt.rcParams.update({
    "text.usetex": True,           # Interprète LaTeX
    "font.family": "serif",        # Police classique de LaTeX
    "font.size": 16,               # Taille de police à 16
    "axes.grid": True,             # Activation de la grille (grid on)
    "grid.color": "#b0b0b0",       # Couleur de grille standard
    "grid.linestyle": "-",
    "grid.linewidth": 0.5,
    "legend.edgecolor": "black",   # Contour de légende noir
    "legend.fancybox": False,      # Bords carrés pour la légende
    "legend.framealpha": 1.0       # Fond blanc opaque pour la légende
})

# --- Graphique 1 : Comparaison de la Tension (UKF) ---
plt.figure(figsize=(10, 6))
# LineWidth 2 et couleurs k (noir) et r (rouge)
plt.plot(time, V_meas, 'k', linewidth=2, label=r'Tension Measured (Expérimentale)')
plt.plot(time, V_model, 'r', linewidth=2, label=r'Tension model (Kalman)')

plt.title(fr'Voltage Comparaison: Measured vs Model Kalman (RMSE = {rmse_V:.4f} V)')
plt.xlabel(r'Temps (s)')
plt.ylabel(r'Tension (V)')
plt.legend(loc='best')
plt.tight_layout()

# --- Graphique 2 : Comparaison du SoC (UKF) ---
plt.figure(figsize=(10, 6))
# LineWidth 2 et couleurs k (noir) et b (bleu)
plt.plot(time, soc_true, 'k', linewidth=2, label=r'SoC REAL (Intégration Théorique)')
plt.plot(time, soc_estimated, 'b', linewidth=2, label=r'SoC Model (Filtre Kalman)')

plt.title(fr'State of charge (SoC) : Real vs Estimated (RMSE = {rmse_soc:.2f} \%)')
plt.xlabel(r'Temps (s)')
plt.ylabel(r'State of Charge (\%)')
plt.legend(loc='best')
plt.tight_layout()

plt.show()