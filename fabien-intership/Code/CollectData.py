import time
import pyvisa as visa
import threading

# ==========================================
# 1. PARAMÈTRES GLOBAUX
# ==========================================
T_SECONDS = 1  # Remplacez par votre intervalle T en secondes
POWER_SUPPLY_ID = 'USB0::0x2EC7::0x6700::802259073777170159::INSTR'
ELECTRONIC_LOAD_ID = 'USB0::0x1AB1::0x0E11::DL3A250700137::INSTR'

# Variable globale pour arrêter la boucle proprement
killThread = False 

# ==========================================
# 2. FONCTION DE LECTURE MATÉRIELLE
# ==========================================
def configMeasureQuery(device, measure):
    """Envoie la commande SCPI pour lire la tension ou le courant."""
    if measure == "VOLT":
        return device.query('MEASure:VOLTage?')
    if measure == "CURR":
        return device.query('MEASure:CURRent?')

# ==========================================
# 3. FONCTION DE COLLECTE DE DONNÉES
# ==========================================
def collect_data(ps_device, el_device):
    """Interroge les deux appareils et récupère les valeurs."""
    try:
        # Lecture Alimentation (Power Supply)
        ps_volt = float(configMeasureQuery(ps_device, "VOLT"))
        ps_curr = float(configMeasureQuery(ps_device, "CURR"))
        
        # Lecture Charge Électronique (Electronic Load)
        el_volt = float(configMeasureQuery(el_device, "VOLT"))
        el_curr = float(configMeasureQuery(el_device, "CURR"))
        
        # Affichage des données
        print(f"Power Supply   -> Tension: {ps_volt:.3f} V | Courant: {ps_curr:.3f} A")
        print(f"Electronic Load-> Tension: {el_volt:.3f} V | Courant: {el_curr:.3f} A")
        print("-" * 50)
        
    except Exception as e:
        print(f"Erreur de lecture: {e}")

# ==========================================
# 4. BOUCLE TEMPORELLE (Le "Counter")
# ==========================================
def data_logger_thread(ps_device, el_device, t_seconds):
    """Exécute la collecte toutes les T secondes sans dérive temporelle."""
    global killThread
    
    # 1. Établir l'ancre temporelle de départ
    next_call = time.time()

    while not killThread:
        # 2. Avancer l'horloge de T secondes
        next_call = next_call + t_seconds
        
        # 3. Collecter les données
        collect_data(ps_device, el_device)
        
        # 4. Calculer le temps restant à dormir
        time_to_sleep = next_call - time.time()
        
        # 5. Dormir uniquement si on n'est pas en retard
        if time_to_sleep > 0:
            time.sleep(time_to_sleep)

# ==========================================
# 5. EXÉCUTION PRINCIPALE
# ==========================================
if __name__ == "__main__":
    rm = visa.ResourceManager()
    
    try:
        # Ouverture des connexions VISA
        ps_device = rm.open_resource(POWER_SUPPLY_ID)
        el_device = rm.open_resource(ELECTRONIC_LOAD_ID)
        
        print(f"Démarrage de la collecte toutes les {T_SECONDS} secondes...")
        print("Appuyez sur Ctrl+C pour arrêter.\n")
        
        # Lancement de la boucle dans un thread séparé
        logger = threading.Thread(target=data_logger_thread, args=(ps_device, el_device, T_SECONDS))
        logger.start()
        
        # Maintenir le programme principal en vie
        while True:
            time.sleep(0.5)
            
    except KeyboardInterrupt:
        # Arrêt propre via Ctrl+C
        print("\nArrêt demandé par l'utilisateur...")
        killThread = True
        logger.join() # Attendre que le thread se termine
        print("Collecte terminée.")
        
    except Exception as err:
        print(f"Erreur d'initialisation VISA : {err}")