import pyvisa

rm = pyvisa.ResourceManager()
ports = ['ASRL1::INSTR', 'ASRL5::INSTR']

for port in ports:
    try:
        # Ouverture de la connexion
        inst = rm.open_resource(port)
        
        # Configuration courante pour les liaisons séries (ajustez le baud_rate si besoin)
        inst.baud_rate = 9600 
        inst.read_termination = '\n'
        inst.write_termination = '\n'
        inst.timeout = 2000  # 2 secondes de délai max

        # Demande d'identification
        idn = inst.query('*IDN?')
        print(f"✅ {port} -> {idn.strip()}")
        inst.close()
    except Exception as e:
        print(f"❌ {port} -> Pas de réponse (Erreur ou mauvais baud rate)")