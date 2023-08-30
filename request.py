import requests
import json

configuration_data = {
    "long": -1.0285,
    "lat": 43.5161,
    "ombriere_surface1": 1,
    "surface_azimuth_ombS1": 90,
    "surface_tilt_ombS1": 30,
    "module_eff": 0.5,
    "a": -4.56,
    "b": -0.085,
    "deltaT": 5,
    "gamma_mp": -0.9,
    "n": 1.526,
    "L": 0.002,
    "K": 4,
    "rho": 0.15,
    "soiling": 0.9,
    "mismatch": 0.03,  # DEFAUT = 0.02
    "connections": 0.015,  # DEFAUT = 0.005
    "DCwiring": 0.025,  # DEFAUT = 0.015
    "inverterloss": 0.025,  # DEFAUT = 0.015
    "inverter_min": 0.03,
    "ACwiring": 0.01,  # DEFAUT = 0.01

    # r_DCAC ratio
    "r_DCAC": 1.2,  # DEFAUT = 1.2

    # Inputs pour Ombrière ou Serre
    "bifac_ombS1": 1,  # Si bifac = 0, module non bifacial, sinon = 1 si bifacial. DEFAUT = 0
    "bifac_ratio_ombS1": 0.65,  # DEFAUT = 0.65
    "H_ombS1": 4,  # DEFAUT = 4m
    "inter_ombS1": 2.5,  # DEFAUT = 2.5m

    # Ratio P90/P50
    "r_P90": 0.94
}


api_endpoint = 'http://api.enervivo.fr/calculate_pv'

response = requests.post(api_endpoint, json=configuration_data)

if response.status_code == 200:

    result = response.json()
    somme = 0
    for elem in (result["P_ac_for_config"]):
        somme += elem

    print(somme)
else:
    print("Error occurred. Status code:", response.status_code)
    print(response.text)
