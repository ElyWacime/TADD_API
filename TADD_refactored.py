from refactored_code import calculate_n_h90_for_config
from refactored_code import calculate_P_ac_for_config, calculate_inverter_load_for_config
from refactored_code import calculate_P_dc_for_config
from refactored_code import modules, module_parameters, module_area, eta_module, calculate_single_res_for_config, calculate_Itot_for_config, calculate_IL_and_IO_and_Rs_and_Rsh_and_nNsVth_for_config, calculate_WS5m_and_Ta_and_Tcell_for_config, calculate_effactive_irr_for_config, calculate_IAM_for_config, calculate_AOI_for_config, calculate_TMY_global_for_config
from pvlib import location
import pandas as pd
import numpy as np
import pvlib

lat = 43.5161
long = -1.0285

TMY_data = pvlib.iotools.get_pvgis_tmy(lat, long, outputformat='csv', usehorizon=True,
                                       userhorizon=None, startyear=2005, endyear=2020, url='https://re.jrc.ec.europa.eu/api/v5_2/',
                                       map_variables=False, timeout=30)
TMY_data = TMY_data[0]
# Remove time zone to save data in excel
TMY_data.set_index(TMY_data.index.tz_localize(None), inplace=True, drop=True)
TMY_global = TMY_data


a = -4.56  # DEFAUT = -3.56
b = -0.085  # DEFAUT = -0.075
deltaT = 5

gamma_mp = -0.9

n = 1.526
L = 0.002
K = 4

rho = 0.15

soiling = 0.9

module_eff = 0.5

mismatch = 0.03  # DEFAUT = 0.02
connections = 0.015  # DEFAUT = 0.005
DCwiring = 0.025  # DEFAUT = 0.015
inverterloss = 0.025  # DEFAUT = 0.015
inverter_min = 0.03

ACwiring = 0.01

r_DCAC = 1.2

bifac = 1
bifac_ratio = 0.65
H = 4
inter = 2.5

r_P90 = 0.94


toiture_surface1 = 0
toiture_surface2 = 0
toiture_surface3 = 0
toiture_surface4 = 0
toiture_surface5 = 0
toiture_surface6 = 0

serre_surface1 = 0
serre_surface2 = 0
serre_surface3 = 0

ombriere_surface1 = 1
ombriere_surface2 = 0
ombriere_surface3 = 0


surface_azimuth_toitureS1 = 0
surface_azimuth_toitureS2 = 0
surface_azimuth_toitureS3 = 0
surface_azimuth_toitureS4 = 0
surface_azimuth_toitureS5 = 0
surface_azimuth_toitureS6 = 0

surface_tilt_toitureS1 = 0
surface_tilt_toitureS2 = 0
surface_tilt_toitureS3 = 0
surface_tilt_toitureS4 = 0
surface_tilt_toitureS5 = 0
surface_tilt_toitureS6 = 0


bifac_toitureS1 = 0
bifac_toitureS2 = 0
bifac_toitureS3 = 0
bifac_toitureS4 = 0
bifac_toitureS5 = 0
bifac_toitureS6 = 0

bifac_ratio_toitureS1 = 0.65
bifac_ratio_toitureS2 = 0.65
bifac_ratio_toitureS3 = 0.65
bifac_ratio_toitureS4 = 0.65
bifac_ratio_toitureS5 = 0.65
bifac_ratio_toitureS6 = 0.65

H_toitureS1 = 4
H_toitureS2 = 4
H_toitureS3 = 4
H_toitureS4 = 4
H_toitureS5 = 4
H_toitureS6 = 4

inter_toitureS1 = 2.5
inter_toitureS2 = 2.5
inter_toitureS3 = 2.5
inter_toitureS4 = 2.5
inter_toitureS5 = 2.5
inter_toitureS6 = 2.5

surface_azimuth_ombS1 = 90
surface_azimuth_ombS2 = 180
surface_azimuth_ombS3 = 180

surface_tilt_ombS1 = 30
surface_tilt_ombS2 = 0
surface_tilt_ombS3 = 0

bifac_ombS1 = 1
bifac_ombS2 = 0
bifac_ombS3 = 0

bifac_ratio_ombS1 = 0.65
bifac_ratio_ombS2 = 0.65
bifac_ratio_ombS3 = 0.65

H_ombS1 = 4
H_ombS2 = 4
H_ombS3 = 4

inter_ombS1 = 2.5
inter_ombS2 = 2.5
inter_ombS3 = 2.5

surface_azimuth_serS1 = 180
surface_azimuth_serS2 = 180
surface_azimuth_serS3 = 180

surface_tilt_serS1 = 0
surface_tilt_serS2 = 0
surface_tilt_serS3 = 0

bifac_serS1 = 0
bifac_serS2 = 0
bifac_serS3 = 0

bifac_ratio_serS1 = 0.65
bifac_ratio_serS2 = 0.65
bifac_ratio_serS3 = 0.65

H_serS1 = 4
H_serS2 = 4
H_serS3 = 4

inter_serS1 = 2.5
inter_serS2 = 2.5
inter_serS3 = 2.5

surfaces_configurations = [
    # Configurations for toiture
    {"name": "toiture_surface1", "surface": toiture_surface1,  "surface_azimuth": surface_azimuth_toitureS1, "surface_tilt": surface_tilt_toitureS1, "bifac": bifac_toitureS1,
        "bifac_ratio": bifac_ratio_toitureS1, "H": H_toitureS1, "inter": inter_toitureS1},
    {"name": "toiture_surface2", "surface": toiture_surface2,  "surface_azimuth": surface_azimuth_toitureS2, "surface_tilt": surface_tilt_toitureS2, "bifac": bifac_toitureS2,
        "bifac_ratio": bifac_ratio_toitureS2, "H": H_toitureS2, "inter": inter_toitureS2},
    {"name": "toiture_surface3", "surface": toiture_surface3,  "surface_azimuth": surface_azimuth_toitureS3, "surface_tilt": surface_tilt_toitureS3, "bifac": bifac_toitureS3,
        "bifac_ratio": bifac_ratio_toitureS3, "H": H_toitureS3, "inter": inter_toitureS3},
    {"name": "toiture_surface4", "surface": toiture_surface4,  "surface_azimuth": surface_azimuth_toitureS4, "surface_tilt": surface_tilt_toitureS4, "bifac": bifac_toitureS4,
        "bifac_ratio": bifac_ratio_toitureS4, "H": H_toitureS4, "inter": inter_toitureS4},
    {"name": "toiture_surface5", "surface": toiture_surface5,  "surface_azimuth": surface_azimuth_toitureS5, "surface_tilt": surface_tilt_toitureS5, "bifac": bifac_toitureS5,
        "bifac_ratio": bifac_ratio_toitureS5, "H": H_toitureS5, "inter": inter_toitureS5},
    {"name": "toiture_surface6", "surface": toiture_surface6,  "surface_azimuth": surface_azimuth_toitureS6, "surface_tilt": surface_tilt_toitureS6, "bifac": bifac_toitureS6,
        "bifac_ratio": bifac_ratio_toitureS6, "H": H_toitureS6, "inter": inter_toitureS6},
    # Configurations for ombrière
    {"name": "ombriere_surface1", "surface": ombriere_surface1, "surface_azimuth": surface_azimuth_ombS1,
        "surface_tilt": surface_tilt_ombS1, "bifac": bifac_ombS1, "bifac_ratio": bifac_ratio_ombS1, "H": H_ombS1, "inter": inter_ombS1},
    {"name": "ombriere_surface2", "surface": ombriere_surface2, "surface_azimuth": surface_azimuth_ombS2,
        "surface_tilt": surface_tilt_ombS2, "bifac": bifac_ombS2, "bifac_ratio": bifac_ratio_ombS2, "H": H_ombS2, "inter": inter_ombS2},
    {"name": "ombriere_surface3", "surface": ombriere_surface3, "surface_azimuth": surface_azimuth_ombS3,
        "surface_tilt": surface_tilt_ombS3, "bifac": bifac_ombS3, "bifac_ratio": bifac_ratio_ombS3, "H": H_ombS3, "inter": inter_ombS3},
    # Configurations for serre
    {"name": "serre_surface1", "surface": serre_surface1, "surface_azimuth": surface_azimuth_serS1,
        "surface_tilt": surface_tilt_serS1, "bifac": bifac_serS1, "bifac_ratio": bifac_ratio_serS1, "H": H_serS1, "inter": inter_serS1},
    {"name": "serre_surface2", "surface": serre_surface2, "surface_azimuth": surface_azimuth_serS2,
        "surface_tilt": surface_tilt_serS2, "bifac": bifac_serS2, "bifac_ratio": bifac_ratio_serS2, "H": H_serS2, "inter": inter_serS2},
    {"name": "serre_surface3", "surface": serre_surface3, "surface_azimuth": surface_azimuth_serS3,
        "surface_tilt": surface_tilt_serS3, "bifac": bifac_serS3, "bifac_ratio": bifac_ratio_serS3, "H": H_serS3, "inter": inter_serS3},
]

times = pd.date_range('2005-01-01', periods=8760, freq='1H', tz='Europe/Paris')
loc = location.Location(lat, long, tz=times.tz)
sp = loc.get_solarposition(times)
# Récupération des angles solaires via le module 'loc' de pvlib
solar_zenith = sp['apparent_zenith'].values
solar_azimuth = sp['azimuth'].values

dni = TMY_global['Gb(n)'].values
ghi = TMY_global['G(h)'].values
dhi = TMY_global['Gd(h)'].values
# elimination des valeurs nulles ou négatives de la DHI pour le modèle de transposition PEREZ
dhi[dhi <= 0] = 0.001

dni_extra = pvlib.irradiance.get_extra_radiation(
    times, solar_constant=1366.1, method='spencer', epoch_year=2023)
AM = pvlib.atmosphere.get_relative_airmass(
    zenith=sp['zenith'], model='kastenyoung1989')

Itot_for_config = calculate_Itot_for_config(solar_zenith=solar_zenith,
                                            solar_azimuth=solar_azimuth,
                                            dni=dni, ghi=ghi, dhi=dhi, dni_extra=dni_extra,
                                            AM=AM, rho=rho,
                                            surfaces_configurations=surfaces_configurations)

AOI_for_config = calculate_AOI_for_config(solar_zenith=solar_zenith, solar_azimuth=solar_azimuth,
                                          surfaces_configurations=surfaces_configurations)

IAM_for_config = calculate_IAM_for_config(AOI_for_config=AOI_for_config,
                                          n=n, K=K, L=L,
                                          surfaces_configurations=surfaces_configurations)

TMY_global_for_config = calculate_TMY_global_for_config(TMY_global=TMY_global,
                                                        Itot_for_config=Itot_for_config,
                                                        IAM_for_config=IAM_for_config,
                                                        surfaces_configurations=surfaces_configurations)

effective_irr_for_config = calculate_effactive_irr_for_config(TMY_global_for_config=TMY_global_for_config,
                                                              soiling=soiling, surfaces_configurations=surfaces_configurations)

(WS5m_for_config, Ta_for_config, Tcell_for_config) = calculate_WS5m_and_Ta_and_Tcell_for_config(TMY_global_for_config=TMY_global_for_config,
                                                                                                effective_irr_for_config=effective_irr_for_config,
                                                                                                a=a, b=b, deltaT=deltaT,
                                                                                                surfaces_configurations=surfaces_configurations)

(IL_for_config, I0_for_config, Rs_for_config, Rsh_for_config, nNsVth_for_config) = calculate_IL_and_IO_and_Rs_and_Rsh_and_nNsVth_for_config(effective_irr_for_config=effective_irr_for_config,
                                                                                                                                            Tcell_for_config=Tcell_for_config,
                                                                                                                                            module_parameters=module_parameters,
                                                                                                                                            surfaces_configurations=surfaces_configurations)

single_res_for_config = calculate_single_res_for_config(IL_for_config=IL_for_config, I0_for_config=I0_for_config,
                                                        Rs_for_config=Rs_for_config, Rsh_for_config=Rsh_for_config,
                                                        nNsVth_for_config=nNsVth_for_config,
                                                        surfaces_configurations=surfaces_configurations)

P_dc_for_config = calculate_P_dc_for_config(single_res_for_config=single_res_for_config,
                                            mismatch=mismatch, connections=connections,
                                            DCwiring=DCwiring, module_area=module_area,
                                            surfaces_configurations=surfaces_configurations)

inverter_load_for_config = calculate_inverter_load_for_config(P_dc_for_config=P_dc_for_config,
                                                              eta_module=eta_module,
                                                              r_DCAC=r_DCAC, inverterloss=inverterloss,
                                                              surfaces_configurations=surfaces_configurations)

P_ac_for_config, sum_P_ac_for_config = calculate_P_ac_for_config(inverterloss=inverterloss, P_dc_for_config=P_dc_for_config,
                                                                 inverter_load_for_config=inverter_load_for_config,
                                                                 inverter_min=inverter_min, r_DCAC=r_DCAC, ACwiring=ACwiring,
                                                                 surfaces_configurations=surfaces_configurations)

n_h90_for_config, gain = calculate_n_h90_for_config(P_ac_for_config=P_ac_for_config, r_P90=r_P90,
                                                    eta_module=eta_module, module_eff=module_eff,
                                                    surfaces_configurations=surfaces_configurations,
                                                    bifac=bifac, bifac_ratio=bifac_ratio, H=H, inter=inter,
                                                    rho=rho)


somme = 0
for elem in sum_P_ac_for_config:
    somme += elem

data = {"n_h90_for_config": n_h90_for_config}

"""pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)"""
print(f'TADD: {somme}')
