#!/usr/bin/env python
# coding: utf-8

# # Initialization

# ## Location

# In[1]:


from refactored_code import calculate_inverter_load_for_config
from refactored_code import calculate_P_ac_for_config
from refactored_code import calculate_P_dc_for_config
from refactored_code import calculate_single_res_for_config
from refactored_code import calculate_IL_and_IO_and_Rs_and_Rsh_and_nNsVth_for_config
from refactored_code import calculate_WS5m_and_Ta_and_Tcell_for_config
from refactored_code import calculate_effactive_irr_for_config
from refactored_code import calculate_TMY_global_for_config
from refactored_code import calculate_IAM_for_config
from refactored_code import calculate_AOI_for_config
from refactored_code import calculate_Itot_for_config
from TADD_refactored import surfaces_configurations
from pvlib import location
import pandas as pd
import numpy as np
import pvlib
lat = 43.5161  # Coordonées du projet
long = -1.0285
# La librairie PVLib possède une fonction qui pemet de récupérer les données PVGIS via son API
TMY_data = pvlib.iotools.get_pvgis_tmy(lat, long, outputformat='csv', usehorizon=True,
                                       userhorizon=None, startyear=2005, endyear=2020, url='https://re.jrc.ec.europa.eu/api/v5_2/',
                                       map_variables=False, timeout=30)
TMY_data = TMY_data[0]
# Remove time zone to save data in excel
TMY_data.set_index(TMY_data.index.tz_localize(None), inplace=True, drop=True)
TMY_global = TMY_data  # Les données météo sont stockées dans la dataframe TMY_global

# ## Input data²

# In[2]:


# Surface Azimuth (deg)
surface_azimuth = 90

# Surface tilt (deg)
surface_tilt = 30

# Module efficiency
module_eff = 0.5  # Default 20%

# Temperature coefficients
a = -4.56  # DEFAUT = -3.56
b = -0.085  # DEFAUT = -0.075
deltaT = 5  # DEFAUT = 3

# Module temperature (@NOCT, %W/°C)
gamma_mp = -0.9  # DEFAUT = -0.45

# Material proprieties
n = 1.526  # DEFAUT = 1.526
L = 0.002  # DEFAUT = 0.002
K = 4  # DEFAUT = 4

# Ground albedo
rho = 0.15  # DEFAUT = 0.2

# Soiling losses
soiling = 0.9  # DEFAUT = 0.03

# DC losses
mismatch = 0.03  # DEFAUT = 0.02
connections = 0.015  # DEFAUT = 0.005
DCwiring = 0.025  # DEFAUT = 0.015
inverterloss = 0.025  # DEFAUT = 0.015
inverter_min = 0.03  # DEFAUT = 0.01

# AC losses
ACwiring = 0.01  # DEFAUT = 0.01

# r_DCAC ratio
r_DCAC = 1.2  # DEFAUT = 1.2

# Inputs pour Ombrière ou Serre
bifac = 1  # Si bifac = 0, module non bifacial, sinon = 1 si bifacial. DEFAUT = 0
bifac_ratio = 0.65  # DEFAUT = 0.65
H = 4  # DEFAUT = 4m
inter = 2.5  # DEFAUT = 2.5m

# Ratio P90/P50
r_P90 = 0.94  # DEFAUT 0.94


# # Import Solar Pos.

# In[3]:


times = pd.date_range('2005-01-01', periods=8760, freq='1H', tz='Europe/Paris')
loc = location.Location(lat, long, tz=times.tz)
sp = loc.get_solarposition(times)
# Récupération des angles solaires via le module 'loc' de pvlib
solar_zenith = sp['apparent_zenith'].values
solar_azimuth = sp['azimuth'].values


# # Get Itot using Perez

# In[4]:


dni = TMY_global['Gb(n)'].values
ghi = TMY_global['G(h)'].values
dhi = TMY_global['Gd(h)'].values
# elimination des valeurs nulles ou négatives de la DHI pour le modèle de transposition PEREZ
dhi[dhi <= 0] = 0.001

dni_extra = pvlib.irradiance.get_extra_radiation(
    times, solar_constant=1366.1, method='spencer', epoch_year=2023)
AM = pvlib.atmosphere.get_relative_airmass(
    zenith=sp['zenith'], model='kastenyoung1989')
# Calcul de la radiation globale et de ses composantes selon le modèle de PEREZ et la librairie pvlib
Itot = pvlib.irradiance.get_total_irradiance(surface_tilt, surface_azimuth, solar_zenith, solar_azimuth,
                                             dni, ghi, dhi, dni_extra, AM, albedo=rho, surface_type=None,
                                             model='perez', model_perez='allsitescomposite1990')
Itot = pd.DataFrame(Itot)


Itot_for_config = calculate_Itot_for_config(solar_azimuth=solar_azimuth, solar_zenith=solar_zenith, dni=dni,
                                            ghi=ghi, dhi=dhi, dni_extra=dni_extra, AM=AM, rho=rho,
                                            surfaces_configurations=surfaces_configurations)

"""for config in surfaces_configurations:
    if "surface_azimuth" in config:
        if (Itot_for_config[config["name"]].reset_index(drop=True)).equals(Itot.reset_index(drop=True)) is False:
            print(f'Found it: {config["name"]}')
        else:
            print("Nope!")"""

# # Get effective irradiance

# In[5]:


"""pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)"""

AOI = pvlib.irradiance.aoi(
    surface_tilt, surface_azimuth, solar_zenith, solar_azimuth)
IAM = pvlib.iam.physical(AOI, n, K, L)
TMY_global['IAM'] = IAM


AOI_for_config = calculate_AOI_for_config(solar_zenith=solar_zenith, solar_azimuth=solar_azimuth,
                                          surfaces_configurations=surfaces_configurations)

"""for config in surfaces_configurations:
    if "surface_azimuth" in config:
        if np.all(AOI_for_config[config["name"]] == AOI) is False:
            print(f'Found it: {config["name"]}')
        else:
            print("Nope!")"""

IAM_for_config = calculate_IAM_for_config(AOI_for_config=AOI_for_config, n=n, K=K, L=L,
                                          surfaces_configurations=surfaces_configurations)
"""for config in surfaces_configurations:
    if "surface_azimuth" in config:
        if np.all(IAM_for_config[config["name"]] == IAM) is False:
            print(f'Found it: {config["name"]}')
        else:
            print("Nope!")"""
# In[6]:


# Compute total effective irradiance
TMY_global['I(b)'] = Itot['poa_direct'].values
TMY_global['I(d)'] = Itot['poa_sky_diffuse'].values
TMY_global['I(r)'] = Itot['poa_ground_diffuse'].values

TMY_global['Itot0'] = IAM * \
    (TMY_global['I(b)']+TMY_global['I(d)']+TMY_global['I(r)'])
effective_irr = TMY_global['Itot0']*(1-soiling)

TMY_global_for_config = calculate_TMY_global_for_config(TMY_global=TMY_global, Itot_for_config=Itot_for_config,
                                                        IAM_for_config=IAM_for_config,
                                                        surfaces_configurations=surfaces_configurations)
"""for config in surfaces_configurations:
    if "surface_azimuth" in config:
        if (TMY_global_for_config[config["name"]].reset_index(drop=True)).equals(TMY_global.reset_index(drop=True)):
            print("Nope!")
        else:
            print(f'Found it: {config["name"]}')"""

effective_irr_for_config = calculate_effactive_irr_for_config(TMY_global_for_config=TMY_global_for_config,
                                                              soiling=soiling,
                                                              surfaces_configurations=surfaces_configurations)

"""for config in surfaces_configurations:
    if "surface_azimuth" in config:
        if (effective_irr_for_config[config["name"]].reset_index(drop=True)).equals(effective_irr.reset_index(drop=True)) is False:
            print(f'Found it: {config["name"]}')
        else:
            print("Nope!")"""

# # Compute Tcell

# In[7]:


# Compute effect of temperature (cell temperature °C)
TMY_global['WS5m'] = TMY_global['WS10m'].values*(0.5**(0.23))
WS5m = TMY_global['WS5m'].values
Ta = TMY_global['T2m'].values
TMY_global['Tc'] = effective_irr*np.exp(a+b*WS5m)+Ta+effective_irr/1000*deltaT
Tcell = TMY_global['Tc'].values


WS5m_for_config, Ta_for_config, Tcell_for_donfig = calculate_WS5m_and_Ta_and_Tcell_for_config(TMY_global_for_config=TMY_global_for_config,
                                                                                              effective_irr_for_config=effective_irr_for_config,
                                                                                              a=a, b=b, deltaT=deltaT, surfaces_configurations=surfaces_configurations)

"""for config in surfaces_configurations:
    if "surface_azimuth" in config:
        if np.all(WS5m_for_config[config["name"]] == WS5m) is False:
            print(f'Found it: {config["name"]}')
        else:
            print("Nope!")"""

# Import module characteristics
modules = pvlib.pvsystem.retrieve_sam('cecmod')
module_parameters = modules['Canadian_Solar_Inc__CS1U_410MS']
# Default module size
module_area = 2.078*0.992  # m2
eta_module = module_parameters['I_mp_ref'] * \
    module_parameters['V_mp_ref']/(module_area*1000)


# # Solve single diode equation

# In[9]:


# Calcul des 5 paramètres du modèle du 'single diode model' de Desoto
IL, I0, Rs, Rsh, nNsVth = pvlib.pvsystem.calcparams_desoto(effective_irradiance=effective_irr, temp_cell=Tcell, alpha_sc=module_parameters['alpha_sc'],
                                                           a_ref=module_parameters['a_ref'], I_L_ref=module_parameters['I_L_ref'],
                                                           I_o_ref=module_parameters[
                                                               'I_o_ref'], R_sh_ref=module_parameters['R_sh_ref'],
                                                           R_s=module_parameters['R_s'],
                                                           EgRef=1.121, dEgdT=- 0.0002677, irrad_ref=1000, temp_ref=25)


IL_for_config, I0_for_config, Rs_for_config, Rsh_for_config, nNsVth_for_config = calculate_IL_and_IO_and_Rs_and_Rsh_and_nNsVth_for_config(
    effective_irr_for_config=effective_irr_for_config, Tcell_for_config=Tcell_for_donfig, module_parameters=module_parameters,
    surfaces_configurations=surfaces_configurations
)

"""for config in surfaces_configurations:
    if "surface_azimuth" in config:
        if (IL_for_config[config["name"]].reset_index(drop=True)).equals(IL.reset_index(drop=True)) is False:
            print(f'Found it: {config["name"]}')
        else:
            print("Nope!")

for config in surfaces_configurations:
    if "surface_azimuth" in config:
        if (I0_for_config[config["name"]].reset_index(drop=True)).equals(I0.reset_index(drop=True)) is False:
            print(f'Found it: {config["name"]}')
        else:
            print("Nope!")

for config in surfaces_configurations:
    if "surface_azimuth" in config:
        if (Rs_for_config[config["name"]].reset_index(drop=True)).equals(Rs.reset_index(drop=True)) is False:
            print(f'Found it: {config["name"]}')
        else:
            print("Nope!")

for config in surfaces_configurations:
    if "surface_azimuth" in config:
        if (Rsh_for_config[config["name"]].reset_index(drop=True)).equals(Rsh.reset_index(drop=True)) is False:
            print(f'Found it: {config["name"]}')
        else:
            print("Nope!")

for config in surfaces_configurations:
    if "surface_azimuth" in config:
        if (nNsVth_for_config[config["name"]].reset_index(drop=True)).equals(nNsVth.reset_index(drop=True)) is False:
            print(f'Found it: {config["name"]}')
        else:
            print("Nope!")"""
# In[10]:


# Resolution du modèle
single_res = pvlib.pvsystem.singlediode(photocurrent=IL,
                                        saturation_current=I0,
                                        resistance_series=Rs,
                                        resistance_shunt=Rsh,
                                        nNsVth=nNsVth,
                                        ivcurve_pnts=None,
                                        method='lambertw')


single_res_for_config = calculate_single_res_for_config(IL_for_config=IL_for_config, I0_for_config=I0_for_config,
                                                        Rs_for_config=Rs_for_config, Rsh_for_config=Rsh_for_config,
                                                        nNsVth_for_config=nNsVth_for_config,
                                                        surfaces_configurations=surfaces_configurations)


"""for config in surfaces_configurations:
    if "surface_azimuth" in config:
        if (single_res_for_config[config["name"]].reset_index(drop=True)).equals(single_res.reset_index(drop=True)) is False:
            print(f'Found it: {config["name"]}')
        else:
            print("Nope!")"""

# Calcul de la sortie AC selon la fonction de l'onduleur dévelopée par EnerVivo
P_dc = single_res['p_mp']*(1-mismatch)*(1-connections)*(1-DCwiring)/module_area
P_dc[P_dc <= 0] = 0
inverter_load = P_dc/(eta_module*1000)/r_DCAC/(1-inverterloss)


P_dc_for_config = calculate_P_dc_for_config(single_res_for_config=single_res_for_config,
                                            mismatch=mismatch, connections=connections, DCwiring=DCwiring,
                                            module_area=module_area, surfaces_configurations=surfaces_configurations)

"""for config in surfaces_configurations:
    if "surface_azimuth" in config:
        if (P_dc_for_config[config["name"]].reset_index(drop=True)).equals(P_dc.reset_index(drop=True)) is False:
            print(f'Found it: {config["name"]}')
        else:
            print("Nope!")"""

inverter_load_for_config = calculate_inverter_load_for_config(P_dc_for_config=P_dc_for_config, eta_module=eta_module,
                                                              r_DCAC=r_DCAC, inverterloss=inverterloss,
                                                              surfaces_configurations=surfaces_configurations)

"""for config in surfaces_configurations:
    if "surface_azimuth" in config:
        if (inverter_load_for_config[config["name"]].reset_index(drop=True)).equals(inverter_load.reset_index(drop=True)) is False:
            print(f'Found it: {config["name"]}')
        else:
            print("Nope!")"""

P_ac = P_dc*(1-inverterloss)*(1-np.exp((-inverter_load-0.04)/0.04))
P_ac[P_ac < inverter_min*eta_module*1000/r_DCAC*(1-ACwiring)] = 0
P_ac[P_ac > eta_module*1000/r_DCAC] = eta_module*1000/r_DCAC

P_ac_for_config = calculate_P_ac_for_config(P_dc_for_config=P_dc_for_config, inverterloss=inverterloss,
                                            inverter_load_for_config=inverter_load_for_config,
                                            inverter_min=inverter_min, r_DCAC=r_DCAC,
                                            ACwiring=ACwiring,
                                            surfaces_configurations=surfaces_configurations)

"""for config in surfaces_configurations:
    if "surface_azimuth" in config:
        if (P_ac_for_config[config["name"]].reset_index(drop=True)).equals(P_ac.reset_index(drop=True)) is False:
            print(f'Found it: {config["name"]}')
        else:
            print("Nope!")"""
# In[12]:


# Calcul du nombre d'heure
n_h90 = sum(P_ac)/1000/eta_module*r_P90*(module_eff/eta_module)  # en kWh/kWc

# # Calcul du gain bifacial

# In[13]:


def gain_bifac(bifac, bifac_ratio, H, inter):
    if bifac == 0:
        return 0
    else:
        return rho*bifac*(1.037*(1-1/(np.sqrt(inter)))*(1-np.exp(-(8.691-H)/inter))+0.125*(1-1/inter**4))


gain = gain_bifac(bifac, bifac_ratio, H, inter)


# # Calcul du nombre d'heures (résultat final)

# In[14]:


# résultat final - nombre d'heures équivalent d'ensoleillement, en kWh/kWc/an
n_h90 = n_h90*(1+gain)


P_ac_list = P_ac.tolist()


somme = 0
for elem in P_ac_list:
    somme += elem

print(sum(P_ac_list))
