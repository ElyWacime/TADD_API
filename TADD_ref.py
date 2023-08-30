from pvlib import location
import pandas as pd
import numpy as np
import pvlib
lat = 44.907  # Coordonées du projet
long = -0.3576
# La librairie PVLib possède une fonction qui pemet de récupérer les données PVGIS via son API
TMY_data = pvlib.iotools.get_pvgis_tmy(lat, long, outputformat='csv', usehorizon=True,
                                       userhorizon=None, startyear=2005, endyear=2020, url='https://re.jrc.ec.europa.eu/api/v5_2/',
                                       map_variables=False, timeout=30)
TMY_data = TMY_data[0]
# Remove time zone to save data in excel
TMY_data.set_index(TMY_data.index.tz_localize(None), inplace=True, drop=True)
TMY_global = TMY_data  # Les données météo sont stockées dans la dataframe TMY_global

# Surface Azimuth (deg)
surface_azimuth = 90

# Surface tilt (deg)
surface_tilt = 30

# Module efficiency
module_eff = 0.20  # Default 20%

# Temperature coefficients
a = -3.56  # DEFAUT = -3.56
b = -0.075  # DEFAUT = -0.075
deltaT = 3  # DEFAUT = 3

# Module temperature (@NOCT, %W/°C)
gamma_mp = -0.45  # DEFAUT = -0.45

# Material proprieties
n = 1.526  # DEFAUT = 1.526
L = 0.002  # DEFAUT = 0.002
K = 4  # DEFAUT = 4

# Ground albedo
rho = 0.2  # DEFAUT = 0.2

# Soiling losses
soiling = 0.03  # DEFAUT = 0.03

# DC losses
mismatch = 0.02  # DEFAUT = 0.02
connections = 0.005  # DEFAUT = 0.005
DCwiring = 0.015  # DEFAUT = 0.015
inverterloss = 0.015  # DEFAUT = 0.015
inverter_min = 0.01  # DEFAUT = 0.01

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


# # Get effective irradiance

# In[5]:


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

AOI = pvlib.irradiance.aoi(
    surface_tilt, surface_azimuth, solar_zenith, solar_azimuth)
IAM = pvlib.iam.physical(AOI, n, K, L)
TMY_global['IAM'] = IAM


# In[6]:


# Compute total effective irradiance
TMY_global['I(b)'] = Itot['poa_direct'].values
TMY_global['I(d)'] = Itot['poa_sky_diffuse'].values
TMY_global['I(r)'] = Itot['poa_ground_diffuse'].values


TMY_global['Itot0'] = IAM * \
    (TMY_global['I(b)']+TMY_global['I(d)']+TMY_global['I(r)'])
effective_irr = TMY_global['Itot0']*(1-soiling)


# # Compute Tcell

# In[7]:


# Compute effect of temperature (cell temperature °C)
TMY_global['WS5m'] = TMY_global['WS10m'].values*(0.5**(0.23))
WS5m = TMY_global['WS5m'].values
Ta = TMY_global['T2m'].values
TMY_global['Tc'] = effective_irr*np.exp(a+b*WS5m)+Ta+effective_irr/1000*deltaT
Tcell = TMY_global['Tc'].values


# # Import module characteristics

# In[8]:


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


# In[10]:


# Resolution du modèle
single_res = pvlib.pvsystem.singlediode(photocurrent=IL,
                                        saturation_current=I0,
                                        resistance_series=Rs,
                                        resistance_shunt=Rsh,
                                        nNsVth=nNsVth,
                                        ivcurve_pnts=None,
                                        method='lambertw')


# In[11]:


# Calcul de la sortie AC selon la fonction de l'onduleur dévelopée par EnerVivo
P_dc = single_res['p_mp']*(1-mismatch)*(1-connections)*(1-DCwiring)/module_area
P_dc[P_dc <= 0] = 0
inverter_load = P_dc/(eta_module*1000)/r_DCAC/(1-inverterloss)
P_ac = P_dc*(1-inverterloss)*(1-np.exp((-inverter_load-0.04)/0.04))
P_ac[P_ac < inverter_min*eta_module*1000/r_DCAC*(1-ACwiring)] = 0
P_ac[P_ac > eta_module*1000/r_DCAC] = eta_module*1000/r_DCAC


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


# In[15]:


n_h90


# In[16]:


gain


# In[ ]:

print(sum(P_ac))
