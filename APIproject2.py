from flask import Flask, request, jsonify
import pvlib
import pandas as pd
import numpy as np
from pvlib import location
from typing import Optional
from pydantic import BaseModel, Field
import logging
from logging.handlers import RotatingFileHandler
from datetime import datetime, timedelta

app = Flask(__name__)

if not app.debug:
    file_handler = RotatingFileHandler(
        'flask.log', maxBytes=1024 * 1024 * 100, backupCount=20)
    file_handler.setLevel(logging.ERROR)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(formatter)
    app.logger.addHandler(file_handler)


class Configuration(BaseModel):

    lat: float
    long: float

    toiture_surface1: float = Field(default=0)
    toiture_surface2: float = Field(default=0)
    toiture_surface3: float = Field(default=0)
    toiture_surface4: float = Field(default=0)
    toiture_surface5: float = Field(default=0)
    toiture_surface6: float = Field(default=0)

    serre_surface1: float = Field(default=0)
    serre_surface2: float = Field(default=0)
    serre_surface3: float = Field(default=0)

    ombriere_surface1: float = Field(default=0)
    ombriere_surface2: float = Field(default=0)
    ombriere_surface3: float = Field(default=0)

    module_eff: float = Field(default=0.2)

    # Surface Azimuth (deg)
    surface_azimuth: float = Field(default=0)

    # Surface tilt (deg)
    surface_tilt: float = Field(default=0)

    # Temperature coefficients
    a: float = Field(default=-3.56)  # DEFAUT = -3.56
    b: float = Field(default=-0.075)  # DEFAUT = -0.075
    deltaT: float = Field(default=3.0)  # DEFAUT = 3

    # Module temperature (@NOCT, %W/°C)
    gamma_mp: float = Field(default=-0.45)  # DEFAUT = -0.45

    # Material proprieties
    n: float = Field(default=1.526)  # DEFAUT = 1.526
    L: float = Field(default=0.002)  # DEFAUT = 0.002
    K: float = Field(default=4)  # DEFAUT = 4

    # Ground albedo
    rho: float = Field(default=0.2)  # DEFAUT = 0.2

    # Soiling losses
    soiling: float = Field(default=0.03)  # DEFAUT = 0.03

    # DC losses
    mismatch: float = Field(default=0.02)  # DEFAUT = 0.02
    connections: float = Field(default=0.005)  # DEFAUT = 0.005
    DCwiring: float = Field(default=0.015)  # DEFAUT = 0.02
    inverterloss: float = Field(default=0.015)  # DEFAUT = 0.035
    inverter_min: float = Field(default=0.01)  # DEFAUT = 0.01

    # AC losses
    ACwiring: float = Field(default=0.01)  # DEFAUT = 0.01

    # r_DCAC ratio
    r_DCAC: float = Field(default=1.2)  # DEFAUT = 1.2

    # Inputs pour Ombrière ou Serre / toiture

    """bifac:int=Field(default=0) #Si bifac = 0, module non bifacial, sinon = 1 si bifacial. DEFAUT = 0
    bifac_ratio: float=Field(default=0.65) # DEFAUT = 0.65
    H:float=Field(default=4) # DEFAUT = 4m
    inter:float=Field(default=2.5) # DEFAUT = 2.5m"""
    surface_azimuth_toitureS1: float = Field(default=0)
    surface_azimuth_toitureS2: float = Field(default=0)
    surface_azimuth_toitureS3: float = Field(default=0)
    surface_azimuth_toitureS4: float = Field(default=0)
    surface_azimuth_toitureS5: float = Field(default=0)
    surface_azimuth_toitureS6: float = Field(default=0)

    surface_tilt_toitureS1: float = Field(default=0)
    surface_tilt_toitureS2: float = Field(default=0)
    surface_tilt_toitureS3: float = Field(default=0)
    surface_tilt_toitureS4: float = Field(default=0)
    surface_tilt_toitureS5: float = Field(default=0)
    surface_tilt_toitureS6: float = Field(default=0)

    bifac_toitureS1: float = Field(default=0)
    bifac_toitureS2: float = Field(default=0)
    bifac_toitureS3: float = Field(default=0)
    bifac_toitureS4: float = Field(default=0)
    bifac_toitureS5: float = Field(default=0)
    bifac_toitureS6: float = Field(default=0)

    bifac_ratio_toitureS1: float = Field(default=0.65)
    bifac_ratio_toitureS2: float = Field(default=0.65)
    bifac_ratio_toitureS3: float = Field(default=0.65)
    bifac_ratio_toitureS4: float = Field(default=0.65)
    bifac_ratio_toitureS5: float = Field(default=0.65)
    bifac_ratio_toitureS6: float = Field(default=0.65)

    H_toitureS1: float = Field(default=4)
    H_toitureS2: float = Field(default=4)
    H_toitureS3: float = Field(default=4)
    H_toitureS4: float = Field(default=4)
    H_toitureS5: float = Field(default=4)
    H_toitureS6: float = Field(default=4)

    inter_toitureS1: float = Field(default=2.5)
    inter_toitureS2: float = Field(default=2.5)
    inter_toitureS3: float = Field(default=2.5)
    inter_toitureS4: float = Field(default=2.5)
    inter_toitureS5: float = Field(default=2.5)
    inter_toitureS6: float = Field(default=2.5)

    surface_tilt_ombS1: float = Field(default=0)
    surface_tilt_ombS2: float = Field(default=0)
    surface_tilt_ombS3: float = Field(default=0)

    surface_azimuth_ombS1: float = Field(default=0)
    surface_azimuth_ombS2: float = Field(default=0)
    surface_azimuth_ombS3: float = Field(default=0)

    bifac_ombS1: float = Field(default=0)
    bifac_ombS2: float = Field(default=0)
    bifac_ombS3: float = Field(default=0)

    bifac_ratio_ombS1: float = Field(default=0.65)
    bifac_ratio_ombS2: float = Field(default=0.65)
    bifac_ratio_ombS3: float = Field(default=0.65)

    H_ombS1: float = Field(default=4)
    H_ombS2: float = Field(default=4)
    H_ombS3: float = Field(default=4)

    inter_ombS1: float = Field(default=2.5)
    inter_ombS2: float = Field(default=2.5)
    inter_ombS3: float = Field(default=2.5)

    surface_tilt_serS1: float = Field(default=0)
    surface_tilt_serS2: float = Field(default=0)
    surface_tilt_serS3: float = Field(default=0)

    surface_azimuth_serS1: float = Field(default=0)
    surface_azimuth_serS2: float = Field(default=0)
    surface_azimuth_serS3: float = Field(default=0)

    bifac_serS1: float = Field(default=0)
    bifac_serS2: float = Field(default=0)
    bifac_serS3: float = Field(default=0)

    bifac_ratio_serS1: float = Field(default=0.65)
    bifac_ratio_serS2: float = Field(default=0.65)
    bifac_ratio_serS3: float = Field(default=0.65)

    H_serS1: float = Field(default=4)
    H_serS2: float = Field(default=4)
    H_serS3: float = Field(default=4)

    inter_serS1: float = Field(default=2.5)
    inter_serS2: float = Field(default=2.5)
    inter_serS3: float = Field(default=2.5)

    # Ratio P90/P50
    r_P90: float = Field(default=0.94)


@app.route('/calculate_pv', methods=['POST'])
def calculate_pv():
    configuration_data = request.get_json()
    configuration = Configuration(**configuration_data)

    # La librairie PVLib possède une fonction qui pemet de récupérer les données PVGIS via son API
    TMY_data = pvlib.iotools.get_pvgis_tmy(configuration.lat, configuration.long, outputformat='csv', usehorizon=True,
                                           userhorizon=None, startyear=2005, endyear=2016, url='https://re.jrc.ec.europa.eu/api/v5_2/',
                                           map_variables=False, timeout=30)
    TMY_data = TMY_data[0]
    # Remove time zone to save data in excel
    TMY_data.set_index(TMY_data.index.tz_localize(None),
                       inplace=True, drop=True)
    TMY_global = TMY_data  # Les données météo sont stockées dans la dataframe TMY_global

    surfaces_configurations = [
        # Configurations for toiture
        {"name": "toiture_surface1", "surface": configuration.toiture_surface1,  "surface_azimuth": configuration.surface_azimuth_toitureS1, "surface_tilt": configuration.surface_tilt_toitureS1, "bifac": configuration.bifac_toitureS1,
            "bifac_ratio": configuration.bifac_ratio_toitureS1, "H": configuration.H_toitureS1, "inter": configuration.inter_toitureS1},
        {"name": "toiture_surface2", "surface": configuration.toiture_surface2,  "surface_azimuth": configuration.surface_azimuth_toitureS2, "surface_tilt": configuration.surface_tilt_toitureS2, "bifac": configuration.bifac_toitureS2,
            "bifac_ratio": configuration.bifac_ratio_toitureS2, "H": configuration.H_toitureS2, "inter": configuration.inter_toitureS2},
        {"name": "toiture_surface3", "surface": configuration.toiture_surface3,  "surface_azimuth": configuration.surface_azimuth_toitureS3, "surface_tilt": configuration.surface_tilt_toitureS3, "bifac": configuration.bifac_toitureS3,
            "bifac_ratio": configuration.bifac_ratio_toitureS3, "H": configuration.H_toitureS3, "inter": configuration.inter_toitureS3},
        {"name": "toiture_surface4", "surface": configuration.toiture_surface4,  "surface_azimuth": configuration.surface_azimuth_toitureS4, "surface_tilt": configuration.surface_tilt_toitureS4, "bifac": configuration.bifac_toitureS4,
            "bifac_ratio": configuration.bifac_ratio_toitureS4, "H": configuration.H_toitureS4, "inter": configuration.inter_toitureS4},
        {"name": "toiture_surface5", "surface": configuration.toiture_surface5,  "surface_azimuth": configuration.surface_azimuth_toitureS5, "surface_tilt": configuration.surface_tilt_toitureS5, "bifac": configuration.bifac_toitureS5,
            "bifac_ratio": configuration.bifac_ratio_toitureS5, "H": configuration.H_toitureS5, "inter": configuration.inter_toitureS5},
        {"name": "toiture_surface6", "surface": configuration.toiture_surface6,  "surface_azimuth": configuration.surface_azimuth_toitureS6, "surface_tilt": configuration.surface_tilt_toitureS6, "bifac": configuration.bifac_toitureS6,
            "bifac_ratio": configuration.bifac_ratio_toitureS6, "H": configuration.H_toitureS6, "inter": configuration.inter_toitureS6},
        # Configurations for ombrière
        {"name": "ombriere_surface1", "surface": configuration.ombriere_surface1, "surface_azimuth": configuration.surface_azimuth_ombS1,
            "surface_tilt": configuration.surface_tilt_ombS1, "bifac": configuration.bifac_ombS1, "bifac_ratio": configuration.bifac_ratio_ombS1, "H": configuration.H_ombS1, "inter": configuration.inter_ombS1},
        {"name": "ombriere_surface2", "surface": configuration.ombriere_surface2, "surface_azimuth": configuration.surface_azimuth_ombS2,
            "surface_tilt": configuration.surface_tilt_ombS2, "bifac": configuration.bifac_ombS2, "bifac_ratio": configuration.bifac_ratio_ombS2, "H": configuration.H_ombS2, "inter": configuration.inter_ombS2},
        {"name": "ombriere_surface3", "surface": configuration.ombriere_surface3, "surface_azimuth": configuration.surface_azimuth_ombS3,
            "surface_tilt": configuration.surface_tilt_ombS3, "bifac": configuration.bifac_ombS3, "bifac_ratio": configuration.bifac_ratio_ombS3, "H": configuration.H_ombS3, "inter": configuration.inter_ombS3},
        # Configurations for serre
        {"name": "serre_surface1", "surface": configuration.serre_surface1, "surface_azimuth": configuration.surface_azimuth_serS1,
            "surface_tilt": configuration.surface_tilt_serS1, "bifac": configuration.bifac_serS1, "bifac_ratio": configuration.bifac_ratio_serS1, "H": configuration.H_serS1, "inter": configuration.inter_serS1},
        {"name": "serre_surface2", "surface": configuration.serre_surface2, "surface_azimuth": configuration.surface_azimuth_serS2,
            "surface_tilt": configuration.surface_tilt_serS2, "bifac": configuration.bifac_serS2, "bifac_ratio": configuration.bifac_ratio_serS2, "H": configuration.H_serS2, "inter": configuration.inter_serS2},
        {"name": "serre_surface3", "surface": configuration.serre_surface3, "surface_azimuth": configuration.surface_azimuth_serS3,
            "surface_tilt": configuration.surface_tilt_serS3, "bifac": configuration.bifac_serS3, "bifac_ratio": configuration.bifac_ratio_serS3, "H": configuration.H_serS3, "inter": configuration.inter_serS3},
    ]
    # Import solar Pos:
    times = pd.date_range('2005-01-01', periods=8760,
                          freq='1H', tz='Europe/Paris')
    loc = location.Location(configuration.lat, configuration.long, tz=times.tz)
    sp = loc.get_solarposition(times)
    # Récupération des angles solaires via le module 'loc' de pvlib
    solar_zenith = sp['apparent_zenith'].values
    solar_azimuth = sp['azimuth'].values

    # Get Itot using perez:
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
    Itot = pvlib.irradiance.get_total_irradiance(configuration.surface_tilt, configuration.surface_azimuth, solar_zenith, solar_azimuth,
                                                 dni, ghi, dhi, dni_extra, AM, albedo=configuration.rho, surface_type=None,
                                                 model='perez', model_perez='allsitescomposite1990')
    Itot = pd.DataFrame(Itot)

    ########### ************#############
    def calculate_Itot_for_each_configuration(surfaces_config):
        if "surface_azimuth" in surfaces_config:
            result = pvlib.irradiance.get_total_irradiance(surfaces_config["surface_tilt"], surfaces_config["surface_azimuth"],
                                                           solar_zenith, solar_azimuth, dni, ghi, dhi, dni_extra, AM,
                                                           albedo=configuration.rho, surface_type=None, model='perez',
                                                           model_perez='allsitescomposite1990')
            return pd.DataFrame(result)
    Itot_for_config = {config["name"]: calculate_Itot_for_each_configuration(
        config) for config in surfaces_configurations}
    ############### ***************############

    # get effective irradiation:
    # pd.set_option('display.max_columns', None)
    # pd.set_option('display.max_rows', None)

    AOI = pvlib.irradiance.aoi(
        configuration.surface_tilt, configuration.surface_azimuth, solar_zenith, solar_azimuth)
    IAM = pvlib.iam.physical(AOI, configuration.n,
                             configuration.K, configuration.L)
    TMY_global['IAM'] = IAM

    def calculate_AOI_for_config(config):
        if "surface_azimuth" in config:
            result = pvlib.irradiance.aoi(
                config["surface_tilt"], config["surface_azimuth"], solar_zenith, solar_azimuth)
            return result

    AOI_for_config = {config["name"]: calculate_AOI_for_config(
        config) for config in surfaces_configurations}

    def calculate_IAM_for_config(config):
        if "surface_azimuth" in config:
            result = pvlib.iam.physical(
                AOI_for_config[config["name"]], configuration.n, configuration.K, configuration.L)
            return result

    IAM_for_config = {config["name"]: calculate_IAM_for_config(
        config) for config in surfaces_configurations}

    TMY_global_for_config = {}
    for config in surfaces_configurations:
        if "surface_azimuth" in config:
            TMY_global_for_config[config["name"]] = TMY_global

    for config in surfaces_configurations:
        config_name = config["name"]
        if config_name in IAM_for_config and config_name in TMY_global_for_config:
            TMY_global_for_config[config_name]['IAM'] = IAM_for_config[config_name]
    ############# **************#############

    # Compute total effective irradiance
    TMY_global['I(b)'] = Itot['poa_direct'].values
    TMY_global['I(d)'] = Itot['poa_sky_diffuse'].values
    TMY_global['I(r)'] = Itot['poa_ground_diffuse'].values

    for config in surfaces_configurations:
        config_name = config["name"]
        if config_name in Itot_for_config and config_name in TMY_global_for_config:
            if Itot_for_config[config_name] is not None:
                TMY_global_for_config[config_name]['I(b)'] = Itot_for_config[config_name]['poa_direct'].values
                TMY_global_for_config[config_name]['I(d)'] = Itot_for_config[config_name]['poa_sky_diffuse'].values
                TMY_global_for_config[config_name]['I(r)'] = Itot_for_config[config_name]['poa_ground_diffuse'].values
    ############# ********************#############

    TMY_global['Itot0'] = IAM * \
        (TMY_global['I(b)']+TMY_global['I(d)']+TMY_global['I(r)'])
    effective_irr = TMY_global['Itot0']*(1-configuration.soiling)

    effective_irr_for_config = {}
    for config in surfaces_configurations:
        config_name = config["name"]
        if config_name in TMY_global_for_config and config_name in IAM_for_config:
            if IAM_for_config[config_name] is not None:
                TMY_global_for_config[config_name]['Itot0'] = IAM_for_config[config_name]*(TMY_global_for_config[config_name]['I(b)']
                                                                                           + TMY_global_for_config[config_name]['I(d)']
                                                                                           + TMY_global_for_config[config_name]['I(r)'])
                effective_irr_for_config[config_name] = TMY_global_for_config[config_name]['Itot0']*(
                    1-configuration.soiling)
    ##################### ******************##############

    # compute Tcell:
    # Compute effect of temperature (cell temperature °C)
    TMY_global['WS5m'] = TMY_global['WS10m'].values*(0.5**(0.23))
    WS5m = TMY_global['WS5m'].values
    Ta = TMY_global['T2m'].values
    TMY_global['Tc'] = effective_irr*np.exp(
        configuration.a+configuration.b*WS5m)+Ta+effective_irr/1000*configuration.deltaT
    Tcell = TMY_global['Tc'].values

    WS5m_for_config = {}
    for config in surfaces_configurations:
        config_name = config["name"]
        if config_name in TMY_global_for_config and config_name in effective_irr_for_config:
            TMY_global_for_config[config_name]["WS5m"] = TMY_global_for_config[config_name]["WS10m"].values*(
                0.5**(0.23))
            WS5m_for_config[config_name] = TMY_global_for_config[config_name]["WS5m"]

            TMY_global_for_config[config_name]['Tc'] = (effective_irr_for_config[config_name].values) * \
                np.exp(configuration.a+configuration.b*(WS5m_for_config[config_name].values)+Ta +
                       effective_irr_for_config[config_name]/1000*configuration.deltaT)
    ############ *******************################
    # Import module charactéristics:
    modules = pvlib.pvsystem.retrieve_sam('cecmod')
    module_parameters = modules['Canadian_Solar_Inc__CS1U_410MS']
    # Default module size
    module_area = 2.078*0.992  # m2
    eta_module = module_parameters['I_mp_ref'] * \
        module_parameters['V_mp_ref']/(module_area*1000)

    # Solve single diode equation:
    # Calcul des 5 paramètres du modèle du 'single diode model' de Desoto
    IL, I0, Rs, Rsh, nNsVth = pvlib.pvsystem.calcparams_desoto(effective_irradiance=effective_irr, temp_cell=Tcell, alpha_sc=module_parameters['alpha_sc'],
                                                               a_ref=module_parameters['a_ref'], I_L_ref=module_parameters['I_L_ref'],
                                                               I_o_ref=module_parameters[
                                                                   'I_o_ref'], R_sh_ref=module_parameters['R_sh_ref'],
                                                               R_s=module_parameters['R_s'],
                                                               EgRef=1.121, dEgdT=- 0.0002677, irrad_ref=1000, temp_ref=25)
    IL_for_config = {}
    I0_for_config = {}
    Rs_for_config = {}
    Rsh_for_config = {}
    nNsVth_for_config = {}
    for config in surfaces_configurations:
        config_name = config["name"]
        if "surface_azimuth" in config:
            if effective_irr_for_config[config_name] is not None:
                IL_for_config[config_name], I0_for_config[config_name], Rs_for_config[config_name], Rsh_for_config[config_name], nNsVth_for_config[config_name] = pvlib.pvsystem.calcparams_desoto(
                    effective_irradiance=effective_irr_for_config[
                        config_name], temp_cell=Tcell, alpha_sc=module_parameters['alpha_sc'],
                    a_ref=module_parameters['a_ref'], I_L_ref=module_parameters['I_L_ref'],
                    I_o_ref=module_parameters[
                        'I_o_ref'], R_sh_ref=module_parameters['R_sh_ref'],
                    R_s=module_parameters['R_s'],
                    EgRef=1.121, dEgdT=- 0.0002677, irrad_ref=1000, temp_ref=25)
    ########## ****************#############

    # Resolution du modèle
    single_res = pvlib.pvsystem.singlediode(photocurrent=IL,
                                            saturation_current=I0,
                                            resistance_series=Rs,
                                            resistance_shunt=Rsh,
                                            nNsVth=nNsVth,
                                            ivcurve_pnts=None,
                                            method='lambertw')

    single_res_for_config = {}
    for config in surfaces_configurations:
        config_name = config["name"]
        if "surface_azimuth" in config:
            if IL_for_config is not None and I0_for_config is not None and Rs_for_config is not None and Rsh_for_config is not None and nNsVth_for_config is not None:
                single_res_for_config[config_name] = pvlib.pvsystem.singlediode(photocurrent=IL_for_config[config_name],
                                                                                saturation_current=I0_for_config[
                                                                                    config_name],
                                                                                resistance_series=Rs_for_config[config_name],
                                                                                resistance_shunt=Rsh_for_config[config_name],
                                                                                nNsVth=nNsVth,
                                                                                ivcurve_pnts=None,
                                                                                method='lambertw')
    ############# ****************##############

    # Calcul de la sortie AC selon la fonction de l'onduleur dévelopée par EnerVivo
    P_dc = single_res['p_mp']*(1-configuration.mismatch)*(
        1-configuration.connections)*(1-configuration.DCwiring)/module_area
    P_dc[P_dc <= 0] = 0
    inverter_load = P_dc/(eta_module*1000) / \
        configuration.r_DCAC/(1-configuration.inverterloss)
    P_ac = P_dc*(1-configuration.inverterloss) * \
        (1-np.exp((-inverter_load-0.04)/0.04))
    P_ac[P_ac < configuration.inverter_min*eta_module*1000 /
         configuration.r_DCAC*(1-configuration.ACwiring)] = 0
    P_ac[P_ac > eta_module*1000/configuration.r_DCAC] = eta_module * \
        1000/configuration.r_DCAC

    P_dc_for_config = {}
    P_ac_for_config = {}
    for config in surfaces_configurations:
        config_name = config["name"]
        if "surface_azimuth" in config and single_res_for_config[config_name] is not None:
            P_dc_for_config[config_name] = single_res_for_config[config_name]['p_mp']*(
                1-configuration.mismatch)*(1-configuration.connections)*(1-configuration.DCwiring)/module_area

            P_dc_for_config[config_name].loc[P_dc_for_config[config_name] <= 0] = 0

            P_ac_for_config[config_name] = P_dc_for_config[config_name] * \
                (1-configuration.inverterloss) * \
                (1-np.exp((-inverter_load-0.04)/0.04))
            P_ac_for_config[config_name].loc[P_ac_for_config[config_name]
                                             < configuration.inverter_min*eta_module*1000/configuration.r_DCAC
                                             * (1-configuration.ACwiring)] = 0
            P_ac_for_config[config_name].loc[P_ac_for_config[config_name]
                                             > eta_module*1000/configuration.r_DCAC] \
                = eta_module*1000/configuration.r_DCAC

    ########### **************#############

    # Calcul du nombre d'heure
    n_h90 = sum(P_ac)/1000/eta_module*configuration.r_P90  # en kWh/kWc
    n_h90_for_config = {config["name"]: (sum(
        P_ac_for_config[config["name"]]/1000/eta_module*configuration.r_P90*(configuration.module_eff/eta_module)))
        for config in surfaces_configurations
        if "surface_azimuth" in config}

    ########### ***************###########

    # Calcul du gain bifacial:
    def gain_bifac(bifac, bifac_ratio, H, inter):
        return configuration.rho*bifac*(1.037*(1-1/(np.sqrt(inter)))*(1-np.exp(-(8.691-H)/inter))+0.125*(1-1/inter**4))

    def calculate_gain_for_surface(surface_config):
        if surface_config["surface"] != 0:
            return gain_bifac(
                surface_config["bifac"],
                surface_config["bifac_ratio"],
                surface_config["H"],
                surface_config["inter"]
            )
        return 0

    gains = {config["name"]: calculate_gain_for_surface(config)
             for config in surfaces_configurations}

    # Résulta final "nombre d'heur":
    n_h90  # résultat final - nombre d'heures équivalent d'ensoleillement, en kWh/kWc/an

    product_for_config = {config["name"]: P_ac_for_config[config["name"]] * config["surface"]
                          for config in surfaces_configurations
                          if config["name"] in P_ac_for_config}

    product_for_config = pd.concat(P_ac_for_config, axis=1)
    sum_P_ac_for_config = product_for_config.sum(axis=1)

    sum_P_ac_for_config = sum_P_ac_for_config.tolist()

    # data = {'hourly_P_ac': P_ac_for_config, 'result': result, 'gain': gains}
    data = {"result": n_h90_for_config, "gain": gains,
            "hourly_P_ac": sum_P_ac_for_config}
    return jsonify(n_h90_for_config)


if __name__ == '__main__':
    app.run(debug=True)

"""Dans "PV configuration", imprimer la variable "result", pour chaque configuration: càd dans les cellules F26, F45, F64, F84, F104, F123, F151, F175, F199, F227, F251 et F275"""
