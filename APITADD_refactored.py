from pvlib import location
import pandas as pd
import numpy as np
import pvlib
from flask import Flask, request, jsonify
from pydantic import BaseModel, Field
import logging
from logging.handlers import RotatingFileHandler
from datetime import datetime, timedelta
from refactored_code import calculate_n_h90_for_config
from refactored_code import calculate_P_ac_for_config, calculate_inverter_load_for_config
from refactored_code import calculate_P_dc_for_config
from refactored_code import modules, module_parameters, module_area, eta_module, calculate_single_res_for_config, calculate_Itot_for_config, calculate_IL_and_IO_and_Rs_and_Rsh_and_nNsVth_for_config, calculate_WS5m_and_Ta_and_Tcell_for_config, calculate_effactive_irr_for_config, calculate_IAM_for_config, calculate_AOI_for_config, calculate_TMY_global_for_config


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
    # AC losses
    ACwiring: float = Field(default=0.01)  # DEFAUT = 0.01

    # r_DCAC ratio
    r_DCAC: float = Field(default=1.2)  # DEFAUT = 1.2

    # Inputs pour Ombrière ou Serre
    # Si bifac = 0, module non bifacial, sinon = 1 si bifacial. DEFAUT = 0
    bifac: float = Field(default=0)
    bifac_ratio: float = Field(default=0.65)  # DEFAUT = 0.65
    H: float = Field(default=4)  # DEFAUT = 4m
    inter: float = Field(default=2.5)  # DEFAUT = 2.5m

    lat: float = Field(default=43.5161)  # DEFAUT = 43.5161
    long: float = Field(default=-1.0285)  # DEFAUT = -1.0285
    surface_azimuth: int = Field(default=180)  # DEFAUT = 180
    surface_tilt: int = Field(default=0)  # DEFAUT = 0
    a: float = Field(default=-3.56)  # DEFAUT = -3.56
    b: float = Field(default=-0.075)  # DEFAUT = -0.075
    deltaT: int = Field(default=3)  # DEFAUT = 3
    gamma_mp: float = Field(default=-0.45)  # DEFAUT = -0.45
    n: float = Field(default=1.526)  # DEFAUT = 1.526
    L: float = Field(default=0.002)  # DEFAUT = 0.002
    K: int = Field(default=4)  # DEFAUT = 4
    rho: float = Field(default=0.2)  # DEFAUT = 0.2
    soiling: float = Field(default=0.03)  # DEFAUT = 0.03
    module_eff: float = Field(default=0.20)  # DEFAUT = 0.20
    mismatch: float = Field(default=0.02)  # DEFAUT = 0.02
    connections: float = Field(default=0.005)  # DEFAUT = 0.005
    DCwiring: float = Field(default=0.015)  # DEFAUT = 0.015
    inverterloss: float = Field(default=0.015)  # DEFAUT = 0.015
    inverter_min: float = Field(default=0.01)  # DEFAUT = 0.01
    r_P90: float = Field(default=0.94)  # DEFAUT = 0.94

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

    surface_azimuth_ombS1: float = Field(default=0)
    surface_azimuth_ombS2: float = Field(default=0)
    surface_azimuth_ombS3: float = Field(default=0)

    surface_tilt_ombS1: float = Field(default=0)
    surface_tilt_ombS2: float = Field(default=0)
    surface_tilt_ombS3: float = Field(default=0)

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

    surface_azimuth_serS1: float = Field(default=0)
    surface_azimuth_serS2: float = Field(default=0)
    surface_azimuth_serS3: float = Field(default=0)

    surface_tilt_serS1: float = Field(default=0)
    surface_tilt_serS2: float = Field(default=0)
    surface_tilt_serS3: float = Field(default=0)

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
    TMY_global = TMY_data

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

    times = pd.date_range('2005-01-01', periods=8760,
                          freq='1H', tz='Europe/Paris')
    loc = location.Location(configuration.lat, configuration.long, tz=times.tz)
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
                                                AM=AM, rho=configuration.rho,
                                                surfaces_configurations=surfaces_configurations)

    AOI_for_config = calculate_AOI_for_config(solar_zenith=solar_zenith, solar_azimuth=solar_azimuth,
                                              surfaces_configurations=surfaces_configurations)

    IAM_for_config = calculate_IAM_for_config(AOI_for_config=AOI_for_config,
                                              n=configuration.n, K=configuration.K, L=configuration.L,
                                              surfaces_configurations=surfaces_configurations)

    TMY_global_for_config = calculate_TMY_global_for_config(TMY_global=TMY_global,
                                                            Itot_for_config=Itot_for_config,
                                                            IAM_for_config=IAM_for_config,
                                                            surfaces_configurations=surfaces_configurations)

    effective_irr_for_config = calculate_effactive_irr_for_config(TMY_global_for_config=TMY_global_for_config,
                                                                  soiling=configuration.soiling, surfaces_configurations=surfaces_configurations)

    (WS5m_for_config, Ta_for_config, Tcell_for_config) = calculate_WS5m_and_Ta_and_Tcell_for_config(TMY_global_for_config=TMY_global_for_config,
                                                                                                    effective_irr_for_config=effective_irr_for_config,
                                                                                                    a=configuration.a, b=configuration.b, deltaT=configuration.deltaT,
                                                                                                    surfaces_configurations=surfaces_configurations)

    (IL_for_config, I0_for_config,
     Rs_for_config, Rsh_for_config, nNsVth_for_config) = calculate_IL_and_IO_and_Rs_and_Rsh_and_nNsVth_for_config(effective_irr_for_config=effective_irr_for_config,
                                                                                                                  Tcell_for_config=Tcell_for_config,
                                                                                                                  module_parameters=module_parameters,
                                                                                                                  surfaces_configurations=surfaces_configurations)

    single_res_for_config = calculate_single_res_for_config(IL_for_config=IL_for_config, I0_for_config=I0_for_config,
                                                            Rs_for_config=Rs_for_config, Rsh_for_config=Rsh_for_config,
                                                            nNsVth_for_config=nNsVth_for_config,
                                                            surfaces_configurations=surfaces_configurations)

    P_dc_for_config = calculate_P_dc_for_config(single_res_for_config=single_res_for_config,
                                                mismatch=configuration.mismatch, connections=configuration.connections,
                                                DCwiring=configuration.DCwiring, module_area=module_area,
                                                surfaces_configurations=surfaces_configurations)

    inverter_load_for_config = calculate_inverter_load_for_config(P_dc_for_config=P_dc_for_config,
                                                                  eta_module=eta_module,
                                                                  r_DCAC=configuration.r_DCAC, inverterloss=configuration.inverterloss,
                                                                  surfaces_configurations=surfaces_configurations)

    P_ac_for_config = calculate_P_ac_for_config(inverterloss=configuration.inverterloss, P_dc_for_config=P_dc_for_config,
                                                inverter_load_for_config=inverter_load_for_config,
                                                inverter_min=configuration.inverter_min, r_DCAC=configuration.r_DCAC,
                                                ACwiring=configuration.ACwiring,
                                                surfaces_configurations=surfaces_configurations)

    n_h90_for_config = calculate_n_h90_for_config(P_ac_for_config=P_ac_for_config, r_P90=configuration.r_P90,
                                                  eta_module=eta_module, module_eff=configuration.module_eff,
                                                  surfaces_configurations=surfaces_configurations,
                                                  bifac=configuration.bifac, bifac_ratio=configuration.bifac_ratio, H=configuration.H, inter=configuration.inter,
                                                  rho=configuration.rho)

    product_for_config = pd.concat(P_ac_for_config, axis=1)
    sum_P_ac_for_config = product_for_config.sum(axis=1)
    sum_P_ac_for_config = sum_P_ac_for_config.tolist()

    data = {"P_ac_for_config": sum_P_ac_for_config,
            "n_h90_for_config": n_h90_for_config}

    return jsonify(data)


if __name__ == '__main__':
    app.run(debug=True)
