from typing import List, Dict
import pandas as pd
import pvlib
import numpy as np


modules = pvlib.pvsystem.retrieve_sam('cecmod')
module_parameters = modules['Canadian_Solar_Inc__CS1U_410MS']
module_area = 2.078*0.992
eta_module = module_parameters['I_mp_ref'] * \
    module_parameters['V_mp_ref']/(module_area*1000)


def calculate_Itot_for_config(solar_zenith: float, solar_azimuth: float,
                              dni: float, ghi: float, dhi: float,
                              dni_extra: float, AM: float, rho: float,
                              surfaces_configurations: List[Dict]) -> Dict:

    try:
        result = {config["name"]: pd.DataFrame(pvlib.irradiance.get_total_irradiance(config["surface_tilt"],
                                                                                     config["surface_azimuth"],
                                                                                     solar_zenith, solar_azimuth,
                                                                                     dni, ghi, dhi, dni_extra, AM,
                                                                                     albedo=rho, surface_type=None, model='perez',
                                                                                     model_perez='allsitescomposite1990'))
                  for config in surfaces_configurations
                  if "surface_azimuth" in config}

    except KeyError as e:
        print(f"Missing key: {e}")
        return None

    except Exception as e:
        print(f"An error occurred: {e}")
        return None

    return (result)


def calculate_AOI_for_config(solar_zenith: float, solar_azimuth: float,
                             surfaces_configurations: List[Dict]) -> Dict:
    try:
        result = {config["name"]: (pvlib.irradiance.aoi(
            config["surface_tilt"], config["surface_azimuth"],
            solar_zenith, solar_azimuth))
            for config in surfaces_configurations
            if "surface_azimuth" in config}

    except KeyError as e:
        print(f"Missing key: {e}")
        return None

    except Exception as e:
        print(f"An error occured: {e}")
        return None

    return result


def calculate_IAM_for_config(AOI_for_config: Dict, n: float, K: float,
                             L: float, surfaces_configurations: List[Dict]) -> Dict:
    try:
        IAM_for_config = {config["name"]: (pvlib.iam.physical(
            AOI_for_config[config["name"]], n, K, L
        ))
            for config in surfaces_configurations
            if "surface_azimuth" in config}
    except KeyError as e:
        print(f"Missing key: {e}")
        return None
    except Exception as e:
        print(f'An error occured: {e}')

    return IAM_for_config


def calculate_TMY_global_for_config(TMY_global: pd.DataFrame,
                                    Itot_for_config: Dict,
                                    IAM_for_config: Dict,
                                    surfaces_configurations: List[Dict]) -> Dict:
    TMY_global_for_config = {}
    try:
        for config in surfaces_configurations:
            config_name = config["name"]
            if "surface_azimuth" in config:
                TMY_global_for_config[config["name"]] = TMY_global
                if config_name in Itot_for_config and config_name in TMY_global_for_config:
                    TMY_global_for_config[config_name]['I(b)'] = Itot_for_config[config_name]['poa_direct'].values
                    TMY_global_for_config[config_name]['I(d)'] = Itot_for_config[config_name][
                        'poa_sky_diffuse'].values
                    TMY_global_for_config[config_name]['I(r)'] = Itot_for_config[config_name][
                        'poa_ground_diffuse'].values
                    TMY_global_for_config[config_name]['Itot0'] = IAM_for_config[config_name]*(TMY_global_for_config[config_name]['I(b)']
                                                                                               + TMY_global_for_config[config_name]['I(d)']
                                                                                               + TMY_global_for_config[config_name]['I(r)'])

    except KeyError as e:
        print(f"Missing key: {e}")
        return None

    except Exception as e:
        print(f'An error occured: {e}')

    return TMY_global_for_config


def calculate_effactive_irr_for_config(TMY_global_for_config: Dict,
                                       soiling: float,
                                       surfaces_configurations: List[Dict]) -> Dict:
    try:
        effective_irr_for_config = {config["name"]: (TMY_global_for_config[config["name"]]['Itot0']*(1-soiling))
                                    for config in surfaces_configurations
                                    if config["name"] in TMY_global_for_config}
    except KeyError as e:
        print(f"Missing key: {e}")
        return None
    except Exception as e:
        print(f"An error occured: {e}")

    return effective_irr_for_config


def calculate_WS5m_and_Ta_and_Tcell_for_config(TMY_global_for_config: Dict,
                                               effective_irr_for_config: Dict,
                                               a: float, b: float, deltaT: float,
                                               surfaces_configurations: List[Dict]):
    WS5m_for_config = {}
    Ta_for_config = {}
    Tcell_for_config = {}
    try:
        for config in surfaces_configurations:
            config_name = config["name"]
            if config_name in TMY_global_for_config and config_name in effective_irr_for_config:
                Ta_for_config[config_name] = TMY_global_for_config[config_name]['T2m'].values
                TMY_global_for_config[config_name]["WS5m"] = TMY_global_for_config[config_name]["WS10m"].values*(
                    0.5**(0.23))
                WS5m_for_config[config_name] = TMY_global_for_config[config_name]["WS5m"].values
                Ta_for_config[config_name] = TMY_global_for_config[config_name]['T2m'].values
                TMY_global_for_config[config_name]['Tc'] = effective_irr_for_config[config_name] *\
                    np.exp(a+b*WS5m_for_config[config_name])+Ta_for_config[config_name] +\
                    effective_irr_for_config[config_name]/1000*deltaT
                Tcell_for_config[config_name] = TMY_global_for_config[config_name]['Tc'].values

    except KeyError as e:
        print(f"Missing_key: {e}")
        return None
    except Exception as e:
        print(f"An error accured: {e}")
        return None

    return (WS5m_for_config, Ta_for_config, Tcell_for_config)


def calculate_IL_and_IO_and_Rs_and_Rsh_and_nNsVth_for_config(effective_irr_for_config: Dict,
                                                             Tcell_for_config: Dict,
                                                             module_parameters,
                                                             surfaces_configurations: List[Dict]):
    IL_for_config = {}
    I0_for_config = {}
    Rs_for_config = {}
    Rsh_for_config = {}
    nNsVth_for_config = {}
    try:
        for config in surfaces_configurations:
            config_name = config["name"]
            if "surface_azimuth" in config:
                if effective_irr_for_config[config_name] is not None:
                    IL_for_config[config_name], I0_for_config[config_name], Rs_for_config[config_name], Rsh_for_config[config_name], nNsVth_for_config[config_name] = pvlib.pvsystem.calcparams_desoto(
                        effective_irradiance=effective_irr_for_config[
                            config_name], temp_cell=Tcell_for_config[config_name], alpha_sc=module_parameters['alpha_sc'],
                        a_ref=module_parameters['a_ref'], I_L_ref=module_parameters['I_L_ref'],
                        I_o_ref=module_parameters[
                            'I_o_ref'], R_sh_ref=module_parameters['R_sh_ref'],
                        R_s=module_parameters['R_s'],
                        EgRef=1.121, dEgdT=- 0.0002677, irrad_ref=1000, temp_ref=25)

    except KeyError as e:
        print(f"Messing key: {e}")
        return None
    except Exception as e:
        print(f"An error accured: {e}")
        return None
    return (IL_for_config, I0_for_config, Rs_for_config, Rsh_for_config, nNsVth_for_config)


def calculate_single_res_for_config(IL_for_config, I0_for_config,
                                    Rs_for_config, Rsh_for_config,
                                    nNsVth_for_config,
                                    surfaces_configurations: List[Dict]):
    try:
        single_res_for_config = {config["name"]: pvlib.pvsystem.singlediode(photocurrent=IL_for_config[config["name"]],
                                                                            saturation_current=I0_for_config[
            config["name"]],
            resistance_series=Rs_for_config[config["name"]],
            resistance_shunt=Rsh_for_config[config["name"]],
            nNsVth=nNsVth_for_config[config["name"]],
            ivcurve_pnts=None,
            method='lambertw')
            for config in surfaces_configurations
            if "surface_azimuth" in config and IL_for_config is not None and I0_for_config is not None and Rs_for_config
            is not None and Rsh_for_config is not None and nNsVth_for_config is not None}
    except KeyError as e:
        print(f"Missing key: {e}")
        return None
    except Exception as e:
        print(f"Error accured: {e}")
        return None
    return single_res_for_config


def calculate_P_dc_for_config(single_res_for_config: Dict[str, any],
                              mismatch: float,
                              connections: float,
                              DCwiring: float,
                              module_area: float,
                              surfaces_configurations: List[Dict]) -> Dict[str, any]:
    try:
        P_dc_for_config = {
            config["name"]: single_res_for_config[config["name"]
                                                  ]['p_mp'] * (1 - mismatch) * (1-connections) * (1 - DCwiring) / module_area
            for config in surfaces_configurations
            if "surface_azimuth" in config
        }

        for config in surfaces_configurations:
            if "surface_azimuth" in config:
                series = P_dc_for_config[config["name"]]
                if isinstance(series, (pd.Series, np.ndarray, List)):
                    series[series <= 0] = 0
                else:
                    P_dc_for_config[config["name"]] = max(
                        0, series)

    except KeyError as e:
        print(f"Missing Key: {e}")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None
    return P_dc_for_config


def calculate_inverter_load_for_config(P_dc_for_config: Dict,
                                       eta_module: float,
                                       r_DCAC: float,
                                       inverterloss: float,
                                       surfaces_configurations: List[Dict]):
    try:
        inverter_load_for_config = {config["name"]: P_dc_for_config[config["name"]]/(eta_module*1000)/r_DCAC/(1-inverterloss)
                                    for config in surfaces_configurations
                                    if "surface_azimuth" in config}

    except KeyError as e:
        print(f"Missing Key: {e}")
        return None
    except Exception as e:
        print(f"Error accured: {e}")
        return None
    return inverter_load_for_config


def calculate_P_ac_for_config(P_dc_for_config: Dict,
                              inverterloss: float,
                              inverter_load_for_config: Dict,
                              inverter_min: float,
                              r_DCAC: float,
                              ACwiring: float,
                              surfaces_configurations: List[Dict]):
    try:
        P_ac_for_config = {config["name"]: P_dc_for_config[config["name"]]*(1-inverterloss)*(1-np.exp((-inverter_load_for_config[config["name"]]-0.04)/0.04))
                           for config in surfaces_configurations
                           if "surface_azimuth" in config}
        for config in surfaces_configurations:
            if "surface_azimuth" in config:
                config_name = config["name"]
                P_ac_for_config[config_name][P_ac_for_config[config_name] < inverter_min *
                                             eta_module*1000/r_DCAC*(1-ACwiring)] = 0
                P_ac_for_config[config_name][P_ac_for_config[config_name]
                                             > eta_module*1000/r_DCAC] = eta_module*1000/r_DCAC

    except KeyError as e:
        print(f"Missing key: {e}")
        return None
    except Exception as e:
        print(f"Error accured: {e}")
        return None
    return P_ac_for_config


def gain_bifac(bifac, bifac_ratio, H, inter, rho):
    if bifac == 0:
        return 0
    else:
        return rho*bifac*(1.037*(1-1/(np.sqrt(inter)))*(1-np.exp(-(8.691-H)/inter))+0.125*(1-1/inter**4))


def calculate_n_h90_for_config(P_ac_for_config: Dict,
                               eta_module: float,
                               r_P90: float,
                               module_eff: float,
                               bifac: float,
                               bifac_ratio: float,
                               H: float,
                               inter: float,
                               rho: float,
                               surfaces_configurations: List[Dict]) -> Dict:
    try:
        n_h90_for_config = {config["name"]: sum(P_ac_for_config[config["name"]])/1000/eta_module *
                            r_P90*(module_eff/eta_module)
                            for config in surfaces_configurations
                            if "surface_azimuth" in config}
        gain = gain_bifac(bifac=bifac, bifac_ratio=bifac_ratio, H=H,
                          inter=inter, rho=rho)
        n_h90_for_config = {config["name"]: n_h90_for_config[config["name"]]*(1+gain)
                            for config in surfaces_configurations
                            if "surface_azimuth" in config}
    except KeyError as e:
        print(f"Missing Key: {e}")
        return None
    except Exception as e:
        print(f"Error accured: {e}")
        return None
    return n_h90_for_config
