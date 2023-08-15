import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import yaml

####-------------------------------------------------------------------------------------------------

# Load parameters from YAML file
with open('config/parameters.yml', 'r') as file:
    param = yaml.safe_load(file)

# Load experimental data
df_exp_1 = pd.read_csv('data/batch_no1/data_combined.csv')
df_exp_2 = pd.read_csv('data/batch_no2/data_combined.csv')

F_glu_1 = df_exp_1['Glucose feed [L/min]']*60 #[L/h]
F_in_1 = df_exp_1['Feed total [L/min]']*60 #[L/h]

F_glu_2 = df_exp_2['Glucose feed [L/min]']*60
F_in_2 = df_exp_2['Feed total [L/min]']*60

####-------------------------------------------------------------------------------------------------
# BATCH NO. 1 - qs calculation

def model_single_timestep_qs(i, qs, c_biomass, c_glucose, vol):
    """
    Args:
        i (int): Integer that contains the time step.
        qs (float): Calculated value from ML model for the glucose uptake rate at time step i.
        c_biomass (float): Value for biomass concentration at time step i.
        c_glucose (float): Value for glucose concentration at time step i.
        vol (float): Value for volume at time step i.
    Returns:
        biomass (float): Value for biomass concentration at time step i+1.
        substrate (float): Value for glucose concentration at time step i+1.
        volume (float): Value for volume at time step i+1.
    """

    # Simulation settings
    t0 = 0
    t_end = 45.8
    delta_t = 2
    dt = delta_t/60
    num_steps = int((t_end - t0) / dt) + 1 # Number of time steps

    c_glu_feed = param['c_glu_feed']

    Yxs = param['set_parameter'][0]
    m_s = param['set_parameter'][3]

    # time steps need to be adapted for experimental data input
    t = i * delta_t
    f_glucose = F_glu_1[t]
    f_total = F_in_1[t]
    
    # since the glucose concentration can't be negative, it is set to zero
    if c_glucose < 0:
        c_glucose = 0

    # Update growth and glucose uptake rate
    mu = qs * Yxs

    # Update biomass and substrate concentrations
    dV_dt = f_total - (0.4*60/num_steps) # [L/h] not complete -- include samples + evaporation
    dX_dt = mu * c_biomass - (c_biomass / vol) * dV_dt # [gx/(Lh)]
    dS_dt = ((f_glucose / vol) * (c_glu_feed - c_glucose)) - ((qs + m_s) * c_biomass) - ((c_glucose / vol) * dV_dt)
    # [gs/(Lh)]

    biomass = c_biomass + dX_dt * dt # [gx/L]
    substrate = c_glucose + dS_dt * dt # [gs/L]
    volume = vol + dV_dt * dt # [L]

    return biomass, substrate, volume

####-------------------------------------------------------------------------------------------------
# BATCH NO. 2 - S calculation

def model_single_timestep_S(i, c_glucose, c_biomass, c_co2, vol):
    """
    Args:
        i (int): Integer that contains the time step.
        c_glucose (float): Calculated value from ML model for glucose concentration at time step i.
        c_biomass (float): Value for biomass concentration at time step i.
        c_co2 (float): Value for co2 concentration at time step i.
        vol (float): Value for volume at time step i.
    Returns:
        biomass (float): Value for biomass concentration at time step i+1.
        substrate (float): Value for glucose concentration at time step i+1.
        co2 (float): Value for glucose concentration at time step i+1.
        volume (float): Value for volume at time step i+1.
    """
    
    # Simulation settings
    t0 = 0
    t_end = 36.8
    delta_t = 2
    dt = delta_t/60
    num_steps = int((t_end - t0) / dt) + 1 # Number of time steps

    set1 = param['set_part1']
    set2 = param['set_part2']

    # time steps need to be adapted for experimental data input
    t = i * delta_t
    f_glucose = F_glu_2[t]
    f_total = F_in_2[t]

    if f_glucose < 0.001:
        p = set1
    else:
        p = set2
    
    # Update growth and glucose uptake rate
    qs = p[2] * c_glucose / (p[3] + c_glucose) * (1 / (np.exp(c_biomass * p[5])))
    mu = qs * p[0]

    # Update biomass and substrate concentrations
    dV_dt = f_total - (0.4*60/num_steps) # [L/h] not complete -- include samples + evaporation
    dX_dt = mu * c_biomass - (c_biomass / vol) * dV_dt # [gx/(Lh)]
    dCO2_dt = p[1] * qs * c_biomass  - (c_co2 * (dV_dt / vol)) # [g/(Lh)]

    biomass = c_biomass + dX_dt * dt # [gx/L]
    substrate = c_glucose
    co2 = c_co2 + dCO2_dt * dt  # [g/L]
    volume = vol + dV_dt * dt # [L]

    return biomass, substrate, co2, volume