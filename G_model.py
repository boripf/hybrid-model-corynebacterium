import pandas as pd
import numpy as np
import yaml

# Load experimental data
df_exp = pd.read_csv('data/data_combined.csv')
# we can not use that because the total feed contains acid and base
F_glu = df_exp['Glucose feed [L/min]']*60 #[L/h]
F_in = df_exp['Feed total [L/min]']*60 #[L/h]

# Load parameters from YAML file
with open('config/parameters.yml', 'r') as file:
    param = yaml.safe_load(file)

def model_single_timestep(i, qs, c_biomass, c_glucose, vol):
    """
    Simulates the fermentation process based on the provided parameters.
    Args:
        params (dict): Dictionary containing the model parameters.
        t0 (float): Initial time.
        t_end (float): End time.
        dt (float): Time step size.
    Returns:
        time (array): Array of time values.
        biomass (array): Array of biomass concentrations.
        substrate (array): Array of substrate concentrations.
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
    f_glucose = F_glu[t]
    f_total = F_in[t]
    
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