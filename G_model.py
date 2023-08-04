import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import yaml

# Load experimental data
df_exp = pd.read_csv('data/batch_no1/data_combined.csv')
# we can not use that because the total feed contains acid and base
F_glu = df_exp['Glucose feed [L/min]']*60 #[L/h]
F_in = df_exp['Feed total [L/min]']*60 #[L/h]

# Load parameters from YAML file
with open('config/parameters.yml', 'r') as file:
    param = yaml.safe_load(file)

def model_single_timestep_qs(i, qs, c_biomass, c_glucose, vol):
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

def model_single_timestep_S(i, c_glucose, c_biomass, c_co2, vol):
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
    t_end = 36.8
    delta_t = 2
    dt = delta_t/60
    num_steps = int((t_end - t0) / dt) + 1 # Number of time steps

    Yxs = param['set_6parameter_no2'][0]
    Yco2s = param['set_6parameter_no2'][1]
    qs_max = param['set_6parameter_no2'][2]
    Ks = param['set_6parameter_no2'][3]
    lag = param['set_6parameter_no2'][5]

    # time steps need to be adapted for experimental data input
    t = i * delta_t
    f_total = F_in[t]
    
    # Update growth and glucose uptake rate
    qs = qs_max * c_glucose / (Ks + c_glucose) * (1 / (np.exp(c_biomass * lag)))
    mu = qs * Yxs

    # Update biomass and substrate concentrations
    dV_dt = f_total - (0.4*60/num_steps) # [L/h] not complete -- include samples + evaporation
    dX_dt = mu * c_biomass - (c_biomass / vol) * dV_dt # [gx/(Lh)]
    dCO2_dt = Yco2s * qs * c_biomass  - (c_co2 * (dV_dt / vol)) # [g/(Lh)]

    biomass = c_biomass + dX_dt * dt # [gx/L]
    substrate = c_glucose
    co2 = c_co2 + dCO2_dt * dt  # [g/L]
    volume = vol + dV_dt * dt # [L]

    return biomass, substrate, co2, volume

####-------------------------------------------------------------------------------------------------

# Load experimental data
df_exp_2 = pd.read_csv('data/batch_no2/data_combined.csv')
# we can not use that because the total feed contains acid and base
F_glu = df_exp_2['Glucose feed [L/min]']*60 #[L/h]
F_in = df_exp_2['Feed total [L/min]']*60 #[L/h]

def model_batch_no2(parameters):
    # Simulation settings
    t0 = 0
    t_end = 36.8
    delta_t = 2
    dt = delta_t/60
    num_steps = int((t_end - t0) / dt) + 1 # Number of time steps

    # Extract parameters
    X0 = param['X0_2']
    S0 = param['S0']
    co20 = param['co20']
    V0 = param['V0_2']
    c_glu_feed = param['c_glu_feed']

    Yxs = parameters[0]
    Yco2s = parameters[1]
    qs_max = parameters[2]
    Ks = parameters[3]
    m_s = parameters[4]
    lag = parameters[5]

    # Arrays to store results
    time = np.linspace(t0, t_end, num_steps)
    biomass = np.zeros(num_steps)
    substrate = np.zeros(num_steps)
    co2 = np.zeros(num_steps)
    volume = np.zeros(num_steps)

    # Set initial values
    biomass[0] = X0
    substrate[0] = S0
    co2[0] = co20
    volume[0] = V0

    # Simulation loop
    for i in range(1, num_steps):
        c_glucose = substrate[i-1]
        c_biomass = biomass[i-1]
        c_co2 = co2[i-1]
        vol = volume[i-1]

        # time steps need to be adapted for experimental data input
        t = i * delta_t
        f_glucose = F_glu[t]
        f_total = F_in[t]
        
        # since the glucose concentration can't be negative, it is set to zero
        if c_glucose < 0:
            c_glucose = 0

        # Update growth and glucose uptake rate
        qs = qs_max * c_glucose / (Ks + c_glucose) * (1 / (np.exp(c_biomass * lag)))
        mu = qs * Yxs

        # Update biomass and substrate concentrations
        dV_dt = f_total - (0.4*60/num_steps) # [L/h] not complete -- include samples + evaporation
        dX_dt = mu * c_biomass - (c_biomass / vol) * dV_dt # [gx/(Lh)]
        dS_dt = ((f_glucose / vol) * (c_glu_feed - c_glucose)) - ((qs + m_s) * c_biomass) - ((c_glucose / vol) * dV_dt)
        # [gs/(Lh)]
        dCO2_dt = Yco2s * qs * c_biomass  - (c_co2 * (dV_dt / vol)) # [g/(Lh)]

        biomass[i] = c_biomass + dX_dt * dt # [gx/L]
        substrate[i] = c_glucose + dS_dt * dt # [gs/L]
        co2[i] = c_co2 + dCO2_dt * dt  # [g/L]
        volume[i] = vol + dV_dt * dt # [L]

    return time, biomass, substrate, co2

def plot_show(time, biomass, substrate, co2):
    fig, ax = plt.subplots()
    ax_2nd = ax.twinx()
    ax_3rd = ax.twinx()

    ax.plot(time, biomass, label='Biomass sim', color='blue')
    ax_2nd.plot(time, substrate, label='Substrate sim', color='orange')
    ax_3rd.plot(time, co2, label='CO2 sim', color='seagreen')

    ax.scatter(df_exp_2['time [h]'], df_exp_2['Biomass [g/L]'], label='Biomass exp', color='dodgerblue')
    ax_2nd.scatter(df_exp_2['time [h]'], df_exp_2['Glucose hplc [g/L]'], label='Glucose conc. exp', color='chocolate')
    ax_3rd.plot(df_exp_2['time [h]'], df_exp_2['Offgas CO2 [g/L]'], label='CO2 exp', color='darkseagreen')
    ax.locator_params(axis='x', nbins=20)

    ax.set_xlabel('time [h]')
    ax.set_ylabel('Biomass [g/L]')
    ax_2nd.set_ylabel('Substrate [g/L]')
    ax_3rd.set_ylabel('CO2 [g/L]', color='seagreen')
    
    ax_3rd.spines['right'].set_position(('outward', 60))

    handles, labels = ax.get_legend_handles_labels()
    handles_2nd, labels_2nd = ax_2nd.get_legend_handles_labels()
    handles_3rd, labels_3rd = ax_3rd.get_legend_handles_labels()
    all_handles = handles + handles_2nd + handles_3rd
    all_labels = labels + labels_2nd + labels_3rd

    # Create a single legend using the combined handles and labels
    ax.legend(all_handles, all_labels, loc='upper center', bbox_to_anchor=(0.5, 1.15), ncols=4)