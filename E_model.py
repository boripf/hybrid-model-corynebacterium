import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml

# Load experimental data
df_exp_1 = pd.read_csv('data/batch_no1/data_combined.csv')
df_exp_2 = pd.read_csv('data/batch_no2/data_combined.csv')

# we can not use that because the total feed contains acid and base
F_glu_1 = df_exp_1['Glucose feed [L/min]']*60 #[L/h]
F_in_1 = df_exp_1['Feed total [L/min]']*60 #[L/h]

F_glu_2 = df_exp_2['Glucose feed [L/min]']*60
F_in_2 = df_exp_2['Feed total [L/min]']*60

# Load parameters from YAML file
with open('config/parameters.yml', 'r') as file:
    param = yaml.safe_load(file)

def model_qs_no1(parameters, delta_t):
    """
    Args:
        parameters (list): List containing the model parameters.
        delta_t (int): Integer that defines the sample frequency

    Returns:
        time (array): Array of time values.
        biomass (array): Array of biomass concentrations.
        substrate (array): Array of substrate concentrations.
        co2 (array): Array of co2 concentrations.
    """

    # Simulation settings
    t0 = 0
    t_end = 45.8
    dt = delta_t/60
    num_steps = int((t_end - t0) / dt) + 1 # Number of time steps

    # Initial values
    X0 = param['X0']
    S0 = param['S0']
    V0 = param['V0']
    c_glu_feed = param['c_glu_feed']

    # Arrays to store results
    time = np.linspace(t0, t_end, num_steps)
    biomass = np.zeros(num_steps)
    substrate = np.zeros(num_steps)
    volume = np.zeros(num_steps)
    uptake_rate = np.zeros(num_steps)
    
    # Set initial values
    biomass[0] = X0
    substrate[0] = S0
    volume[0] = V0
    uptake_rate[0] = np.nan

    # Set parameters from input
    Yxs = parameters[0]
    qs_max = parameters[1]
    Ks = parameters[2]
    m_s = parameters[3]
    lag = parameters[4]

    # Simulation loop
    for i in range(1, num_steps):
        c_glucose = substrate[i-1]
        c_biomass = biomass[i-1]
        vol = volume[i-1]

        # time steps need to be adapted for experimental data input
        t = i * delta_t
        f_glucose = F_glu_1[t]
        f_total = F_in_1[t]

        # Update growth and glucose uptake rate
        qs = qs_max * c_glucose / (Ks + c_glucose) * (1 / (np.exp(c_biomass * lag)))
        mu = qs * Yxs

        # Update biomass and substrate concentrations
        dV_dt = f_total - (0.4*60/num_steps) # [L/h] not complete -- include samples + evaporation
        dX_dt = mu * c_biomass - (c_biomass / vol) * dV_dt # [gx/(Lh)]

        dS_dt = ((f_glucose / vol) * (c_glu_feed - c_glucose)) - ((qs + m_s) * c_biomass) - ((c_glucose / vol) * dV_dt)
        # [gs/(Lh)]

        biomass[i] = c_biomass + dX_dt * dt # [gx/L]
        substrate[i] = c_glucose + dS_dt * dt # [gs/L]
        volume[i] = vol + dV_dt * dt # [L]
        uptake_rate[i] = qs

        # since the glucose concentration can't be negative, it is set to zero
        if substrate[i] < 0:
            substrate[i] = 0

        df = pd.DataFrame()
        df[['time', 'biomass', 'glucose', 'qs']] = pd.DataFrame({
            'time': time, 
            'biomass': biomass, 
            'glucose': substrate,
            'qs': uptake_rate
            })

    return df

def model_S_no2(set1, set2, delta_t):
    """
    Args:
        set1 (list): List containing the model parameters for batch phase.
        set2 (list): List containing the model parameters for fed-batch phase.
        delta_t (int): Integer that defines the sample frequency

    Returns:
        df (dataframe): Containing time, biomass, substrate and co2.
    """
    # Simulation settings
    t0 = 0
    t_end = 36.8
    dt = delta_t/60
    num_steps = int((t_end - t0) / dt) + 1 # Number of time steps

    # Initial values
    X0 = param['X0_2']
    S0 = param['S0']
    co20 = param['co20']
    V0 = param['V0_2']
    c_glu_feed = param['c_glu_feed']

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
        dS_dt = ((f_glucose / vol) * (c_glu_feed - c_glucose)) - ((qs + p[4]) * c_biomass) - ((c_glucose / vol) * dV_dt)
        # [gs/(Lh)]
        dCO2_dt = p[1] * qs * c_biomass  - (c_co2 * (dV_dt / vol)) # [g/(Lh)]

        biomass[i] = c_biomass + dX_dt * dt # [gx/L]
        substrate[i] = c_glucose + dS_dt * dt # [gs/L]
        co2[i] = c_co2 + dCO2_dt * dt  # [g/L]
        volume[i] = vol + dV_dt * dt # [L]

        # since the glucose concentration can't be negative, it is set to zero
        if substrate[i] < 0:
            substrate[i] = 0

        df = pd.DataFrame()
        df[['time', 'biomass', 'glucose', 'co2']] = pd.DataFrame({
            'time': time, 
            'biomass': biomass, 
            'glucose': substrate,
            'co2': co2
            })

    return df

def show_plot_no1(df):
    fig, ax = plt.subplots()
    ax_2nd = ax.twinx()
    ax_3rd = ax.twinx()

    ax.plot(df['time'], df['biomass'], label='Biomass sim', color='blue')
    ax_2nd.plot(df['time'], df['glucose'], label='Substrate sim', color='orange')
    ax_3rd.plot(df['time'], df['co2'], label='co2 sim', color='seagreen')
    
    ax.scatter(df_exp['time [h]'], df_exp['Biomass [g/L]'], label='Biomass exp', color='dodgerblue')
    ax_2nd.scatter(df_exp['time [h]'], df_exp['Glucose [g/L]'], label='Glucose conc. exp', color='chocolate')
    ax_3rd.plot(df_exp['time [h]'], df_exp['Offgas CO2 [g/L]- smoothed'], label='CO2 exp', color='darkseagreen')  

    ax.set_xlabel('time [h]')
    ax.set_ylabel('Biomass [g/L]', color='blue')
    ax_2nd.set_ylabel('Substrate [g/L]', color='orange')
    ax_3rd.set_ylabel('CO2 [g/L]', color='seagreen')
    
    ax_3rd.spines['right'].set_position(('outward', 60))
    
    handles, labels = ax.get_legend_handles_labels()
    handles_2nd, labels_2nd = ax_2nd.get_legend_handles_labels()
    handles_3rd, labels_3rd = ax_3rd.get_legend_handles_labels()
    all_handles = handles + handles_2nd + handles_3rd
    all_labels = labels + labels_2nd + labels_3rd

    # Create a single legend using the combined handles and labels
    ax.legend(all_handles, all_labels, loc='upper center', bbox_to_anchor=(0.5, 1.15), ncols=4)
    plt.show()

def show_plot_no2(df):
    fig, ax = plt.subplots()
    ax_2nd = ax.twinx()
    ax_3rd = ax.twinx()

    ax.plot(df['time'], df['biomass'], label='Biomass sim', color='blue')
    ax_2nd.plot(df['time'], df['glucose'], label='Substrate sim', color='orange')
    ax_3rd.plot(df['time'], df['co2'], label='co2 sim', color='seagreen')
    
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
    plt.show()