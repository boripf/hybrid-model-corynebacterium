import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml

# Load experimental data
df_exp = pd.read_csv('data/data_combined.csv')

# Load parameters from YAML file
with open('config/parameters.yml', 'r') as file:
    param = yaml.safe_load(file)

def model(parameters):
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
    t_end = 46.1
    dt = 1/60
    num_steps = int((t_end - t0) / dt) + 1 # Number of time steps

    # Extract experimental data
    df_exp = pd.read_csv('data/data_combined.csv')
    # we can not use that because the total feed contains acid and base
    F_glu = df_exp['Glucose feed [L/min]']*60 #[L/h]
    F_in = df_exp['Feed total [L/min]']*60 #[L/h]

    # Extract parameters
    X0 = param['X0']
    S0 = param['S0']
    V0 = param['V0']
    c_glu_feed = param['c_glu_feed']

    Yxs = parameters[0]
    qs_max = parameters[1]
    Ks = parameters[2]
    m_s = parameters[3]
    lag = parameters[4]

    # Arrays to store results
    time = np.linspace(t0, t_end, num_steps)
    biomass = np.zeros(num_steps)
    substrate = np.zeros(num_steps)
    volume = np.zeros(num_steps)
    S_met = np.zeros(num_steps)
    dS_dt = np.zeros(num_steps)

    # Set initial values
    biomass[0] = X0
    substrate[0] = S0
    volume[0] = V0

    # Simulation loop
    for i in range(1, num_steps):
        c_glucose = substrate[i-1]
        c_biomass = biomass[i-1]
        V = volume[i-1]
        f_glucose = F_glu[i]
        
        # since the glucose concentration can't be negative, it is set to zero
        if c_glucose < 0:
            c_glucose = 0

        # Update growth and glucose uptake rate
        qs = qs_max * c_glucose / (Ks + c_glucose) * (1 / (np.exp(c_biomass * lag)))
        mu = qs * Yxs

        S_met[i] = qs * c_biomass * dt #[gs/(gx h) * gx/L * h = gs/L]
        
        # Update biomass and substrate concentrations
        dV_dt = F_in[i] - (0.4*60/num_steps) # [L/h] not complete -- include samples + evaporation
        dX_dt = mu * c_biomass - (c_biomass / V) * dV_dt # [gx/(Lh)]

        dS_dt[i] = ((f_glucose / V) * (c_glu_feed - c_glucose)) - (qs + m_s) * c_biomass - (c_glucose / V) * dV_dt
        # [gs/(Lh)]

        biomass[i] = c_biomass + dX_dt * dt # [gx/L]
        substrate[i] = c_glucose + dS_dt[i] * dt # [gs/L]
        volume[i] = V + dV_dt * dt # [L]

        df = pd.DataFrame()
        df[['time', 'biomass', 'glucose']] = pd.DataFrame({'time': time, 'biomass': biomass, 'glucose': substrate})

    return df

def show_plot(df):
    fig, ax = plt.subplots()
    ax_2nd = ax.twinx()

    ax.plot(df['time'], df['biomass'], label='Biomass sim', color='blue')
    ax_2nd.plot(df['time'], df['glucose'], label='Substrate sim', color='orange')
    ax.scatter(df_exp['time [h]'], df_exp['Biomass [g/L]'], label='Biomass exp', color='dodgerblue')
    ax_2nd.scatter(df_exp['time [h]'], df_exp['Glucose [g/L]'], label='Glucose conc. exp', color='chocolate')
    
    ax.set_xlabel('time [h]')
    ax.set_ylabel('Biomass [g/L]')
    ax_2nd.set_ylabel('Substrate [g/L]')

    handles, labels = ax.get_legend_handles_labels()
    handles_2nd, labels_2nd = ax_2nd.get_legend_handles_labels()
    all_handles = handles + handles_2nd
    all_labels = labels + labels_2nd

    # Create a single legend using the combined handles and labels
    ax.legend(all_handles, all_labels, loc='upper center', bbox_to_anchor=(0.5, 1.15), ncols=4)
    plt.show()