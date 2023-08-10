import numpy as np
import pandas as pd
import yaml
import matplotlib.pyplot as plt

# import experimental data
df_exp = pd.read_csv('data/batch_no1/data_combined.csv')
df_exp_2 = pd.read_csv('data/batch_no2/data_combined.csv')

F_glu = df_exp['Glucose feed [L/min]']*60 # [L/h]
F_in = df_exp['Feed total [L/min]']*60

F_glu_2 = df_exp_2['Glucose feed [L/min]']*60
F_in_2 = df_exp_2['Feed total [L/min]']*60

# Load parameters from YAML file
with open('config/parameters.yml', 'r') as file:
    param = yaml.safe_load(file)

# Calculate growth rate and substrate uptake rate
def mu_eq(qs, Yxs):
    # -- MONOD / insert: mu_max, c_glucose, Ks
    #mu = mu_max * c_glucose / (c_glucose + Ks)
    # -- LOGISTIC / insert: mu_max, c_biomass, X_max
    #mu = mu_max * (1 - (c_biomass/ X_max)) 
    # -- MONOD + LOGISTIC / insert: mu_max, c_glucose, Ks, c_biomass, X_max
    # mu = mu_max * (c_glucose / (c_glucose + Ks)) * (1 - (c_biomass/ X_max))
    # -- MONOD + LOGISTIC + INHIBITION / insert: mu_max, c_glucose, Ks, Ki, c_biomass, X_max
    #mu = mu_max * (c_glucose / (c_glucose + Ks + (c_glucose**2/ Ki))) * (1 - (c_biomass/ X_max))
    mu = qs * Yxs
    return mu

def qs_eq(qs_max, c_glucose, Ks_qs, c_biomass, lag):
    # -- MONOD / insert: qs_max, c_glucose, Ks_qs
    #qs = qs_max * c_glucose / (Ks_qs + c_glucose)
    # -- MONOD + NON COMPETITIVE INHIBITION / insert: qs_max, c_glucose, Ks_qs, Ki, glu_met
    #s = qs_max * c_glucose / (Ks_qs + c_glucose) * (Ki / (Ki + glu_met))‚
    # -- YIELD / insert: mu, Yxs
    # qs = mu/Yxs 
    # -- MONOD + METABOLIZED GLU / insert: qs_max, c_glucose, Ks_qs, c_biomass, lag
    qs = qs_max * c_glucose / (Ks_qs + c_glucose) * (1 / (np.exp(c_biomass * lag)))
    return qs

def co2_eq(YxCO2, mu, dt_dx):
    co2 = YxCO2 * mu * dt_dx
    return co2

def model(delta_t):
    """
    Simulates the fermentation process based on the provided parameters.
    Args:
        delta_t (float): Time step size in seconds.
    Returns:
        time (array): Array of time values.
        biomass (array): Array of biomass concentrations.
        substrate (array): Array of substrate concentrations.
    """

    # Simulation settings
    t0 = 0
    t_end = 45.8
    dt = delta_t / 60  # Convert delta_t to minutes
    num_steps = int((t_end - t0) / dt) + 1  # Number of time steps

    # Extract parameters
    Yxs = param['Yxs']
    Yco2s = param['Yco2_x']
    mu_max = param['mu_max']
    X_max = param['X_max']
    Ks = param['Ks']
    Ks_qs = param['Ks_qs']
    Ki = param['Ki']
    qs_max = param['qs_max']
    m_s = param['m_s']
    lag = param['lag']
    X0 = param['X0']
    S0 = param['S0']
    V0 = param['V0']
    co20 = param['co20']
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
        f_glucose = F_glu[t]
        f_total = F_in[t]

        # since the glucose concentration can't be negative, it is set to zero
        if c_glucose < 0:
            c_glucose = 0

        # Update growth and glucose uptake rate
        qs = qs_eq(qs_max, c_glucose, Ks_qs, c_biomass, lag)
        mu = mu_eq(qs, Yxs)

        # Update biomass and substrate concentrations
        dV_dt = f_total - (0.4 * 60 / num_steps)  # [L/h] not complete -- include samples + evaporation
        dX_dt = mu * c_biomass - (c_biomass * (dV_dt/ vol))   # [gx/(Lh)]
        dS_dt = ((f_glucose / vol) * (c_glu_feed - c_glucose)) - ((qs + m_s) * c_biomass) - (c_glucose * (dV_dt / vol))  # [gs/(Lh)]
        dCO2_dt = Yco2s * qs * c_biomass  - (c_co2 * (dV_dt / vol)) # [g/(Lh)]
        
        biomass[i] = c_biomass + dX_dt * dt  # [gx/L]
        substrate[i] = c_glucose + dS_dt * dt  # [gs/L]
        co2[i] = c_co2 + dCO2_dt * dt  # [g/L]
        volume[i] = vol + dV_dt * dt  # [L]

    return time, biomass, substrate, co2

def plot_simulation(time, biomass, substrate, co2, title):
    fig, ax = plt.subplots(figsize=(12,6))
    ax_2nd = ax.twinx()
    ax_3rd = ax.twinx()

    ax.plot(time, biomass, label='Biomass sim', color='blue')
    ax_2nd.plot(time, substrate, label='Glucose sim', color='orange')
    ax_3rd.plot(time, co2, label='CO2 sim', color='seagreen')
    ax.locator_params(axis='x', nbins=20)
    
    ax.scatter(df_exp['time [h]'], df_exp['Biomass [g/L]'], label='Biomass exp', color='dodgerblue')
    ax_2nd.scatter(df_exp['time [h]'], df_exp['Glucose [g/L]'], label='Glucose exp', color='chocolate')
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
    ax.legend(all_handles, all_labels)

    plt.title(title)
    plt.show()

def model_2parts(delta_t):
    # Simulation settings
    t0 = 0
    t_end = 36.8
    dt = delta_t/60
    num_steps = int((t_end - t0) / dt) + 1 # Number of time steps

    # Extract parameters
    X0 = param['X0_2']
    S0 = param['S0']
    co20 = param['co20']
    V0 = param['V0_2']
    c_glu_feed = param['c_glu_feed']

    set1 = param['set_part1']
    set2 = param['set_part2']

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
        
        # since the glucose concentration can't be negative, it is set to zero
        if c_glucose < 0:
            c_glucose = 0

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

    return time, biomass, substrate, co2

def plot_simulation_2parts(time, biomass, substrate, co2, title):
    fig, ax = plt.subplots()
    ax_2nd = ax.twinx()
    ax_3rd = ax.twinx()

    ax.plot(time, biomass, label='Biomass sim', color='blue')
    ax_2nd.plot(time, substrate, label='Glucose sim', color='orange')
    ax_3rd.plot(time, co2, label='CO2 sim', color='seagreen')
    ax.locator_params(axis='x', nbins=20)
    
    ax.scatter(df_exp_2['time [h]'], df_exp_2['Biomass [g/L]'], label='Biomass exp', color='dodgerblue')
    ax_2nd.scatter(df_exp_2['time [h]'], df_exp_2['Glucose hplc [g/L]'], label='Glucose exp', color='chocolate')
    ax_3rd.plot(df_exp_2['time [h]'], df_exp_2['Offgas CO2 [g/L]'], label='CO2 exp', color='darkseagreen')  

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
    ax.legend(all_handles, all_labels)

    plt.title(title)
    plt.show()