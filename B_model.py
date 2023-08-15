import numpy as np
import pandas as pd
import yaml
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error

####-------------------------------------------------------------------------------------------------

# Load parameters from YAML file
with open('config/parameters.yml', 'r') as file:
    param = yaml.safe_load(file)

# Import experimental data
df_exp_1 = pd.read_csv('data/batch_no1/data_combined.csv')
df_exp_2 = pd.read_csv('data/batch_no2/data_combined.csv')

F_glu_1 = df_exp_1['Glucose feed [L/min]']*60 # [L/h]
F_in_1 = df_exp_1['Feed total [L/min]']*60 # [L/h]

F_glu_2 = df_exp_2['Glucose feed [L/min]']*60 # [L/h]
F_in_2 = df_exp_2['Feed total [L/min]']*60 # [L/h]

# experimental data from batch No. 2
biomass_exp_2 = df_exp_2['Biomass [g/L]']
substrate_exp_2 = df_exp_2['Glucose hplc [g/L]']
co2_exp_2 = df_exp_2['Offgas CO2 [g/L]']

####-------------------------------------------------------------------------------------------------
# BATCH NO. 1 - Library of reaction kinetics

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

def model_batch_no1(delta_t):
    """
    Args:
        delta_t (float): Time step size in minutes.
    Returns:
        time (array): Array of time values.
        biomass (array): Array of biomass concentrations.
        substrate (array): Array of substrate concentrations.
        co2 (array): Array of co2 concentrations.
    """

    # Simulation settings
    t0 = 0
    t_end = 45.8
    dt = delta_t / 60  # Convert delta_t to hours
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

    # Initial values
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
        f_glucose = F_glu_1[t]
        f_total = F_in_1[t]

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

        # since the glucose concentration can't be negative, it is set to zero
        if substrate[i] < 0:
            substrate[i] = 0

    return time, biomass, substrate, co2

def plot_simulation_no1(time, biomass, substrate, co2, title):
    fig, ax = plt.subplots()
    ax_2nd = ax.twinx()
    ax_3rd = ax.twinx()

    ax.plot(time, biomass, label='Biomass sim', color='blue')
    ax_2nd.plot(time, substrate, label='Glucose sim', color='orange')
    ax_3rd.plot(time, co2, label='CO2 sim', color='seagreen')
    ax.locator_params(axis='x', nbins=20)
    
    ax.scatter(df_exp_1['time [h]'], df_exp_1['Biomass [g/L]'], label='Biomass exp', color='dodgerblue')
    ax_2nd.scatter(df_exp_1['time [h]'], df_exp_1['Glucose [g/L]'], label='Glucose exp', color='chocolate')
    ax_3rd.plot(df_exp_1['time [h]'], df_exp_1['Offgas CO2 [g/L]- smoothed'], label='CO2 exp', color='darkseagreen')  

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

    ax.legend(all_handles, all_labels)

    plt.title(title)
    plt.show()

####-------------------------------------------------------------------------------------------------
# BATCH NO. 2 - Segment fermentation into 2 phases

def model_batch_no2(delta_t):
    """
    Args:
        delta_t (float): Time step size in minutes.
    Returns:
        time (array): Array of time values.
        biomass (array): Array of biomass concentrations.
        substrate (array): Array of substrate concentrations.
        co2 (array): Array of co2 concentrations.
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

    # Define 2 parameter sets - previously found through part C
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
        
        # switch from parameter set 1 to set 2 when glucose is fed into the reactor
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

    return time, biomass, substrate, co2

def plot_simulation_no2(time, biomass, substrate, co2, title):
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
    ax.legend(all_handles, all_labels, loc='upper center')

    plt.title(title)
    plt.show()

def objective_function_no2(delta_t):
    # Solve the model using the optimal parameters
    time_pred, biomass_pred, substrate_pred, co2_pred = model_batch_no2(delta_t)  # Solve the model using the current parameters
    biomass = pd.concat([biomass_exp_2, pd.Series(biomass_pred)], axis=1, keys=['biomass_exp', 'biomass_pred']).dropna()
    biomass_exp_ = biomass['biomass_exp'].values
    biomass_pred_ = biomass['biomass_pred'].values
    mse_x = mean_squared_error(biomass_exp_, biomass_pred_)  # Calculate mean squared error for biomass

    glucose = pd.concat([substrate_exp_2, pd.Series(substrate_pred)], axis=1, keys=['substrate_exp', 'substrate_pred']).dropna()
    substrate_exp_ = glucose['substrate_exp'].values
    substrate_pred_ = glucose['substrate_pred'].values
    mse_s = mean_squared_error(substrate_exp_, substrate_pred_)  # Calculate mean squared error for substrate

    co2 = pd.concat([co2_exp_2, pd.Series(co2_pred)], axis=1, keys=['co2_exp', 'co2_pred']).dropna()
    co2_exp_ = co2['co2_exp'].values
    co2_pred_ = co2['co2_pred'].values
    mse_co2 = mean_squared_error(co2_exp_, co2_pred_)  # Calculate mean squared error for substrate
    
    # Calculate the combined rmse
    mse = (mse_x + mse_s + mse_co2) / 3
    rmse = np.sqrt(mse)  # Calculate root mean squared error
    return rmse

####-------------------------------------------------------------------------------------------------
# BATCH NO. 2 - without segmentation

def objective_function_no2_est(parameters):
    # Solve the model using the optimal parameters
    time_pred, biomass_pred, substrate_pred, co2_pred = model_no2_est(parameters)  # Solve the model using the current parameters
    biomass = pd.concat([biomass_exp_2, pd.Series(biomass_pred)], axis=1, keys=['biomass_exp', 'biomass_pred']).dropna()
    biomass_exp_ = biomass['biomass_exp'].values
    biomass_pred_ = biomass['biomass_pred'].values
    mse_x = mean_squared_error(biomass_exp_, biomass_pred_)  # Calculate mean squared error for biomass

    glucose = pd.concat([substrate_exp_2, pd.Series(substrate_pred)], axis=1, keys=['substrate_exp', 'substrate_pred']).dropna()
    substrate_exp_ = glucose['substrate_exp'].values
    substrate_pred_ = glucose['substrate_pred'].values
    mse_s = mean_squared_error(substrate_exp_, substrate_pred_)  # Calculate mean squared error for substrate

    co2 = pd.concat([co2_exp_2, pd.Series(co2_pred)], axis=1, keys=['co2_exp', 'co2_pred']).dropna()
    co2_exp_ = co2['co2_exp'].values
    co2_pred_ = co2['co2_pred'].values
    mse_co2 = mean_squared_error(co2_exp_, co2_pred_)  # Calculate mean squared error for substrate
    
    # Calculate the combined rmse
    mse = (mse_x + mse_s + mse_co2) / 3
    rmse = np.sqrt(mse)  # Calculate root mean squared error
    return rmse, time_pred, biomass_pred, substrate_pred, co2_pred

def model_no2_est(parameters):
    """
    Args:
        parameters (list): List containing the model parameters.

    Returns:
        time (array): Array of time values.
        biomass (array): Array of biomass concentrations.
        substrate (array): Array of substrate concentrations.
        co2 (array): Array of co2 concentrations.
    """

    # Simulation settings
    t0 = 0
    t_end = 36.8
    delta_t = 1
    dt = delta_t/60
    num_steps = int((t_end - t0) / dt) + 1 # Number of time steps

    # Initial values
    X0 = param['X0_2']
    S0 = param['S0']
    co20 = param['co20']
    V0 = param['V0_2']
    c_glu_feed = param['c_glu_feed']

    # Set parameters from input
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
        f_glucose = F_glu_2[t]
        f_total = F_in_2[t]
        
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

        # since the glucose concentration can't be negative, it is set to zero
        if substrate[i] < 0:
            substrate[i] = 0

    return time, biomass, substrate, co2