import numpy as np
import pandas as pd
import yaml
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from pyDOE import lhs
import os

# Extract experimental data
df_exp = pd.read_csv('data/data_combined.csv')

## load biomass and substrate from experiment
biomass_exp = df_exp['Biomass [g/L]']
substrate_exp = df_exp['Glucose [g/L]']

## we have 2 inputs for the model
## feed of glucose and total feed that includes acid&base
F_glu = df_exp['Glucose feed [L/min]']*60
F_in = df_exp['Feed total [L/min]']*60

# Load parameters from YAML file
with open('config/parameters.yml', 'r') as file:
    param = yaml.safe_load(file)

####-------------------------------------------------------------------------------------------------

def get_LHS_samples(num_samples, num_parameters, parameter_bounds):
    # Generate Latin hypercube samples
    lhs_samples = lhs(num_parameters, samples=num_samples, criterion='maximin')

    # Rescale the samples to the specified parameter ranges
    rescaled_samples = []
    for i in range(num_parameters):
        # rescaling = lb + lhs_value * (up - lb)
        rescaling = parameter_bounds[i][0] + lhs_samples[:, i] * (parameter_bounds[i][1] - parameter_bounds[i][0])
        rescaled_samples.append(rescaling)
    samples = np.array(rescaled_samples).T
    return samples

# Root mean squared error is the objective function
def objective_function(parameters, qs_eq, num_qs):
    # Solve the model using the optimal parameters
    time_pred, biomass_pred, substrate_pred, volume_pred = model(parameters, qs_eq, num_qs)  # Solve the model using the current parameters
    biomass = pd.concat([biomass_exp, pd.Series(biomass_pred)], axis=1, keys=['biomass_exp', 'biomass_pred']).dropna()
    biomass_exp_ = biomass['biomass_exp'].values
    biomass_pred_ = biomass['biomass_pred'].values
    mse_x = mean_squared_error(biomass_exp_, biomass_pred_)  # Calculate mean squared error for biomass

    glucose = pd.concat([substrate_exp, pd.Series(substrate_pred)], axis=1, keys=['substrate_exp', 'substrate_pred']).dropna()
    substrate_exp_ = glucose['substrate_exp'].values
    substrate_pred_ = glucose['substrate_pred'].values
    mse_s = mean_squared_error(substrate_exp_, substrate_pred_)  # Calculate mean squared error for substrate
    
    # Calculate the combined rmse
    mse = (mse_x + mse_s)/2
    rmse = np.sqrt(mse)  # Calculate root mean squared error
    return rmse, time_pred, biomass_pred, substrate_pred, volume_pred

def model(parameters, qs_eq, num_qs):
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
    delta_t = 1
    dt = delta_t/60
    num_steps = int((t_end - t0) / dt) + 1 # Number of time steps

    # Extract parameters
    X0 = param['X0']
    S0 = param['S0']
    V0 = param['V0']
    c_glu_feed = param['c_glu_feed']

    Yxs = parameters[0]
    qs_max = parameters[1]
    Ks = parameters[2]
    Ki = parameters[3]
    m_s = parameters[4]
    lag = parameters[5]

    # Arrays to store results
    time = np.linspace(t0, t_end, num_steps)
    biomass = np.zeros(num_steps)
    substrate = np.zeros(num_steps)
    volume = np.zeros(num_steps)
    S_met = np.zeros(num_steps)

    # Set initial values
    biomass[0] = X0
    substrate[0] = S0
    volume[0] = V0

    # Simulation loop
    for i in range(1, num_steps):
        c_glucose = substrate[i-1]
        c_biomass = biomass[i-1]
        glu_met = S_met[i-1]
        vol = volume[i-1]

        # time steps need to be adapted for experimental data input
        t = i * delta_t
        f_glucose = F_glu[t]
        f_total = F_in[t]        
        
        # since the glucose concentration can't be negative, it is set to zero
        if c_glucose < 0:
            c_glucose = 0

        # Update growth and glucose uptake rate
        if num_qs == 0:
            qs = qs_eq(qs_max, c_glucose, Ks)
        elif num_qs == 1:
            qs = qs_eq(qs_max, c_glucose, Ks, Ki, glu_met)
        elif num_qs == 2:
            qs = qs_eq(qs_max, c_glucose, Ks, c_biomass, lag)
        
        mu = qs * Yxs

        S_met[i] = qs * c_biomass * dt #[gs/(gx h) * gx/L * h = gs/L]
        
        # Update biomass and substrate concentrations
        dV_dt = f_total - (0.4*60/num_steps) # [L/h] not complete -- include samples + evaporation
        dX_dt = mu * c_biomass - (c_biomass / vol) * dV_dt # [gx/(Lh)]

        dS_dt = ((f_glucose / vol) * (c_glu_feed - c_glucose)) - ((qs + m_s) * c_biomass) - ((c_glucose / vol) * dV_dt)
        # [gs/(Lh)]

        biomass[i] = c_biomass + dX_dt * dt # [gx/L]
        substrate[i] = c_glucose + dS_dt * dt # [gs/L]
        volume[i] = vol + dV_dt * dt # [L]

    return time, biomass, substrate, volume

def plot_save(time, biomass, substrate, volume, title, plot_name, set_num):
    # import experimental data
    df_exp = pd.read_csv('data/data_combined.csv')

    fig, ax = plt.subplots()
    ax_2nd = ax.twinx()

    ax.plot(time, biomass, label='Biomass sim', color='blue')
    ax.plot(time, volume, label='Volume sim', color='red')
    ax_2nd.plot(time, substrate, label='Substrate sim', color='orange')
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
    plt.title(title)
    
    # Save images in the corresponding folder
    directory = f'data/estimation/0207_6'

    # Create the directory if it doesn't exist
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Save the plot in the created directory
    plt.savefig(os.path.join(directory, plot_name))
    plt.close()

####-------------------------------------------------------------------------------------------------

def objective_function_estimation(parameters):
    # Solve the model using the optimal parameters
    time_pred, biomass_pred, substrate_pred, volume_pred = model_estimation(parameters)  # Solve the model using the current parameters
    biomass = pd.concat([biomass_exp, pd.Series(biomass_pred)], axis=1, keys=['biomass_exp', 'biomass_pred']).dropna()
    biomass_exp_ = biomass['biomass_exp'].values
    biomass_pred_ = biomass['biomass_pred'].values
    mse_x = mean_squared_error(biomass_exp_, biomass_pred_)  # Calculate mean squared error for biomass

    glucose = pd.concat([substrate_exp, pd.Series(substrate_pred)], axis=1, keys=['substrate_exp', 'substrate_pred']).dropna()
    substrate_exp_ = glucose['substrate_exp'].values
    substrate_pred_ = glucose['substrate_pred'].values
    mse_s = mean_squared_error(substrate_exp_, substrate_pred_)  # Calculate mean squared error for substrate
    
    # Calculate the combined rmse
    mse = (mse_x + mse_s)/2
    rmse = np.sqrt(mse)  # Calculate root mean squared error
    return rmse, time_pred, biomass_pred, substrate_pred, volume_pred

def model_estimation(parameters):
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
    delta_t = 1
    dt = delta_t/60
    num_steps = int((t_end - t0) / dt) + 1 # Number of time steps

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

    # Set initial values
    biomass[0] = X0
    substrate[0] = S0
    volume[0] = V0

    # Simulation loop
    for i in range(1, num_steps):
        c_glucose = substrate[i-1]
        c_biomass = biomass[i-1]
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

        biomass[i] = c_biomass + dX_dt * dt # [gx/L]
        substrate[i] = c_glucose + dS_dt * dt # [gs/L]
        volume[i] = vol + dV_dt * dt # [L]

    return time, biomass, substrate, volume

def plot_show(time, biomass, substrate):
    # import experimental data
    df_exp = pd.read_csv('data/data_combined.csv')

    fig, ax = plt.subplots()
    ax_2nd = ax.twinx()

    ax.plot(time, biomass, label='Biomass sim', color='blue')
    ax_2nd.plot(time, substrate, label='Substrate sim', color='orange')
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