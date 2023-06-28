import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import mean_squared_error

# Load experimental data
df_exp = pd.read_csv('data/data_combined.csv')
biomass_exp = df_exp['Biomass [g/L]']
substrate_exp = df_exp['Glucose [g/L]']

# Load parameters from YAML file
with open('config/parameters.yml', 'r') as file:
    param = yaml.safe_load(file)


def model_sensitivity_analysis(param, parameters):
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
    F_glu = df_exp['Glucose feed [L/min]']*60
    F_in = df_exp['Feed total [L/min]']*60

    # Extract parameters
    X0 = param['X0']
    S0 = param['S0']
    V0 = param['V0']
    c_glu_feed = param['c_glu_feed']

    mu_max = parameters[0]
    X_max = parameters[1]
    Yxs = parameters[2]
    #m_s = parameters[3]
    m_s = 0.18

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
        glu_met = S_met[i-1]
        f_glucose = F_glu[i]
        
        
        # since the glucose concentration can't be negative, it is set to zero
        if c_glucose < 0:
            c_glucose = 0

        # Define equations for growth rate and glucose uptake rate
        mu = mu_max * (1 - (c_biomass/ X_max))
        qs = mu / Yxs

        S_met[i] = qs * c_biomass * dt #[gs/(gx h) * gx/L * h = gs/L]
        
        # Update biomass and substrate concentrations
        dV_dt = F_in[i] - (0.4*60/num_steps) # [L/h] not complete -- include samples + evaporation
        dX_dt = mu * c_biomass - (c_biomass / V) * dV_dt # [gx/(Lh)]
        dS_dt[i] = ((f_glucose / V) * (c_glu_feed - c_glucose)) - (qs + m_s) * c_biomass - (c_glucose / V) * dV_dt
        # [gs/(Lh)]

        biomass[i] = c_biomass + dX_dt * dt # [gx/L]
        substrate[i] = c_glucose + dS_dt[i] * dt # [gs/L]
        volume[i] = V + dV_dt * dt # [L]

    return time, biomass, substrate

# Root mean squared error is the objective function
def objective_function(parameters):
    # Solve the model using the optimal parameters
    time_pred, biomass_pred, substrate_pred = model_sensitivity_analysis(param, parameters)  # Solve the model using the current parameters
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
    return rmse