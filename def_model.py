import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def model_basic(param, t0, t_end, dt):
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
    # Extract data from experiment
    df_exp = pd.read_csv('fermentation raw data/data_combined.csv')
    F_glu = df_exp['Glucose feed [L/h]']
    Volume = df_exp['Volume [L]']

    # Extract parameters
    c_gluc_feed = param['c_glu_feed']
    mu_max = param['mu_max']
    Ks = param['Ks']
    Ks_qs = param['Ks_qs']
    qs_max = param['qs_max']
    Yxs = param['Yxs']
    X0 = param['X0']
    S0 = param['S0']
    lag = param['lag']

    # Number of time steps
    num_steps = int((t_end - t0) / dt) + 1

    # Arrays to store results
    time = np.linspace(t0, t_end, num_steps)
    biomass = np.zeros(num_steps)
    substrate = np.zeros(num_steps)

    # Set initial values
    biomass[0] = X0
    substrate[0] = S0

    # Simulation loop
    for i in range(1, num_steps):

        # Calculate growth rate and substrate uptake rate
        c_glucose = substrate[i-1]
        if c_glucose < 0:
            c_glucose = 0
        mu = mu_max * c_glucose / (c_glucose + Ks)
        qs = qs_max * c_glucose / (Ks_qs + c_glucose) * (1/np.exp(biomass[-1] * lag))

        # Update biomass and substrate concentrations
        dX_dt = mu * biomass[i-1] # [g/(Lh)]
        dS_dt = (F_glu[i] * c_gluc_feed/ Volume[i]) - qs * biomass[i-1] # [g/(Lh)]
        biomass[i] = biomass[i-1] + dX_dt * dt
        substrate[i] = substrate[i-1] + dS_dt * dt

    return time, biomass, substrate

def model_maintenance(param, t0, t_end, dt):
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

    # Extract parameters
    m_s = param['m_s']
    mu_max = param['mu_max']
    Ks = param['Ks']
    Ks_qs = param['Ks_qs']
    qs_max = param['qs_max']
    X0 = param['X0']
    S0 = param['S0']

    # Number of time steps
    num_steps = int((t_end - t0) / dt) + 1

    # Arrays to store results
    time = np.linspace(t0, t_end, num_steps)
    biomass = np.zeros(num_steps)
    substrate = np.zeros(num_steps)

    # Set initial values
    biomass[0] = X0
    substrate[0] = S0

    # Simulation loop
    for i in range(1, num_steps):
        # Calculate growth rate and substrate uptake rate
        c_glucose = substrate[i-1]
        mu = mu_max * c_glucose / (c_glucose + Ks)
        qs = qs_max * c_glucose / (Ks_qs + c_glucose)

        # Update biomass and substrate concentrations
        dX_dt = mu * biomass[i-1]
        dS_dt = -qs * biomass[i-1] - (m_s * biomass[i-1])
        biomass[i] = biomass[i-1] + dX_dt * dt
        substrate[i] = substrate[i-1] + dS_dt * dt

    return time, biomass, substrate

def model_inhibition(param, t0, t_end, dt):
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

    # Extract parameters
    mu_max = param['mu_max']
    Ks = param['Ks']
    Ks_qs = param['Ks_qs']
    Ki = param['Ki']
    qs_max = param['qs_max']
    X0 = param['X0']
    S0 = param['S0']

    # Number of time steps
    num_steps = int((t_end - t0) / dt) + 1

    # Arrays to store results
    time = np.linspace(t0, t_end, num_steps)
    biomass = np.zeros(num_steps)
    substrate = np.zeros(num_steps)

    # Set initial values
    biomass[0] = X0
    substrate[0] = S0

    # Simulation loop
    for i in range(1, num_steps):
        # Calculate growth rate and substrate uptake rate
        c_glucose = substrate[i-1] / biomass[i-1]
        mu = mu_max * c_glucose / (c_glucose + Ks + (c_glucose**2 / Ki))
        qs = qs_max * c_glucose / (Ks_qs + c_glucose)

        # Update biomass and substrate concentrations
        dX_dt = mu * biomass[i-1]
        dS_dt = -qs * biomass[i-1]
        biomass[i] = biomass[i-1] + dX_dt * dt
        substrate[i] = substrate[i-1] + dS_dt * dt

    return time, biomass, substrate

def plot_simulation(time, biomass, substrate):
    # import experimental data
    df_exp = pd.read_csv('fermentation raw data/data_combined.csv')

    plt.plot(time, biomass, label='Biomass sim')
    plt.plot(time, substrate, label='Substrate sim')
    plt.scatter(df_exp['time [h]'], df_exp['Biomass [g/L]'], label='Biomass exp', color='purple')
    plt.xlabel('Time [h]')
    plt.ylabel('Concentration [g/L]')
    plt.legend()
    plt.show()