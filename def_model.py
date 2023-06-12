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

    # Extract parameters
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
    dS_dt = np.zeros(num_steps)

    # Set initial values
    biomass[0] = X0
    substrate[0] = S0

    # Simulation loop
    for i in range(1, num_steps):
        
        # Calculate growth rate and substrate uptake rate
        mu = mu_max * substrate[i-1] / (substrate[i-1] + Ks)
        qs = qs_max * substrate[i-1] / (Ks_qs + substrate[i-1])

        # Update biomass and substrate concentrations
        dX_dt = mu * biomass[i-1]
        dS_dt[i] = -qs * biomass[i-1]
        biomass[i] = biomass[i-1] + dX_dt * dt
        substrate[i] = substrate[i-1] + dS_dt[i] * dt

    return time, biomass, substrate, dS_dt

def model_basic_logistic_mu(param, t0, t_end, dt):
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
    qs_max = param['qs_max']
    X0 = param['X0']
    S0 = param['S0']
    X_max = param['X_max']

    # Number of time steps
    num_steps = int((t_end - t0) / dt) + 1

    # Arrays to store results
    time = np.linspace(t0, t_end, num_steps)
    biomass = np.zeros(num_steps)
    substrate = np.zeros(num_steps)
    dS_dt = np.zeros(num_steps)

    # Set initial values
    biomass[0] = X0
    substrate[0] = S0

    # Simulation loop
    for i in range(1, num_steps):

        # Calculate growth rate and substrate uptake rate
        mu = mu_max * (1 - (biomass[i-1]/ X_max))
        qs = qs_max * substrate[i-1] / (Ks_qs + substrate[i-1])

        # Update biomass and substrate concentrations
        dX_dt = mu * biomass[i-1]
        dS_dt[i] = -qs * biomass[i-1]
        biomass[i] = biomass[i-1] + dX_dt * dt
        substrate[i] = substrate[i-1] + dS_dt[i] * dt

    return time, biomass, substrate, dS_dt

def model_monod_logistic_mu(param, t0, t_end, dt):
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
    qs_max = param['qs_max']
    X0 = param['X0']
    S0 = param['S0']
    X_max = param['X_max']

    # Number of time steps
    num_steps = int((t_end - t0) / dt) + 1

    # Arrays to store results
    time = np.linspace(t0, t_end, num_steps)
    biomass = np.zeros(num_steps)
    substrate = np.zeros(num_steps)
    dS_dt = np.zeros(num_steps)

    # Set initial values
    biomass[0] = X0
    substrate[0] = S0

    # Simulation loop
    for i in range(1, num_steps):

        # Calculate growth rate and substrate uptake rate
        mu = mu_max * (substrate[i-1] / (substrate[i-1] + Ks)) * (1 - (biomass[i-1]/ X_max))
        qs = qs_max * substrate[i-1] / (Ks_qs + substrate[i-1])

        # Update biomass and substrate concentrations
        dX_dt = mu * biomass[i-1]
        dS_dt[i] = -qs * biomass[i-1]
        biomass[i] = biomass[i-1] + dX_dt * dt
        substrate[i] = substrate[i-1] + dS_dt[i] * dt

    return time, biomass, substrate, dS_dt

def model_monod_inhib_logistic_mu(param, t0, t_end, dt, Ki):
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
    #Ki = param['Ki']
    qs_max = param['qs_max']
    X0 = param['X0']
    S0 = param['S0']
    X_max = param['X_max']

    # Number of time steps
    num_steps = int((t_end - t0) / dt) + 1

    # Arrays to store results
    time = np.linspace(t0, t_end, num_steps)
    biomass = np.zeros(num_steps)
    substrate = np.zeros(num_steps)
    dS_dt = np.zeros(num_steps)

    # Set initial values
    biomass[0] = X0
    substrate[0] = S0

    # Simulation loop
    for i in range(1, num_steps):

        # Calculate growth rate and substrate uptake rate
        mu = mu_max * (substrate[i-1] / (substrate[i-1] + Ks + ((substrate[i-1]**2) / Ki))) * (1 - (biomass[i-1]/ X_max))
        qs = qs_max * substrate[i-1] / (Ks_qs + substrate[i-1])

        # Update biomass and substrate concentrations
        dX_dt = mu * biomass[i-1]
        dS_dt[i] = -qs * biomass[i-1]
        biomass[i] = biomass[i-1] + dX_dt * dt
        substrate[i] = substrate[i-1] + dS_dt[i] * dt

    return time, biomass, substrate, dS_dt

def model_fedbatch_basic(param, t0, t_end, dt):
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
    # Extract experimental data
    df_exp = pd.read_csv('fermentation raw data/data_combined.csv')
    F_glu = df_exp['Glucose feed [L/min]']
    Volume = df_exp['Volume [L]']

    # Extract parameters
    c_gluc_feed = param['c_glu_feed']
    mu_max = param['mu_max']
    X0 = param['X0']
    S0 = param['S0']
    X_max = param['X_max']
    Yxs = param['Yxs']

    # Number of time steps
    num_steps = int((t_end - t0) / dt) + 1

    # Arrays to store results
    time = np.linspace(t0, t_end, num_steps)
    biomass = np.zeros(num_steps)
    substrate = np.zeros(num_steps)
    dS_dt = np.zeros(num_steps)

    # Set initial values
    biomass[0] = X0
    substrate[0] = S0

    # Simulation loop
    for i in range(1, num_steps):
        # Calculate growth rate and substrate uptake rate
        if substrate[i-1] < 0:
             substrate[i-1] = 0

        mu = mu_max * (1 - (biomass[i-1]/ X_max))
        qs = 1/Yxs

        # Update biomass and substrate concentrations
        dX_dt = mu * biomass[i-1]
        dS_dt[i] = ((F_glu[i]*60 / Volume[i]) * (c_gluc_feed - substrate[i-1])) - qs * biomass[i-1]
        biomass[i] = biomass[i-1] + dX_dt * dt
        substrate[i] = substrate[i-1] + dS_dt[i] * dt

    return time, biomass, substrate, dS_dt

def model_fedbatch_monod_glu(param, t0, t_end, dt):
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
    # Extract experimental data
    df_exp = pd.read_csv('fermentation raw data/data_combined.csv')
    F_glu = df_exp['Glucose feed [L/min]']
    Volume = df_exp['Volume [L]']

    # Extract parameters
    c_gluc_feed = param['c_glu_feed']
    mu_max = param['mu_max']
    Ks = param['Ks']
    Ks_qs = param['Ks_qs']
    qs_max = param['qs_max']
    X0 = param['X0']
    S0 = param['S0']
    X_max = param['X_max']

    # Number of time steps
    num_steps = int((t_end - t0) / dt) + 1

    # Arrays to store results
    time = np.linspace(t0, t_end, num_steps)
    biomass = np.zeros(num_steps)
    substrate = np.zeros(num_steps)
    dS_dt = np.zeros(num_steps)

    # Set initial values
    biomass[0] = X0
    substrate[0] = S0

    # Simulation loop
    for i in range(1, num_steps):
        # Calculate growth rate and substrate uptake rate
        if substrate[i-1] < 0:
             substrate[i-1] = 0

        mu = mu_max * (1 - (biomass[i-1]/ X_max))
        qs = qs_max * substrate[i-1] / (Ks_qs + substrate[i-1]) # 1/Yxs

        # Update biomass and substrate concentrations
        dX_dt = mu * biomass[i-1]
        dS_dt[i] = ((F_glu[i]*60 / Volume[i]) * (c_gluc_feed - substrate[i-1])) - qs * biomass[i-1]
        biomass[i] = biomass[i-1] + dX_dt * dt
        substrate[i] = substrate[i-1] + dS_dt[i] * dt

    return time, biomass, substrate, dS_dt

def model_fedbatch_maintenance(param, t0, t_end, dt):
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
    # Extract experimental data
    df_exp = pd.read_csv('fermentation raw data/data_combined.csv')
    F_glu = df_exp['Glucose feed [L/min]']
    Volume = df_exp['Volume [L]']

    # Extract parameters
    c_gluc_feed = param['c_glu_feed']
    mu_max = param['mu_max']
    Ks = param['Ks']
    Ks_qs = param['Ks_qs']
    qs_max = param['qs_max']
    X0 = param['X0']
    S0 = param['S0']
    X_max = param['X_max']
    m_s = param['m_s']

    # Number of time steps
    num_steps = int((t_end - t0) / dt) + 1

    # Arrays to store results
    time = np.linspace(t0, t_end, num_steps)
    biomass = np.zeros(num_steps)
    substrate = np.zeros(num_steps)
    dS_dt = np.zeros(num_steps)

    # Set initial values
    biomass[0] = X0
    substrate[0] = S0

    # Simulation loop
    for i in range(1, num_steps):
        # Calculate growth rate and substrate uptake rate
        if substrate[i-1] < 0:
             substrate[i-1] = 0

        mu = mu_max * (1 - (biomass[i-1]/ X_max))
        qs = qs_max * substrate[i-1] / (Ks_qs + substrate[i-1]) # 1/Yxs

        # Update biomass and substrate concentrations
        dX_dt = mu * biomass[i-1]
        dS_dt[i] = ((F_glu[i]*60 / Volume[i]) * (c_gluc_feed - substrate[i-1])) - qs * biomass[i-1] - m_s * biomass[-1]
        biomass[i] = biomass[i-1] + dX_dt * dt
        substrate[i] = substrate[i-1] + dS_dt[i] * dt

    return time, biomass, substrate, dS_dt

def model_fedbatch_inhib(param, t0, t_end, dt):
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
    # Extract experimental data
    df_exp = pd.read_csv('fermentation raw data/data_combined.csv')
    F_glu = df_exp['Glucose feed [L/min]']
    Volume = df_exp['Volume [L]']

    # Extract parameters
    c_gluc_feed = param['c_glu_feed']
    mu_max = param['mu_max']
    Ki = param['Ki']
    Ks = param['Ks']
    Ks_qs = param['Ks_qs']
    qs_max = param['qs_max']
    X0 = param['X0']
    S0 = param['S0']
    X_max = param['X_max']

    # Number of time steps
    num_steps = int((t_end - t0) / dt) + 1

    # Arrays to store results
    time = np.linspace(t0, t_end, num_steps)
    biomass = np.zeros(num_steps)
    substrate = np.zeros(num_steps)
    dS_dt = np.zeros(num_steps)

    # Set initial values
    biomass[0] = X0
    substrate[0] = S0

    # Simulation loop
    for i in range(1, num_steps):
        # Calculate growth rate and substrate uptake rate
        if substrate[i-1] < 0:
             substrate[i-1] = 0

        mu = mu_max * (1 - (biomass[i-1]/ X_max))
        qs = qs_max * substrate[i-1] / (Ks_qs + substrate[i-1] + (substrate[i-1]**2 / Ki)) # 1/Yxs

        # Update biomass and substrate concentrations
        dX_dt = mu * biomass[i-1]
        dS_dt[i] = ((F_glu[i]*60 / Volume[i]) * (c_gluc_feed - substrate[i-1])) - qs * biomass[i-1]
        biomass[i] = biomass[i-1] + dX_dt * dt
        substrate[i] = substrate[i-1] + dS_dt[i] * dt

    return time, biomass, substrate, dS_dt

def model_inhibition(param, t0, t_end, dt, Ki):
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
        qs = qs_max * c_glucose / (Ks_qs + c_glucose + c) #* (1/np.exp(biomass[-1] * lag))
        mu = qs * Yxs

        # Update biomass and substrate concentrations
        dX_dt = (mu - F_glu[i-1]/Volume[i-1]) * biomass[i-1] # [g/(Lh)]
        dS_dt = (F_glu[i-1] * (c_gluc_feed - c_glucose)/ Volume[i]) - qs * biomass[i-1] # [g/(Lh)]
        biomass[i] = biomass[i-1] + dX_dt * dt
        substrate[i] = substrate[i-1] + dS_dt * dt

    return time, biomass, substrate

def plot_simulation(time, biomass, substrate, dS_dt, title):
    # import experimental data
    df_exp = pd.read_csv('fermentation raw data/data_combined.csv')
    fig, ax = plt.subplots()

    ax_2nd = ax.twinx()
    ax.plot(time, biomass, label='Biomass sim', color='blue')
    ax_2nd.plot(time, substrate, label='Substrate sim', color='orange')
    ax.plot(time, dS_dt, color='green', label='dS/dt [g/L]')
    ax.scatter(df_exp['time [h]'], df_exp['Biomass [g/L]'], label='Biomass exp', color='purple')
    plt.xlabel('Time [h]')
    ax.set_ylabel('Biomass [g/L] & dS/dt [g/L]')
    ax_2nd.set_ylabel('Substrate [g/L]')
    ax.legend(loc='upper left')
    ax_2nd.legend(loc='lower left')
    plt.title(title)
    plt.show()