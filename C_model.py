import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyDOE import lhs
import os

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

def model_optimization(param, parameters, qs_eq, num_qs):
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
        dV_dt = F_in[i] - (0.4*60/num_steps) # [L/h] not complete -- include samples + evaporation
        dX_dt = mu * c_biomass - (c_biomass / V) * dV_dt # [gx/(Lh)]

        ## (qs + m_s)
        dS_dt[i] = ((f_glucose / V) * (c_glu_feed - c_glucose)) - (qs + m_s) * c_biomass - (c_glucose / V) * dV_dt
        # [gs/(Lh)]

        biomass[i] = c_biomass + dX_dt * dt # [gx/L]
        substrate[i] = c_glucose + dS_dt[i] * dt # [gs/L]
        volume[i] = V + dV_dt * dt # [L]

    return time, biomass, substrate, volume

def plot_estimation(time, biomass, substrate, volume, title, plot_name, set_num):
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
    directory = f'data/estimation/0207_1'

    # Create the directory if it doesn't exist
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Save the plot in the created directory
    plt.savefig(os.path.join(directory, plot_name))
    plt.close()