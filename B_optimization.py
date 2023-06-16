import pandas as pd
import numpy as np
import yaml
from scipy.optimize import minimize
from sklearn.metrics import mean_squared_error
from B_model import mu_eq, qs_eq, dXdt, dSdt, model_optimization, plot_simulation

# Load parameters from YAML file
with open('config/parameters.yml', 'r') as file:
    param = yaml.safe_load(file)

# Load experimental data
df_exp = pd.read_csv('fermentation raw data/data_combined.csv')
biomass_exp = df_exp['Biomass [g/L]']
substrate_exp = df_exp['Glucose [g/L]']

mu_max = param['mu_max']
X_max = param['X_max']
Ks = param['Ks']
Ks_qs = param['Ks_qs']
Ki = param['Ki']
Yxs = param['Yxs']
qs_max = param['qs_max']
m_s = param['m_s']
lag = param['lag']

# first define mu, qs, dXdt, dSdt -- later extension that we have a list of different mu,qs each and all are tested automatically
mu_all = mu_eq() # mu_max, c_glucose, Ks, Ki, c_biomass, X_max
qs_all = qs_eq() # qs_max, c_glucose, Ks_qs, Ki, glu_met, Yxs, f_glucose, V, lag

init_parameters = [mu_max, X_max, Ks, Ks_qs, Ki, Yxs, qs_max, m_s, lag]
parameter_bounds = [(0,0), (0,0), (0,0), (0,0), (0,0), (0,0), (0,0), (0,0), (0,0)]

# the predictions are made by using the model
##  only input is the feed (add acid and glucose)
# difference between experimental data and simulation is calculated
# define parameter bounds
# parameter estimation with minimize function

def objective_function(parameters):
    # Solve the model using the optimal parameters
    biomass_pred, substrate_pred = model_optimization(parameters)  # Solve the model using the current parameters
    mse_x = mean_squared_error(biomass_exp, biomass_pred)  # Calculate mean squared error for biomass
    mse_s = mean_squared_error(substrate_exp, substrate_pred)  # Calculate mean squared error for substrate
    mse = (mse_x + mse_s)/2
    rmse = np.sqrt(mse)  # Calculate root mean squared error
    return rmse

def optimize_parameters(init_parameters):
    result = minimize(objective_function, init_parameters, bounds=parameter_bounds)
    optimal_parameters = result.x
    print(result.x)
    return optimal_parameters

# Optimize the parameters
optimal_parameters = optimize_parameters()

# Solve the model using the optimal parameters
time_pred, biomass_pred, substrate_pred = model_optimization(optimal_parameters)

# Calculate and print the RMSE
rmse = objective_function(optimal_parameters)

'''
rmse_overview = ['mu', 'qs', 'mu_max', 'X_max', 'Ks', 'Ks_qs', 'Ki', 'Yxs', 'qs_max', 'm_s', 'lag', 'rmse']
for i in range(len(mu_all)):
    for j in range(len(qs_all)):
        # Optimize the parameters
        optimal_p = optimize_parameters()

        # Solve the model using the optimal parameters
        time_pred, biomass_pred, substrate_pred = model_optimization(optimal_parameters)

        # Calculate and print the RMSE
        rmse = objective_function(optimal_parameters)
        append_list=[i, j, opt_p[0], opt_p[1], opt_p[2], opt_p[3], opt_p[4], opt_p[5], opt_p[6], opt_p[7]]
        rmse_overview = rmse_overview.append(append_list)
        
'''