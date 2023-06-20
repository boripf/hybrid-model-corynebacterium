import pandas as pd
import numpy as np
import yaml
from scipy.optimize import minimize
from sklearn.metrics import mean_squared_error
from C_functions_opti import model_optimization, plot_simulation

# Load experimental data
df_exp = pd.read_csv('fermentation raw data/data_combined.csv')
biomass_exp = df_exp['Biomass [g/L]']
substrate_exp = df_exp['Glucose [g/L]']

# Load parameters from YAML file
with open('config/parameters.yml', 'r') as file:
    param = yaml.safe_load(file)

mu_max = param['mu_max']
X_max = param['X_max']
Ks = param['Ks']
Ks_qs = param['Ks_qs']
Ki = param['Ki']
Yxs = param['Yxs']
qs_max = param['qs_max']
m_s = param['m_s']
lag = param['lag']

# Calculate growth rate and substrate uptake rate
def mu_eq(mu_max, c_glucose, Ks, Ki, c_biomass, X_max):
    # -- MONOD / insert: mu_max, c_glucose, Ks
    mu0 = mu_max * c_glucose / (c_glucose + Ks)
    # -- LOGISTIC / insert: mu_max, c_biomass, X_max
    mu1 = mu_max * (1 - (c_biomass/ X_max)) 
    # -- MONOD + LOGISTIC / insert: mu_max, c_glucose, Ks, c_biomass, X_max
    mu2 = mu_max * (c_glucose / (c_glucose + Ks)) * (1 - (c_biomass/ X_max))
    # -- MONOD + LOGISTIC + INHIBITION / insert: mu_max, c_glucose, Ks, Ki, c_biomass, X_max
    mu3 = mu_max * (c_glucose / (c_glucose + Ks + (c_glucose**2/ Ki))) * (1 - (c_biomass/ X_max))
    return mu0, mu1, mu2, mu3

def qs_eq(qs_max, c_glucose, Ks_qs, Ki, glu_met, Yxs, f_glucose, V, lag):
    # -- MONOD / insert: qs_max, c_glucose, Ks_qs
    qs0 = qs_max * c_glucose / (Ks_qs + c_glucose)
    # -- MONOD + NON COMPETITIVE INHIBITION / insert: qs_max, c_glucose, Ks_qs, Ki, glu_met
    qs1 = qs_max * c_glucose / (Ks_qs + c_glucose) * (Ki / (Ki + glu_met))
    # -- YIELD / insert: Yxs, f_glucose, V
    qs2 = 1/Yxs * f_glucose / V #NOT SURE IF CORRECT because g_S/g_X but not per h
    # -- MONOD + METABOLIZED GLU / insert: qs_max, c_glucose, Ks_qs, glu_met, lag
    qs3 = qs_max * c_glucose / (Ks_qs + c_glucose) * (1 / (np.exp(glu_met * lag)))
    return qs0, qs1, qs2, qs3

# first define mu, qs, dXdt, dSdt -- later extension that we have a list of different mu,qs each and all are tested automatically
mu_all = mu_eq # mu_max, c_glucose, Ks, Ki, c_biomass, X_max
qs_all = qs_eq # qs_max, c_glucose, Ks_qs, Ki, glu_met, Yxs, f_glucose, V, lag

init_parameters = [mu_max, X_max, Ks, Ks_qs, Ki, Yxs, qs_max, m_s, lag]
# random bounds now -- maybe range of 30% of literature value
parameter_bounds = [(0.0,0.6), (15,30), (5,10), (5,10), (0.0,0.1), (0.0,1.0), (0.0,1.0), (0.0,0.5), (0.0,0.1)]

# the predictions are made by using the model
##  only input is the feed (add acid and glucose)
# difference between experimental data and simulation is calculated
# define parameter bounds
# parameter estimation with minimize function

def objective_function(parameters):
    # Solve the model using the optimal parameters
    time_pred, biomass_pred, substrate_pred = model_optimization(param, optimal_parameters)  # Solve the model using the current parameters
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
optimal_parameters = optimize_parameters(init_parameters)

# Solve the model using the optimal parameters
time_pred, biomass_pred, substrate_pred = model_optimization(param, optimal_parameters)

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