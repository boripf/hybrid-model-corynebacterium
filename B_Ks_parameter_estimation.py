import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 1try)  to load the substrat and mu data directly 

# Load the data from the A_preprocessing.ipynb file
# doesnt work: df_Yxs = pd.read_csv('A_preprocessing.ipynb') 

# df_Yxs.to_csv('data.csv', index=False)

# Experimental data
# substrate_concentration = df_Yxs['Glucose [g/L]'].values  # Substrate concentration
# growth_rate = df_Yxs['mu [1/h]'].values  # Corresponding growth rate

# 2try)  to load the substrat and mu data directly 

# Read data from Excel file
df = pd.read_excel('data/Yxs_table.xlsx')

# Extract growth rate and substrate concentration data
growth_rate = df['mu 2.2 [1/h]'].values
substrate_concentration = df['Glucose [g/L]'].values

# Check and handle missing values and inf values
missing_values = np.isnan(substrate_concentration) | np.isnan(growth_rate)
inf_values = np.isinf(substrate_concentration) | np.isinf(growth_rate)
valid_values = ~(missing_values | inf_values)

# Filter out invalid values
substrate_concentration = substrate_concentration[valid_values]
growth_rate = growth_rate[valid_values]

# Define the Monod equation
def monod_equation(substrate, μmax, Ks):
    return μmax * substrate / (Ks + substrate)

initial_guess = [0.58, 1.5]  # Initial values for μmax and Ks

# Perform curve fitting
params, params_covariance = curve_fit(monod_equation, substrate_concentration, growth_rate)

# Extract the fitted parameters
μmax_fit = params[0]  # Fitted maximum specific growth rate
Ks_fit = params[1]  # Fitted substrate saturation constant

# Generate a range of substrate concentrations for plotting
substrate_range = np.linspace(0, np.max(substrate_concentration), 100)

# Plot the experimental data and the fitted curve
plt.scatter(substrate_concentration, growth_rate, label='Experimental Data')
plt.plot(substrate_range, monod_equation(substrate_range, μmax_fit, Ks_fit), label='Fitted Curve')
plt.xlabel('Substrate Concentration')
plt.ylabel('Growth Rate')
plt.legend()
plt.show()

# Validate the fitted model
predicted_growth_rate = monod_equation(substrate_concentration, μmax_fit, Ks_fit)
mean_squared_error = np.mean((predicted_growth_rate - growth_rate) ** 2)
r_squared = 1 - (mean_squared_error / np.var(growth_rate))

#print('Fitted Parameters:')
#print('μmax:', μmax_fit)
#print('Ks:', Ks_fit)
#print('R-squared:', r_squared)