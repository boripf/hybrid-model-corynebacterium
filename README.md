# hybrid-model-corynebacterium
Hybrid semi-parametric model that includes a mechanistic model and a data-driven model. Application for the fermentation of Corynebacterium glutamicum.

### Structure
- config
    - parameters.yml (fermentation parameters are defined)
- fermentation raw data (contains all the raw experimental data)
    - data_combined.csv (contains the preprocessed data with calculated parameters and is used for further processing)
- images (contains plots of simulations)
- 01_preprocessing.ipynb (preprocessing of the raw data, output: data_combined.csv)
- 02_model.py (contains all ODEs for mu, qs, dX/dt, dS/dt and the model itself)
- 03_simulation.ipynb (running the simulation based on the imported model from def_model.py, output: plots in folder images)
- matlab_sim_fedbatch.m (matlab code for a fedbatch simulation)
- requirements.txt (pip install -r requirements.txt - to install all required packages)

For changing the model in def_model.py, one need to define the mu, qs, dX/dt and dSdt. In order to do that choose on of the equations and adapt the input of each function based on your chosen equation (input is given in the description of each equation). Afterwards, one need to change the input within the model() function too. Make sure that the input in the model is the same as in the above equations.
