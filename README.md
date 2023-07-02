# hybrid-model-corynebacterium
Hybrid semi-parametric model that includes a mechanistic model and a data-driven model. Application for the fermentation of Corynebacterium glutamicum.

### Workflow with Git
Everyone is working in their own branch (debbi & marc). From there changes can be merged with the main branch. Follow the step by step guide:
- commit and push all changes to your own branch
- open a new **Git Bash** terminal
- switch from your branch to the main branch: **git checkout main** 
- get the latest updates: **git pull**
- merge the 2 branches: **git merge <branch_name>**
- in the vs code window the merge conflict will pop up
    - for each conflict in a file you can see your own and the main branch version
    - you can decide which version you would like to keep
    - you can show the comparison for a better visualization
- stage the solved conflicts
- commit the changes to the main branche
- to continue working in your own branch
    - checkout to your own branch: **git checkout <branch_name>**
    - pull the copy of the main branch: **git pull origin main**
    

### Structure
Each file starts with a letter which should represent the steps of development within the developing of the model. Each step contains a .py file where the model is defined and other relevant functions. The functions are then imported to the jupyter notebook, where the output is further processed.

- config
    - parameters.yml (fermentation parameters are defined)
- data
    - fermentation raw data (contains all the raw experimental data)
    - pre treated (contains modified data sets)
    - data_combined.csv (contains the preprocessed data with calculated parameters and is used for further processing)
    - Yxs_table (contains biomass & glucose concentrations with calculated parameters)
- images (contains plots of simulations)

**A** is looking fundamentally into the data which is firstly preprocessed to be able to work with it and secondly to understand what is happening during the process.
- A_experiment_analysis (analysis of the fermentation based on experimental data)
- A_preprocessing.ipynb (preprocessing of the raw data, output: data_combined.csv)

**B** contains the first most simple model.
- B_model.py (contains all ODEs for mu, qs, dX/dt, dS/dt and the model itself)
- B_simulation.ipynb

**C** is a development of B and contains parameter sampling based on LHS, an automated testing of different equations and an objective function.
- C_model.py 
- C_simulation.ipynb (running the simulation based on the imported model from def_model.py, output: plots in folder images)

**D** captures the sensitivity analysis.
- D_model.py
- D_sensitivity_analysis.ipynb

- matlab_sim_fedbatch.m (matlab code for a fedbatch simulation)
- requirements.txt (pip install -r requirements.txt - to install all required packages)

For changing the model in def_model.py, one need to define the mu, qs, dX/dt and dSdt. In order to do that choose on of the equations and adapt the input of each function based on your chosen equation (input is given in the description of each equation). Afterwards, one need to change the input within the model() function too. Make sure that the input in the model is the same as in the above equations.
