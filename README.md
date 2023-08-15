# Hybrid semi-parametric model for fermentation of Corynebacterium glutamicum

### Introduction
Hybrid semi-parametric model that includes a mechanistic model and a data-driven model. Application for the fermentation of Corynebacterium glutamicum.
- advantages of this approach
- purpose of the project

![Structure of hybrid semi-parametric Model](images\hybrid_model_structure.png)

### Structure

![Project timeline](images/TimeLine.png)

Each file starts with a letter which should represent the steps of development within the developing of the model. Each step contains a .py file where the model is defined and other relevant functions. The functions are then imported to the jupyter notebook, where the output is further processed.

- config
    - parameters.yml (fermentation parameters are defined)
- data
    - batch_no1 (raw data, preprocessed dataset, parameter estimation, generated batches)
    - batch_no2 (raw data, preprocessed dataset, parameter estimation, generated batches)
- images

**A** Data Preprocessing.
**B** Mechanistic Models.
**C** Parameter Estimation.
**D** Sensitivity Analysis.
**E** Data Generation
**F** Random Forest
**G** Hybrid Model
- requirements.txt (pip install -r requirements.txt - to install all required packages)

Each section has its .py file that contains the mechanistic models and other functions and secondly jupyter Notebooks (.ipynb) that contain the preprocessing and the main code.

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