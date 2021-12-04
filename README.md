[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/alexliberzonlab/Offshore_Ulva_reactor_model/HEAD)

# Offshore Ulva st. reactor model
Predicts the growth rate and internal nitrate content of Ulva in open sea offshore reactors by Meiron Zollmann, Tel Aviv University

## How to use it
1. Run two Jupyter notebooks in this specific order to calibrate key parameters: 
- [Notebook 1](/notebooks/1_Indoor_system_optimisation_HighRes_nofSfT.ipynb)
- [Notebook 2](/notebooks/2_Brine_system_optimization_new_constant_miu.ipynb)

2. Run this Jupyter notebook to create the offshore prediction
- [Notebook 3](/notebooks/3_Reverse_engineering_Offshore_Next.ipynb)

## Use pre-defined conda environment
Use: 

    conda create --name --file environment.yml

or:

    conda create --name offshore_model --file conda_requirements.txt

## Future directions
https://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html#





