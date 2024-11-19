# P-model and the new soil moisture stress function
# Code to reproduce the work by Mengoli et al. submitted to GCB journal in November 2024.

# P model (Wang et al., 2017; Stocker et al., 2020) at sub-daily timescale (Mengoli et al., 2022). Repository: P-model_subDaily
There are two methods to run the P model at half-hourly timestep: the running mean and the weighted mean, which produce almost the same results. The user can opt to use one of these approaches by running one of the two main scripts: main_script_running_mean; main_script_weighted_mean. 
Both main scripts require a given number of functions and a settings file. While the functions are saved in the 'R' folder, the settings files are stored in the 'settings files' folder.
All the required functions contain their own descriptions in the R scripts.
In the settings files folder there are example files to use according to the different averaging approaches investigated by Mengoli et al. (2022).To run the P model on another site, it is necessary to change the site name, site elevation and data paths in the settings file.

# P-model_new_water_stress_function repository
In the settings_main_script_rm folder there are example files to have a template to use for the 'noon' approach, which is the one that gives better performance (Mengoli et al. 2022).  
The two example file settings are for two different sites used to perform the running mean approach. 
Use the data (FLUXNET for meteo and MODIS for fAPAR) in the folder to obtain the simulated GPP at half-hourly timestep.
Use the main_running_mean_GMD paper file, in the code folder, to obtain the simulations of GPP at half-hourly timestep. This file is the updated one. This code will refer to the updated functions in the R folder (.R_for_main_running_mean_GMD) and the external settings file, one for each site.
The main_for_GMD paper file is the other part of the code to use to reproduce the work submitted to the GCB journal. 
The code, main_for_GMD paper, is a simple code divided into sections. All the information about how to run the code is embedded in the code itself, in the comments. 
There is a first section where the user needs to run the P model with the main_running_mean_GMD paper to obtain the simulations of GPP. The input data to run the P model are FLUXNET data (https://www.nature.com/articles/s41597-020-0534-3) and MODIS data(http://doi.org/10.5281/zenodo.4392703). Then, it is necessary to run the SPLASHv1 model (Devis et al. 2017) at: https://doi.org/10.5281/zenodo.376293, to obtain the soil water content and the aridity index. After having aggregated the sub-daily values of GPP to daily values (using an external function called dailyAggr.R) a quality check filter is applied to remove the observed GPP with poor quality. Then, there is the section to compute the beta theta ratio (observed GPP/predicted GPP) to thus compute the breakpoint regression to obtain the values of the maximum level and the critical threshold of the ratio. A non-linear regression analysis is performed to obtain the parameters of the new function that needs to be applied to the well-watered GPP.

# Author contact: Giulia Mengoli, gmengoli@ic.ac.uk or giulia.mengoli@gmail.com
# Acknowledgement: this project is funded by the ERC-funded project REALM (grant number 787203)

# Please contact the author prior to any usage of the code
