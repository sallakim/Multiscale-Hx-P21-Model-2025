# Multiscale-Hx-P21-Model
Companion code to "Multiscale Computational Modeling of the Cardiopulmonary Consequences of Postnatal Hyperoxia with Implications for Preterm Born Children"

!Flags to note!:

- **line 75 in model_sol**: reduces computational time for sensitivity analysis and optimization by bypassing the steady-state algorithm. Comment OUT to allow for the steady-state algorithm to run. 

- **line 64 and lines 467-470 in model_sol**: try-catch is useful when running the sensitivity analysis and optimization. However, try-catch can mask bugs and errors when simply running the model via RUN_Cardiovascular_Mechanics, in which case it should be commented OUT.

**RUN_Cardiovascular_Mechanics**: runs the selected animals for nominal or optimized parameters. Calls plot_PV_loops to plot the model PV-loop output with the PV-loop data.

**RUN_sensitivity**: runs the local sensitivity analysis.

**RUN_optimization**: runs the local and/or global optimization code.

**RUN_force_Ca**: runs the model and plots the force-calcium curve.

**plot_combined_sensitivities**: Plots boxplots of the sensitivities for each parameter for all Nx and Hx animals. Corresponds to Figure 2 in the companion manuscript.

**plot_optimized_pars**: plots boxplot comparisons of the optimized parameters for Nx vs Hx. Corresponds to Figure 5 in the companion manuscript.

**plot_figure_6_8**: plots myofiber power comparisons as boxplots and scatter plots to understand correlations. Corresponds to Figure 6 and 8 in the companion manuscript.

**plot_curvature**: plots the septal curvature and the relative shortening for the LV, septum, and RV. Corresponds to Figure 7 in the companion manuscript.

Functions: 

- **make_data_structure_P21**: intakes animal id and excel data and generates a corresponding data structure for use in the model.

- **plot_PV_loops**: intakes the experimental data and simulation results to output figures with the experimental data and model simulations together. Corresponds to Figure 4 in the companion manuscript. 

- **plot_pv_timecourse**: similar to plot_PV_loops but plots the pressure and volume signals.  

- **load_outputs**: unpacks the outputs data structure.

Data: 

- **P21_data_input.xlsx**: excel table containing the data for parameterizing the model, extracted from the experimental data. 

- **raw_data** (folder): raw data as .txt files, data smoothing code, smoothed data 

- **data** (folder): processed data for use in model calibration

Folders: 

- **model**: contains functions called by RUN_Cardiovascular_Mechanics that run and solve the model.

- **opt_pars**: stores the optimized parameters results from model calibration in RUN_optimization.

- **outputs**: stores the simulation results from RUN_Cardiovascular_Mechanics 

- **parameters**: contains functions for nominal estimation of the parameters based on the data from P21_data_input.xlsx.

- **sens_results**: stores the sensitivity results from the local sensitivity analysis in RUN_sensitivity.

- **sensitivity**: contains functions used for the local sensitivity analysis in RUN_sensitivity.

- **Figures**: stores figures.


Environment: MATLAB R2024b
