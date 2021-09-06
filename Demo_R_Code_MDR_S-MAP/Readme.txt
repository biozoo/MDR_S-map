###################################################################################################################
### Steps for applying MDR S-map to reconstruct large interaction networks with chaotic dynamics  #################
1. Detailed algorithm of MDR S-map and all necessary functions were in the R file recorded in 'Demo_MDR_function.R'
2. To reproduce the main results based on multi-species Ricker model (Fig. 2), all necessary codes and documents were provided in the file, 'MainResults_RickerModels'
3. We parameterized and generated time sereis data by the R files, 'RickerParameter.R', 'TS_RickerModel.R', respectively (in the file 'TimeSeriesGeneration').
4. To execute MDR S-map, we can start from the R script, 'Demo_MDR_Smap_20210625.R'.  
5. To calculate multiview distance, we firstly identified the causal variables for each network nodes by CCM. The results were recorded in 'ccm_rho_XXX' & 'ccm_sig_XXX'
6. We reconstructed numerous multiview SSR (recorded in 'eseleLag_XXX.csv') and then calulated multiview distance via the function 'mvdist.demo'
7. With multiview distance, we performed cross-validation to select the optimal parameter sets (theta, alpha, lambda) that minimizes rMSE (recorded in the file XXX_OptimalCV_Nmvx_Rallx.csv) (This is the most time-consuming step)
8. Based on the optimal parameter set, we performed MDR S-map by using either classical elastic-net or adaptive elastic-net algorithm (default). 
9. Finally, we obtained the interaction Jacobians at each time point ('XXX_jcof_Nmvx_Rallx_demo.csv') as well as the one-step forward forecast skills for each node ('XXX_rout_Nmvx_Rallx_demo.csv') based on MDR S-map  
10. By executing the R sript, 'NetworkComparison_20210625.R', we compared the reconstructed interaction Jacobian networks with the theoretical expectations derived from network models (i.e., Fig. 2 in main text)  

###################################################################################################################
### Tips for applying MDR S-map in more complicated models ########################################################
# The steps 1-10 mention aboved can also be applied to reconstruct interaciton Jacobians in more complicated models, including 
# (ii) discrete Lotka-Volterra competition model (data name=LV), (iii) Ricker-Beverton-Holt model (data name=RBH), (iv) host-parasitoid Nicholson-Bailey model (data name=HP)
# In addition, models incoporating time-varying interaciton coefficients were also included (data name='RK_TV' (Ricker model with time-varying interaction coefficients), 'LV_TV', 'RBH_TV', and 'HP_TV').
# All the data and R script associated with these models were deposited in the file 'Additional models'.
# R codes for parameterizing models and generating time series were included in the file 'TimeSeriesGeneration'.
# Specifically, the parameters obtained from random simulations were included in the file 'Parameters', time series data of population abundance and theoretical interaction Jacobians were included in the file, 'Data' 
# R codes for analyzing time series by MDR S-map were included in the file 'MDR_S-MAP' 
# Finally, we compared the reconstructed interaction Jacobians with their theoretical expectations by executing the R script deposited in the file 'InferenceSkills' based on the results of MDR S-map deposited in the file 'Output'
