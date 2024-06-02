This subdirectory contains the code for the simulations performed in Moen et al. (2024),	arXiv:2306.04702. 

The user may have to change the vectors "ns" and "ps" to replicate the simulations from the paper. 
The values of these vectors used in the simulations in the paper are commented out. 

- single_changepoint_simulation.R performs the simulation study in Section 4.1. To recreate the simulation from Section E.2. 
  for a single changepoint, set even_spread = FALSE in this script. 
- multiple_changepoint_simulation.R and multiple_changepoint_simulation_runtime.R replicate the simulation study from Section 4.2. 
  to recreate the simulation study from Section E.2. for multiple changepoints, use multiple_changepoint_simulation_unevenmeanchange.R
- single_changepoint_simulation_misspecifiedmodel.R performs the simulation study in Section 4.3. 
- multiple_changepoint_simulation_variants_of_ESAC.R performs the simulation study in Appendix D
