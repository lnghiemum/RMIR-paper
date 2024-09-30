Code to reproduce simulation results, figures, and real data analysis results from the paper "Random effects model-based sufficient dimension reduction for independent clustered data" by Linh Nghiem and Francis K.C.Hui.

# Organization

## RMIR

This folder contains an implementation of the RPFC and RMIR methods as described in Section 3 and Section 6 of the paper. 

The `ancillary_functions.R` file contains functions that compute a projection matrix associated with a matrix, generate samples from a singular matrix normal distribution, and perform exponential map and inverse exponential map on a Grassmann manifold. 

The `IntrinsicMeanGrassMan.R` file contains functions to perform geodesic and compute Frechet mean on a Grassmann manifold. 

The two R functions `RPFC_MCEM_code_*.R` implement the MCEM algorithm in the proposed two-stage estimation method for the RPFC model, each corresponding to one case of covariance matrix $\mathbf{\Sigma}$ (isotropic or unstructured). Both functions require a call to the Cpp file `RPFC_MCEM_update.cpp`. 

The R function `RMIR_MCEM_update_allcontinuous+qbinary.R` implements the MCEM algorithm in the proposed two-stage estimation method for the RMIR model with a mix of continuous and (time-invariant/time-varying) binary predictors. This function calls the Cpp file `RMIR_MCEM_update_allcontinuous+qbinary.cpp`.

For the time-invariant binary predictor, a global  Ising model is implemented in the file `global_Ising_with_covariates. R`, which contains a call to `globalIsing.cpp`. 

## simulation code

This folder contains simulation code throughout the paper and Supplementary Material.

### RPFC folder

The R file `SimulationRPFC_Section5_p=7_unstructuredSigma.R` corresponds to the simulation described in Section 5 of the main paper where the covariance matrix $\mathbf{\Sigma}$ is fitted in an unstructured manner.
  - Running this code with `id_artemis = 1-180` gives the results in the folder /simulation_results/RPFC_unstrSigma/smalln/
  - Commenting the lines xxx-xxx and running the lines xxx-xxx with `id_artemis = 1-180`  gives the results in the folder /simulation_results/RPFC_unstrSigma/smalln/

The R file `SimulationRPFC_Section5_p=7_isotropicSigma.R` corresponds to the simulation described in Section 5 of the main paper where the true covariance matrix $\mathbf{\Sigma}$ is isotropic and is assumed to be known to have that structure
  - Running this code with `id_artemis = 1-180` gives the results in the folder /simulation_results/RPFC_isotropicSigma/smalln/
  - Commenting the lines xxx-xxx and running the lines xxx-xxx with `id_artemis = 1-180`  gives the results in the folder /simulation_results/RPFC_isotropicSigma/bign/

### RMIR folder

The R file `Simulation_allcontinuous+qbinarytimeinvariant_Ising.R` corresponds to the simulation described in section 6.3 of the main paper with time-invariant binary covariates.
  - Running this file with `id_artemis = 1-180` gives the results in the folder /simulation_results/RMIR/timeinvariant

The R file `Simulation_allcontinuous+qbinarytimeinvariant_GLMM.R` corresponds to the simulation described in section 6.3 of the main paper with time-varying binary covariates.
  - Running this file with `id_artemis = 1-180` gives the results in the folder /simulation_results/RMIR/timevariant

### Choice of d** folder

The R file `Simulation_dselection.R` corresponds to the simulation described in the Supplementary S4 for selecting the number of dimensions.
  - Running this file with `id_artemis = 1-180` gives the results in the folder /simulation_results/selection-of-d/results

    
## simulation results

## data application

The **data** folder contains the raw dataset used by the analysis, which is downloaded from [Gapminder public repository](https://www.gapminder.org/data/).

The `gapminder_datawrangler+fit.R` file contains codes to process raw data, fit the model, and save the results to the `gapminder_applications_withbinarycovariates.Rdata` object. 

The `gapminder_results_visualization.R` file contains code to reproduce numerical values reported in Section 7 and Figure 2 of the paper and Figures S2 and S3 in the appendix. These figures are also included in the folder. 




