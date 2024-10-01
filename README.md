Code to reproduce simulation results, figures, and real data analysis results from the paper "Random effects model-based sufficient dimension reduction for independent clustered data" by Linh Nghiem and Francis K.C.Hui.

# Organization

It is important to open and run the code in the R project (RMIR.Rproj), since all the paths included in the files are relative to the project folder.

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
  - Running this code with `id_artemis = 1-120` gives the results in the folder /simulation_results/RPFC_unstrSigma/smalln/
  - Commenting the lines 181-206 and running the lines 211-233 with `id_artemis = 1-120`  gives the results in the folder /simulation_results/RPFC_unstrSigma/bign/

The R file `SimulationRPFC_Section5_p=7_isotropicSigma.R` corresponds to the simulation described in Section 5 of the main paper where the true covariance matrix $\mathbf{\Sigma}$ is isotropic and is assumed to be known to have that structure
  - Running this code with `id_artemis = 1-120` gives the results in the folder /simulation_results/RPFC_isotropicSigma/smalln/
  - Commenting the lines 188-215 and running the lines 221-247 with `id_artemis = 1-120`  gives the results in the folder /simulation_results/RPFC_isotropicSigma/bign/

### RMIR folder

The R file `Simulation_allcontinuous+qbinarytimeinvariant_Ising.R` corresponds to the simulation described in section 6.3 of the main paper with time-invariant binary covariates.
  - Running this file with `id_artemis = 1-180` gives the results in the folder /simulation_results/RMIR/timeinvariant

The R file `Simulation_allcontinuous+qbinarytimeinvariant_GLMM.R` corresponds to the simulation described in section 6.3 of the main paper with time-varying binary covariates.
  - Running this file with `id_artemis = 1-180` gives the results in the folder /simulation_results/RMIR/timevariant

### Choice of d folder

The R file `Simulation_dselection.R` corresponds to the simulation described in the Supplementary S4 for selecting the number of dimensions.
  - Running this file with `id_artemis = 1-180` gives the results in the folder /simulation_results/selection-of-d/results

    
## simulation results

The simulation results folder contain all raw results, as described previously. Additionally it contains the code to reproduce the numbers in the tables. 

- Running `RPFC/RPFC-unstrSigma/analyze_results_unstructuredRE.R` reproduces results in Table 1 in the main paper.
- Running `RPFC-isotropicSigma/analyze_results_isotropic.R` reproduces results in Table S1 in the Supplementary Materials.
- Running `RMIR/analyze_results.R` reproduces results in Table 2 in the main paper.
- Running `selection-of-d/analyse_results_choice_of_d.R` reproduces Figure S1 in the Supplementary Materials.

## data application

The **data** folder contains the raw dataset used by the analysis, which is downloaded from [Gapminder public repository](https://www.gapminder.org/data/).

The `gapminder_datawrangler+fit.R` file contains codes to process raw data, fit the model, and save the results to the `gapminder_applications_withbinarycovariates.Rdata` object. 

The `gapminder_results_visualization.R` file contains code to reproduce numerical values reported in Section 7 and Figure 2 of the paper and Figures S2 and S3 in the appendix. These figures are also included in the folder. 

## Figure 1

This folder contains code to generate Figure 1 in the main paper. Four files (`Sphere1.png`, `Sphere2.png`, `TangentSpace1.png`, `TangentSpace2.png`) are generated from running `plots_sphere_tangent.R`. Then all these files are placed together using Google Slides and save it to `SphereTangentSpace.pdf`.


## packages

The codes use the packages `GrassmannOptim` and `ldr` that can't be downloaded directly from CRAN (only previous versions are available there). Before running any codes in this folder, install these packages from sources as

```{r}
install.packages("packages/GrassmannOptim_2.0.1.tar.gz", type = "sources", repos = NULL)
install.packages("packages/ldr_1.3.tar.gz", type = "sources", repos = NULL)
```





