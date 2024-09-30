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




## simulation results

## data application



