# make table with standard error
library(tidyverse)

## Analyse results for the simulation
analyseResultsMeanAndSE <- function(folder, numCases, numFilesPerCases = NULL){
  if(is.null(numFilesPerCases)){
    numFilesPerCases <- rep(length(folder)/numCases, numCases)
  }
  # Read all files in folder
  allFiles <- lapply(folder, function(j) get(load(j)))    
  resultsALL <- lapply(1:numCases, function(jj){
    print(jj)
    if (jj==1){
      b_index = 1
    } else{
      b_index = cumsum(numFilesPerCases)[jj-1] + 1
    }
    e_index <- b_index + numFilesPerCases[jj] - 1
    files_one_case <- unlist(allFiles[b_index:e_index], recursive = FALSE)
    index_error <- which(sapply(files_one_case, class) != "list")
    if (length(index_error) > 0){
      files_one_case <- files_one_case[-index_error]
    }
    ### Mean 
    MSEGamma <- rowMeans(sapply(files_one_case, function(kk) colMeans(kk$estGammaError, na.rm = TRUE))) 
    MSEMeanCentralSubspace <-  rowMeans(sapply(files_one_case, function(kk) kk$estMeanSpaceError))
    MSESigma <-  rowMeans(sapply(files_one_case, function(kk) kk$estSigmaError))
    ### SD Error
    SEGamma <- apply(sapply(files_one_case, function(kk) colMeans(kk$estGammaError, na.rm = TRUE)), 1, sd) 
    SEMeanCentralSubspace <-  apply(sapply(files_one_case, function(kk) kk$estMeanSpaceError), 1, sd)
    SESigma <-  apply(sapply(files_one_case, function(kk) kk$estSigmaError), 1, sd)
    
    results <- list(MSEGamma = MSEGamma, MSEMeanCentralSubspace = MSEMeanCentralSubspace,
                    MSESigma = MSESigma, SEGamma = SEGamma, SEMeanCentralSubspace = SEMeanCentralSubspace, 
                    SESigma = SESigma)
    return(results)
  }) 
  AllMSEGamma <- t(sapply(resultsALL, "[[", 1))
  AllMSEMeanCentralSubspace <- t(sapply(resultsALL, "[[", 2))
  AllMSESigma <- t(sapply(resultsALL, "[[", 3))
  
  AllSEGamma <- t(sapply(resultsALL, "[[", 4))
  AllSEMeanCentralSubspace <- t(sapply(resultsALL, "[[", 5))
  AllSESigma <- t(sapply(resultsALL, "[[", 6))
  
  df_MSE <- data.frame(AllMSEMeanCentralSubspace, AllMSESigma, AllMSEGamma) 
  df_SE <- data.frame(AllSEMeanCentralSubspace, AllSESigma, AllSEGamma)
  
  return(list(df_MSE, df_SE))
}

replace_repeat <- function(v){
  u <- v
  for (i in 2:length(v)) if(v[i] == v[i-1]) u[i] = NA
  return(u)
}

### unequal cluster size p = 7 unstructured random effect
folder = dir("simulation-results/RPFC-isotropicSigma/smalln", pattern = ".Rdata", full.names = TRUE)
folder = gtools::mixedsort(folder) # sorting with embedded numbers
sim_inverse_isotropic <- analyseResultsMeanAndSE(folder, 12) 

### paste mean and se
sim_inverse_isotropic_all <- matrix(NA, nrow(sim_inverse_isotropic[[1]]), ncol(sim_inverse_isotropic[[1]]))
tmp_matrix <- lapply(sim_inverse_isotropic, as.matrix)
sim_inverse_isotropic_all[] <- paste(format(round(tmp_matrix[[1]], 2), nsmall = 2), " (",
                                   format(round(tmp_matrix[[2]], 2), nsmall = 2), ")", 
                                   sep = "")
colnames(sim_inverse_isotropic_all) <-
  c(paste0(c("RPFC", "GPFC", "SPFC"), "Mean"),
                                   paste0(c("RPFC", "SPFC", "GPFC"), "Sigma"),
                                   "RPFC_Gamma", "SPFCGamma")

Settings <- expand.grid(n = c(100, 500, 1000),
                        d = 1:2,
                        sigmastar2 = c(0.04, 0.10, 0.50)) %>% 
  filter(n!=1000)


df3 = cbind(Settings, sim_inverse_isotropic_all) %>%
  select(sigmastar2, d, n, everything()) %>%
  mutate(d = factor(d, labels = c("M1", "M2"))) %>%
  arrange(sigmastar2, d, n) # %>% 

# Results with n = 1000
folder = dir("simulation-results/RPFC-isotropicSigma/bign", full.names = TRUE)
folder = gtools::mixedsort(folder) # sorting with embedded numbers
sim_inverse_isotropic_n1000 <- analyseResultsMeanAndSE(folder, 6) 

### paste mean and se
sim_inverse_isotropic_all_n1000 <- matrix(NA, nrow(sim_inverse_isotropic_n1000[[1]]), ncol(sim_inverse_isotropic_n1000[[1]]))
tmp_matrix <- lapply(sim_inverse_isotropic_n1000, as.matrix)
sim_inverse_isotropic_all_n1000[] <- paste(format(round(tmp_matrix[[1]], 2), nsmall = 2), " (",
                                   format(round(tmp_matrix[[2]], 2), nsmall = 2), ")", 
                                   sep = "")
colnames(sim_inverse_isotropic_all_n1000) <-
  c(paste0(c("RPFC", "GPFC", "SPFC"), "Mean"),
    paste0(c("RPFC", "SPFC", "GPFC"), "Sigma"),
    "RPFC_Gamma", "SPFCGamma")

Settings2a <- expand.grid(n = 1000,
                        d = 1:2,
                        sigmastar2 = c(0.04, 0.10, 0.50))

df3a = cbind(Settings2a, sim_inverse_isotropic_all_n1000) %>%
  select(sigmastar2, d, n, everything()) %>%
  mutate(d = factor(d, labels = c("M1", "M2"))) %>%
  arrange(sigmastar2, d, n)

df_p7_isotropic <- rbind(df3, df3a) %>%
  arrange(sigmastar2, d, n) %>%
  mutate_at(1:3, replace_repeat) %>%
  mutate_at(1:3, as.character)
library(kableExtra)
options(knitr.kable.NA = '')

kbl(df_p7_isotropic, format = "latex", escape = F, digits = 2, row.names = F, 
    booktabs = T, linesep = c("", "", "\\addlinespace"))




