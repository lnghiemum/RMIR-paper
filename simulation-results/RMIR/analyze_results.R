# make table with standard error
library(tidyverse)
library(kableExtra)
## Analyse results for the simulation
analyseResultsMeanAndSE <- function(allFiles, numCases, numFilesPerCases = NULL){
  if(is.null(numFilesPerCases)){
    numFilesPerCases <- rep(length(allFiles)/numCases, numCases)
  }
  resultsALL <- lapply(1:numCases, function(jj){
    print(jj)
    if (jj==1){
      b_index = 1
    } else{
      b_index = cumsum(numFilesPerCases)[jj-1] + 1
    }
    e_index <- b_index + numFilesPerCases[jj] - 1
    files_one_case <- unlist(allFiles[b_index:e_index],  recursive = FALSE)
    index_error <- which(sapply(files_one_case, class) == "character")
    if (length(index_error) > 0){
      files_one_case <- files_one_case[-index_error]
    }
    ### Mean 
    rawResults <- sapply(files_one_case, unlist)
    #browser()
    meanResults <- apply(rawResults, 1, mean, na.rm = TRUE) 
    SDresults <- apply(rawResults, 1, sd, na.rm = TRUE)
    
    results <- list(meanResults, SDresults)
    return(results)
  }) 
  
  df_MSE <- t(sapply(resultsALL, "[[", 1))
  df_SE <- t(sapply(resultsALL, "[[", 2))
  
  return(list(df_MSE, df_SE))
}

replace_repeat <- function(v){
  u <- v
  for (i in 2:length(v)) if(v[i] == v[i-1]) u[i] = NA
  return(u)
}

#### Results for binary time-invariant (upper of Table 2)

Settings <- expand.grid(n = c(100, 500, 1000),
                        SigmaTildeType = 1:3)

folder = dir("simulation-results/RMIR/timeinvariant/", pattern = "Rdata", full.names = TRUE)
folder = gtools::mixedsort(folder) # sorting with embedded numbers
# Read all files in folder
allFiles <- lapply(1:length(folder), function(j){
  print(j)
  get(load(folder[j]))
})

sim_inverse_binarytimeinvariant <- analyseResultsMeanAndSE(allFiles, 9) 
### paste mean and se
sim_binarytimeinvariant_final <- matrix(NA, nrow(sim_inverse_binarytimeinvariant[[1]]), ncol(sim_inverse_binarytimeinvariant[[1]]))
tmp_matrix <- lapply(sim_inverse_binarytimeinvariant, as.matrix)
sim_binarytimeinvariant_final[] <- paste0(format(round(tmp_matrix[[1]], 2), nsmall = 2), " (",
                                     format(round(tmp_matrix[[2]], 2), nsmall = 2), ")")
sim_binarytimeinvariant_final <- sim_binarytimeinvariant_final[, -c(4, 5)] # Remove the columns corresponding to binary space that is not of interest


colnames(sim_binarytimeinvariant_final) <-
  c(paste0(c("RPFC", "GPFC", "SPFC"), "MeanContinuous"),
    paste0(c("RPFC", "SPFC"), "Sigma"),
    "RPFC_Gamma", "SPFCGamma")



# Results q binary GLMM
folder = dir("simulation-results/RMIR/timevariant/", 
             pattern = ".Rdata", full.names = TRUE)

folder = gtools::mixedsort(folder) # sorting with embedded numbers
allFiles <- lapply(1:length(folder), function(j){
  print(j)
  get(load(folder[j]))
})


sim_binarytimevariant <- analyseResultsMeanAndSE(allFiles, 9)
sim_binarytimevariant[[1]] <- sim_binarytimevariant[[1]] %>%
  data.frame() %>%
  select(estMeanSpaceContinuous.RMIR, estMeanSpaceContinuous.GMIR, everything())
  
sim_binarytimevariant[[2]] <- sim_binarytimevariant[[2]] %>%
  data.frame() %>%
  select(estMeanSpaceContinuous.RMIR, estMeanSpaceContinuous.GMIR, everything())



### paste mean and se
sim_binarytimevariant_final <- matrix(NA, nrow(sim_binarytimevariant[[1]]), ncol(sim_binarytimevariant[[1]]))
tmp_matrix <- lapply(sim_binarytimevariant, as.matrix)
sim_binarytimevariant_final[] <- paste(format(round(tmp_matrix[[1]], 2), nsmall = 2), " (",
                                     format(round(tmp_matrix[[2]], 2), nsmall = 2), ")", 
                                     sep = "")

colnames(sim_binarytimevariant_final) <-
  c(paste(c("GPFC", "RPFC", "SPFC"), "MeanContinuous", sep = "_"),
    paste(c("RPFC", "SPFC"), "Sigma", sep = "_"),
    "RPFC_Gamma", "SPFC_Gamma")


options(knitr.kable.NA = '')

cbind(Settings, sim_binarytimeinvariant_final) %>% 
  select(SigmaTildeType, n, everything()) %>% 
  mutate(SigmaTildeType = factor(SigmaTildeType, 
                                 labels = c("Diagonal", "AR(1)", "Exchangeable"))) %>%
  mutate(across(1, replace_repeat)) %>%
  kbl(format = "latex", row.names = F, 
      booktabs = T, linesep = c("", "", "\\addlinespace")
  )

cbind(Settings, sim_binarytimevariant_final) %>% 
  select(SigmaTildeType, n, everything()) %>% 
  mutate(SigmaTildeType = factor(SigmaTildeType, 
                                 labels = c("Diagonal", "AR(1)", "Exchangeable"))) %>%
  mutate(across(1, replace_repeat)) %>%
  kbl(format = "latex", row.names = F, 
      booktabs = T, linesep = c("", "", "\\addlinespace")
  )