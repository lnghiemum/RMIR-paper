.libPaths("./package")
library(pracma)
library(ldr)
library(abind)
library(Rcpp)
library(RcppArmadillo)
rm(list = ls())

source("RMIR/IntrinsicMeanGrassMan.R")
source("RMIR/ancillary_functions.R")
source("RMIR/RPFC_MCEM_code_unstrSigma.R")
sourceCpp("RMIR/RPFC_MCEM_update.cpp")

### Data generation 
sim1 <- function(n = 300, p = 4, d = 2,
                 minclustersize = 12, 
                 maxclustersize = 20,
                 Sigmatype = 1,
                 Gamma0 = NULL, 
                 seed = 1000, mc.cores = detectCores()/2){
  set.seed(seed)
  if (is.null(Gamma0)){
    A <- matrix(rnorm(p^2), p, p, byrow=TRUE)
    gs <- pracma::gramSchmidt(A)
    Gamma0 <- gs$Q[, 1:d, drop = FALSE]  
  }
  K <- diag(p) - tcrossprod(Gamma0)
  if (Sigmatype == 1){
    Sigma <- 0.5*diag(p)
  }else if (Sigmatype == 2){
    Sigma <- 0.5*0.5^abs(outer(1:p, 1:p, "-"))
  }else if (Sigmatype == 3){
    Sigma <- matrix(0.1, p, p)
    diag(Sigma) <- 0.5
  }
  Sigmastar <- K %*% Sigma %*% K
  theta <- lapply(1:n, function(k){
    v <- mvtnorm::rmvnorm(d, rep(0, p), Sigma)
    K %*% t(v)
  })
  
  Gamma <- sapply(theta, expmap, Gamma0, simplify = F)
  
  # Generate cluster size
  clusterSize <- sample(minclustersize:maxclustersize, 
                        n, replace = TRUE)
  
  if(d==1){
    Delta <- 0.5^abs(outer(1:p, 1:p, "-"))
    datcluster <- lapply(1:n, function(j){
      m <- clusterSize[j]
      y <- rnorm(m)
      y <- data.frame(cluster = j, time = 1:m, yresp=y)
      v_y <- cbind(y$yresp + 1/2*y$yresp^2 + 1/3* y$yresp^3)
      X <- v_y[y$cluster == j, ] %*% t(Gamma[[j]]) + mvtnorm::rmvnorm(m,  sigma = Delta)
      return(list(X=X, Y=y))
    }) 
    trueMeansubspace <- solve(Delta, projection(Gamma0))
  } else if(d==2){
    Delta <- 0.5^abs(outer(1:p, 1:p, "-"))
    datcluster <- lapply(1:n, function(j){
      m <- clusterSize[j]
      y <- rnorm(m)
      y <- data.frame(cluster = j, time = 1:m, yresp=y)
      v_y <- cbind(y$yresp + 1/2*y$yresp^2 + 1/3* y$yresp^3,
                   y$yresp)
      X <- v_y[y$cluster == j, ] %*% t(Gamma[[j]]) + mvtnorm::rmvnorm(m,  sigma = Delta)
      return(list(X=X, Y=y))
    })
    trueMeansubspace <- solve(Delta, projection(Gamma0))
  }
  
  X <- Reduce(rbind, lapply(datcluster, "[[", 1))
  y <- Reduce(rbind, lapply(datcluster, "[[", 2))
  
  fy <- bf(y$yresp, degree = 4) 
  
  ### Global PFC fit
  globalPFCfit <- pfc(X, y$yresp, fy, numdir = d, structure = "unstr")
  globalPFCdir <- solve(globalPFCfit$Deltahat, globalPFCfit$Gammahat)
  
  ### Getting initial estimates for RPFC from separate PFC
  separatePFCfit <- lapply(1:n, function(jj){
    try({
      Xi  <- X[y$cluster == jj, ]
      yi <- y$yresp[y$cluster==jj]
      fyi <- fy[y$cluster == jj, ]
      pfcfit <- ldr::pfc(Xi, yi, fyi, numdir = d, structure = "unstr")
    }, silent = T)  
  })
  
  ind_tmp <- sapply(separatePFCfit, class) == "try-error"
  separatePFCfit[ind_tmp] = list(globalPFCfit) # replace error by global fit
  sepPFCFitDelta <- lapply(separatePFCfit, function(jj) jj$Deltahat)
  sepPFCFitDir <- lapply(separatePFCfit, function(jj){
    a <- solve(jj$Deltahat, projection(jj$Gammahat))
    B <- qr.Q(qr(a))[, 1:d, drop = FALSE]
  })
  
  
  InitEst <- list()
  InitEst$Delta <- globalPFCfit$Deltahat
  InitEst$Gamma0 <- globalPFCfit$Gammahat
  InitEst$C <- globalPFCfit$Betahat
  separatePFCfit_Theta <- 
    lapply(sepPFCFitDir, invexpmap, Gamma0 = InitEst$Gamma0)
  
  Kest <-  diag(p) - tcrossprod(InitEst$Gamma0)
  q <- sapply(separatePFCfit_Theta, tcrossprod, simplify = "array")
  InitEst$Sigma <- MLESigmaSingular(apply(q, 1:2, mean)/d, p-d)
  InitEst$Gamma <- lapply(separatePFCfit, function(jj) jj$Gammahat)
  
  # Running MCEM
  RPFC <- MCEM_update_fixedGamma0(X, y, fy, InitEst, 
                  niter = 50, tol = .05, burnin = 5, mc.cores = mc.cores)
  
  ### Predicting cluster specific central subspaces
  ### Update the last time 
  B <- 5000
  estepfinal <- Estep_allclusters(X, y, fy, RPFC$estGamma0, 
                                  RPFC$estC,
                                  RPFC$estDelta, 
                                  RPFC$estSigma, 
                                  B = B, mc.cores = mc.cores)
  
  
  estGamma <- lapply(1:n, function(ii){
    ## Average on the tangent space
    avgthetai <- Reduce(`+`,Map(`*`, estepfinal$Theta, 
                                estepfinal$wb[ii, ]))
    ## Then project it on the sphere
    gammai<- expmap(avgthetai, RPFC$estGamma0)
    tmp <- solve(RPFC$estDelta, projection(gammai))
    ### Compute projection error
    A1 <- projection(gammai) - projection(Gamma[[ii]])
    B1 <- projection(sepPFCFitDir[[ii]]) - projection(Gamma[[ii]]) 
    
    allerror <- sapply(list(A1, B1), norm, "F")
  })
  estGammaError <- data.frame(Reduce("rbind", estGamma))
  names(estGammaError) <- c("RPFC", "SPFC")
  MSEGamma <- colMeans(estGammaError, na.rm = TRUE)
  
  ## Estimation error on mean central subspace [Gamma0]
  
  RPFC_estMeanSubspace <- solve(RPFC$estDelta, projection(RPFC$estGamma0))
  GPFC_estMeanSubspace <- solve(globalPFCfit$Deltahat, projection(globalPFCfit$Gammahat))
  SPFC_estGamma0  <- frechetMeanGr(sepPFCFitDir)
  SPFC_estMeanSubspace <- projection(SPFC_estGamma0)
  allestMeanSubspace <- list(RPFC = RPFC_estMeanSubspace,
                             GPFC = GPFC_estMeanSubspace,
                             SPFC = SPFC_estMeanSubspace)
  
  
  allestMeanSubspaceError <- 
    sapply(allestMeanSubspace, function(jj) (norm(jj - trueMeansubspace, "F")))
  
  ### Estimation error on Sigma
  SPFC_estV <- lapply(sepPFCFitDir, invexpmap, Gamma0 = SPFC_estGamma0)
  SPFC_estSigma <- Reduce("+", lapply(SPFC_estV, tcrossprod))/length(SPFC_estV)
  
  estSigma <- list(RPFC = RPFC$estSigma, 
                   SPFC = SPFC_estSigma, GPFC = matrix(0, p, p))
  estSigmaError <- sapply(estSigma, function(k) norm(k - Sigma, "F"))
  
  set.seed(NULL)
  
  return(list(Gamma0 = Gamma0, Gamma = Gamma, 
              estMeanSpace = allestMeanSubspace,
              estMeanSpaceError = allestMeanSubspaceError,
              estSigmaError = estSigmaError,
              estGammaError = estGammaError,
              MSEGamma = MSEGamma, 
              estGamma = estGamma))
}  


# For settings with n=100 and n=500 clusters, each setting is run across 10 batches, 
# each batch consisting of 20 samples

Settings <- expand.grid(n = c(100, 500),
                        d = 1:2,
                        Sigmatype = 1:3)


id_artemis <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
if(is.na(id_artemis)){
  id_artemis = 1
}

index_table <- (id_artemis -1) %/% 10 + 1

n <- Settings$n[index_table]
p <- 7
Sigmatype <- Settings$Sigmatype[index_table]
d <- Settings$d[index_table]
nsim <- 20

Sim_onesigma <- lapply(1:nsim, function(seed){
  print(seed)
  try(sim1(n = n, minclustersize = 10, maxclustersize = 15,
                    p = p, Sigmatype = Sigmatype, d=d,
                    Gamma0 = NULL, seed = seed*id_artemis))
})

save(Sim_onesigma, file = sprintf("simulation-results/RPFC-unstrSigma/smalln/results_sim_index%.2d.Rdata", id_artemis))

# # For settings with n=1000, each setting is run across 20 batches,
# # each batch consisting of 10 samples
# 
# Settings <- expand.grid(n = 1000,
#                           d = 1:2,
#                           Sigmatype = 1:3)
# 
# id_artemis <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
# if(is.na(id_artemis)){
#   id_artemis = 1
# }
# index_table <- (id_artemis-1) %/% 20 + 1
# 
# n <- Settings$n[index_table]
# p <- 7
# Sigmatype <- Settings$Sigmatype[index_table]
# d <- Settings$d[index_table]
# nsim <- 10
# 
# Sim_onesigma <- lapply(1:nsim, function(seed){
#   print(seed)
#   try(sim1(n = n, minclustersize = 10, maxclustersize = 15,
#                     p = p, Sigmatype = Sigmatype, d=d,
#                     Gamma0 = NULL, seed = seed*id_artemis))
# })
#save(Sim_onesigma, file = sprintf("simulation-results/RPFC-unstrSigma/bign/results_sim_index%.2d.Rdata", id_artemis))
