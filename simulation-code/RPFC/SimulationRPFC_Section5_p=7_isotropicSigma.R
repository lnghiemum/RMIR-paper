.libPaths("./package")
library(pracma)
library(ldr)
library(abind)
library(Rcpp)
library(RcppArmadillo)
rm(list = ls())

source("RMIR/IntrinsicMeanGrassMan.R")
source("RMIR/ancillary_functions.R")
source("RMIR/RPFC_MCEM_code_isotropicSigma.R")
sourceCpp("RMIR/RPFC_MCEM_update.cpp")

### Data generation 
sim_iso_sigma <- function(n = 300, p = 4, d = 2,
                 minclustersize = 12, 
                 maxclustersize = 20,
                 sigma2 = .25,
                 Gamma0 = NULL, seed = 1000, mc.cores = detectCores()/2){
  set.seed(seed)
  if (is.null(Gamma0)){
    A <- matrix(rnorm(p^2), p, p, byrow=TRUE)
    gs <- pracma::gramSchmidt(A)
    Gamma0 <- gs$Q[, 1:d, drop = FALSE]  
  }
  K <- diag(p) - tcrossprod(Gamma0)

  Sigma <- K * sigma2 
  
  theta <- lapply(1:n, function(k){
    v <- mvtnorm::rmvnorm(d, rep(0, p), sigma2*diag(p))
    K %*% t(v)
  })
  a1 <- sapply(theta, tcrossprod, simplify = "array")
  # q <- apply(a1, 1:2, mean)
  # sum(diag(K %*% q))/(p-1)
  # 
  Gamma <- sapply(theta, expmap, Gamma0, simplify = F)
  
  # Generate cluster size
  if (minclustersize == maxclustersize){
    clusterSize <- rep(minclustersize, n)
  } else{
    clusterSize <- sample(minclustersize:maxclustersize, n, replace = TRUE)
  }
  if(d==1){
    Delta <- 0.5^abs(outer(1:p, 1:p, "-"))
    datcluster <- lapply(1:n, function(j){
      m <- clusterSize[j]
      y <- rnorm(m)
      y <- data.frame(cluster = j, time = 1:m, yresp=y)
      v_y <- cbind(y$yresp + 1/2*y$yresp^2 + 1/3* y$yresp^3)
      X <- v_y[y$cluster == j, ,drop = FALSE] %*% t(Gamma[[j]]) + mvtnorm::rmvnorm(m,  sigma = Delta)
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
  
  ### GPFC 
  globalPFCfit <- pfc(X, y$yresp, fy, numdir = d, structure = "unstr")
  globalPFCdir <- solve(globalPFCfit$Deltahat, globalPFCfit$Gammahat)
  
  ### Getting initial estimates from separate PFC
  separatePFCfit <- lapply(1:n, function(jj){
    a <- try({
    Xi  <- X[y$cluster == jj, ]
    yi <- y$yresp[y$cluster==jj]
    fyi <- fy[y$cluster == jj, ]
    pfcfit <- ldr::pfc(Xi, yi, fyi, numdir = d, structure = "unstr")
    }, silent = TRUE)
  })
  indexError_sepPFC <- sapply(separatePFCfit, class) == "try-error"
 
  separatePFCfit[indexError_sepPFC] <- list(globalPFCfit)
  sepPFCFitDelta <- lapply(separatePFCfit, function(jj) jj$Deltahat)
  sepPFCFitDir <- lapply(separatePFCfit, function(jj){
    a <- solve(jj$Deltahat, projection(jj$Gammahat))
    B <- qr.Q(qr(a))[, 1:d, drop = FALSE]
  })
  
  ### RPFC fit
  InitEst <- list()
  InitEst$Delta <- globalPFCfit$Deltahat
  InitEst$Gamma0 <- globalPFCfit$Gammahat
  InitEst$C <- globalPFCfit$Betahat
  
  SPFC_estGamma0  <- frechetMeanGr(sepPFCFitDir[!indexError_sepPFC])
  SPFCinvexpmap <- lapply(sepPFCFitDir[!indexError_sepPFC], invexpmap, Gamma0 = SPFC_estGamma0)
  q <- sapply(SPFCinvexpmap, tcrossprod, simplify = "array")
  Kest = diag(p) - tcrossprod(SPFC_estGamma0)
  SPFCestsigma2 <- sum(diag(Kest %*% apply(q, 1:2, mean)))/(d*(p-d))
  SPFC_estSigma <- SPFCestsigma2 * (diag(p) - tcrossprod(SPFC_estGamma0))
  
  separatePFCfit_Theta <- 
    lapply(sepPFCFitDir[!indexError_sepPFC], invexpmap, Gamma0 = InitEst$Gamma0)
  Kest <-  diag(p) - tcrossprod(InitEst$Gamma0)
  q <- sapply(separatePFCfit_Theta, tcrossprod, simplify = "array")
  
  InitEst$sigma2 <- sum(diag(Kest %*% apply(q, 1:2, mean)))/(d*(p-d))
  InitEst$Gamma <- lapply(separatePFCfit, function(jj) jj$Gammahat)
  
  RPFC <- MCEM_update_fixedGamma0(X, y, fy, InitEst, 
                                          niter = 100, tol = .05, burnin = 5, mc.cores = mc.cores)
  
  ### Predicting cluster-specific central subspaces
  ### Update the last time 
  B <- 5000
  estepfinal <- Estep_allclusters(X, y, fy, RPFC$estGamma0, 
                                              RPFC$estC, RPFC$estDelta, 
                                              RPFC$estsigma2,
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
  
  estsigma2 <- c(RPFC = RPFC$estsigma2, 
                     SPFC = SPFCestsigma2, GPFC = 0)
  estsigma2Error <- (estsigma2 - sigma2)^2
  
  estSigma <- list(RPFC = RPFC$estSigma, 
                       SPFC = SPFC_estSigma, GPFC = matrix(0, p, p))
  estSigmaError <- sapply(estSigma, function(k) norm(k - Sigma, "F"))
  
  set.seed(NULL)
  
  return(list(estMeanSpaceError = allestMeanSubspaceError,
              estsigma2Error = estsigma2Error,
              estSigmaError = estSigmaError,
              estGammaError = estGammaError))
}  


# For settings with n=100 and n=500 clusters, each setting is run across 10 batches, 
# each batch consisting of 20 samples

Settings <- expand.grid(n = c(100, 500),
                        d = 1:2,
                        sigma2 = c(0.04, 0.10, 0.50))


id_artemis <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
if(is.na(id_artemis)){
  id_artemis = 1
}

id_artemis_list <- (1:120)
index_table <- (id_artemis_list[id_artemis] -1) %/% 10 + 1

n <- Settings$n[index_table]
p <- 7
sigma2 <- Settings$sigma2[index_table]
d <- Settings$d[index_table]
nsim <- 20

Sim_onesigma <- lapply(1:nsim, function(seed){
  print(seed)
  try(sim_iso_sigma(n = n, minclustersize = 10, maxclustersize = 15,
                    p = p, sigma2 = sigma2, d=d,
                    Gamma0 = NULL, seed = seed*id_artemis))
})

save(Sim_onesigma, 
     file = sprintf("simulation-results/RPFC-isotropicSigma/smalln/results_sim_index%.2d.Rdata", id_artemis_list[id_artemis]))


# # # For settings with n=1000, each setting is run across 20 batches, 
# # # each batch consisting of 10 samples
# 
# Settings <- expand.grid(n = 1000,
#                         d = 1:2,
#                         sigma2 = c(0.04, 0.10, 0.50))
# 
# 
# id_artemis <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
# if(is.na(id_artemis)){
#   id_artemis = 1
# }
# id_artemis_list <- 1:120
# index_table <- (id_artemis_list[id_artemis] -1) %/% 20 + 1
# 
# n <- Settings$n[index_table]
# p <- 7
# sigma2 <- Settings$sigma2[index_table]
# d <- Settings$d[index_table]
# nsim <- 10
# 
# Sim_onesigma <- lapply(1:nsim, function(seed){
#   print(seed)
#   try(sim_iso_sigma(n = n, minclustersize = 10, maxclustersize = 15,
#                     p = p, sigma2 = sigma2, d=d,
#                     Gamma0 = NULL, seed = seed*id_artemis))
# })
# 
# save(Sim_onesigma, file = 
#        sprintf("simulation-results/RPFC-isotropicSigma/bign/results_sim_index%.2d.Rdata", id_artemis_list[id_artemis]))
