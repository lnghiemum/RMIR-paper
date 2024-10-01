library(dplyr)
.libPaths("../Simulation_Section3/package")
library(pracma)
library(ldr)
library(abind)

rm(list = ls())

source("RMIR/IntrinsicMeanGrassMan.R")
source("RMIR/ancillary_functions.R")

### Data generation 
sim_data_gen <- function(n = 300, p = 4, d = 2,
                         minclustersize = 12, 
                         maxclustersize = 20,
                         Sigmatype = 1, seed = 1000,
                         sigmay2 = 1,
                         Gamma0 = NULL){
  
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
    Sigma <- matrix(0.5, p, p)
    diag(Sigma) <- 0.1
  }
  Sigmastar <- K %*% Sigma %*% K
  
  theta <- lapply(1:n, function(k){
    v <- mvtnorm::rmvnorm(d, rep(0, p), Sigma)
    K %*% t(v)
  })
  
  Gamma <- sapply(theta, expmap, Gamma0, simplify = F)
  
  # Generate cluster size
  clusterSize <- sample(minclustersize:maxclustersize, n, replace = TRUE)
  if (d==1){
    Delta <- 0.2^abs(outer(1:p, 1:p, "-"))
    datcluster <- lapply(1:n, function(j){
      m <- clusterSize[j]
      y <- sqrt(sigmay2)*rnorm(m)
      y <- data.frame(cluster = j, time = 1:m, yresp=y)
      v_y <- cbind(y$yresp + 1/2*y$yresp^2 + 1/3* y$yresp^3)
      X <- v_y[y$cluster == j, ] %*% t(Gamma[[j]]) + mvtnorm::rmvnorm(m,  sigma = Delta)
      return(list(X=X, Y=y))
    }) 
    trueMeansubspace <- solve(Delta, projection(Gamma0))
  } else if(d==2){
    Delta <- 0.2^abs(outer(1:p, 1:p, "-"))
    datcluster <- lapply(1:n, function(j){
      m <- clusterSize[j]
      y <- sqrt(sigmay2)*rnorm(m)
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
  return(list(X = X, y = y))
} 

numdir_selection <- function(X, y, d_max = 5){
  d_cand = 1:d_max
  fy <- bf(y$yresp, case = "poly", degree = d_max) 
  n <- max(y$cluster)
  
  clusterSize = c(table(y$cluster))
  
  gpfc.fit <- pfc(X, y, fy, structure = "unstr", numdir = max(d_cand), numdir.test = TRUE)
  gpfc.fit$loglik2 <- colSums(mapply(function(Gamma, Delta, beta){
    Xc <- scale(X, center = TRUE, scale = FALSE)
    lkh <- Rfast::dmvnorm(Xc - fy %*% t(beta) %*% t(Gamma), mu = rep(0, p), sigma = Delta, logged = TRUE)
  }, gpfc.fit$Gammahat, gpfc.fit$Deltahat, gpfc.fit$Betahat), na.rm = TRUE)
  gpfc.fit$aic2 <- -2*gpfc.fit$loglik2 + 2*gpfc.fit$numpar[-1]
  gpfc.fit$bic2 <- -2*gpfc.fit$loglik2 + log(nrow(X))*gpfc.fit$numpar[-1]
  
  # selecting d based on either AIC or BIC
  d_gpfc_aic <- d_cand[which.min(gpfc.fit$aic2)]
  d_gpfc_bic <- d_cand[which.min(gpfc.fit$bic2)]
  
  sep.pfc <- lapply(1:n, function(ii){
    Xi <- X[y$cluster == ii, ]
    yi <- y[y$cluster == ii, ]
    fii <- fy[y$cluster ==ii, ]
    fit.pfc <- try({
      pfc.fit <- pfc(Xi, yi$yresp, fii, numdir = max(d_cand),
                     structure = "unstr", numdir.test = TRUE)
      # Getting the full log lik (for d > 0)
      pfc.fit$loglik2 <-  mapply(function(Gamma, Delta, beta){
        Xic <- scale(Xi, center = TRUE, scale = FALSE)
        lkh <- mvtnorm::dmvnorm(Xic - fii %*% t(beta) %*% t(Gamma), mean = rep(0, p), sigma = Delta, log = TRUE)
      }, pfc.fit$Gammahat, pfc.fit$Deltahat, pfc.fit$Betahat) %>% colSums()
    }, silent = T)
    
    if (class(fit.pfc) == "try-error"){
          Gamma.hat <- gpfc.fit$Gammahat
          Delta.hat <- gpfc.fit$Deltahat
          beta.hat <- gpfc.fit$Betahat
          Xic <- scale(Xi, center = TRUE, scale = FALSE)
          llh <- colSums(mapply(function(Gamma, Delta, beta){
            lkh <- Rfast::dmvnorm(Xic - fii %*% t(beta) %*% t(Gamma), mu = rep(0, p), sigma = Delta, logged = TRUE)
          }, Gamma.hat, Delta.hat, beta.hat)) 
          return(llh)
    } else{
          return(pfc.fit$loglik2)
    }
  })  
  sumlogllh <- Reduce("+", sep.pfc)
  
  numpar_cov <- p*(p+1)/2
  numpar <- p + (p-d_cand)*d_cand + d_cand*ncol(fy) + numpar_cov
  aic <- -2*sumlogllh +  2*n*numpar
  bic <- -2*sumlogllh + sum(log(clusterSize))*numpar
  
  # Assuming d >=1
  d_spfc_aic <- d_cand[which.min(aic)]
  d_spfc_bic <- d_cand[which.min(bic)]
  
  dselected = c(gaic = d_gpfc_aic, 
                gbic = d_gpfc_bic, saic = d_spfc_aic, 
                sbic = d_spfc_bic)
  selectedsepPFC <- list(
                         LKH = sumlogllh, 
                         numpar = numpar,
                         gaic = gpfc.fit$aic2,
                         gbic = gpfc.fit$bic2,
                         saic = aic, 
                         sbic = bic,
                         dselected = dselected,
                         gpfc = gpfc.fit,
                         spfc = sep.pfc)
  
  return(selectedsepPFC)
}  

Settings <- expand.grid(n = c(100, 500, 750, 1000),
                        d = 1:2,
                        Sigmatype = c(1, 2, 3),
                        sigmay = 2)


id_artemis <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX")) # from 1-80
if (is.na(id_artemis)){
  id_artemis = 1
}
index_table = id_artemis

p <- 7
nsim <- 10

n <- Settings$n[index_table]
Sigmatype <- Settings$Sigmatype[index_table]
typeofModel <- Settings$typeofModel[index_table]
d <- Settings$d[index_table]
sigmay2 <- Settings$sigmay[index_table]

Sim_onesigma <- lapply(1:nsim, function(seed){
    print(seed)
    mydat <- sim_data_gen(n = n, p = p, d = d, 
                          Sigmatype = Sigmatype, 
                          seed = id_artemis*seed,
                          sigmay2 = sigmay2,
                          minclustersize = 15, 
                          maxclustersize = 20)
    d1 <- numdir_selection(mydat$X, mydat$y, d_max = 4)
}) 
  
save(Sim_onesigma, file = sprintf("simulation-results/selection-of-d/results/choice_of_d_index%.2d.Rdata", id_artemis))
