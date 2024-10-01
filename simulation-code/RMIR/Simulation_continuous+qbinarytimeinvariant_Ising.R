# Sufficient reductions with mixed types
# continuous predictors X and categorical predictor time-invariant W
rm(list = ls())
library(purrr)
library(RRRR)
library(dplyr)
library(Rcpp)
library(IsingSampler)
library(RcppArmadillo)
sourceCpp("RMIR/globalIsing.cpp")
sourceCpp("RMIR/RMIR_MCEM_update_allcontinuous+qbinary.cpp")
source("RMIR/global_Ising_with_covariates.R")
source("RMIR/ancillary_functions.R")
source("RMIR/IntrinsicMeanGrassMan.R")
source("../../RMIR/RMIR_MCEM_code_allcontinuous+qbinary.R")


# time invariant predictors and y across all clusters 
gen_gamma <- function(n, p, d, SigmaTilde, Gamma0 = NULL){
  if(is.null(Gamma0)){
    Gamma0 <- pracma::randortho(p)[, 1:d, drop = FALSE]    
  } 
  K <- diag(p) - tcrossprod(Gamma0)
  theta <- lapply(1:n, function(k){
    v <- mvtnorm::rmvnorm(d, rep(0, p), SigmaTilde)
    K %*% t(v)
  })
  Gamma <- sapply(theta, expmap, Gamma0, simplify = FALSE) 
  Sigma <- K %*% SigmaTilde %*% K
  return(list(Gamma0 = Gamma0, Gamma = Gamma, Sigma = Sigma))
}

gen_onecluster <- function(m, Gammai, beta, Delta, # arguments for X | W, Y 
                           A = NULL, C1, 
                           C2,  # arguments for W | Y
                           sigmay = 1,
                           muY = 0, muW = 0){
    # n i
    p <- nrow(Gammai);
    d <- ncol(Gammai); # d: dimensions for continuous parts
    Y <- rnorm(m, mean = muY, sd = sqrt(sigmay))
    fy <- cbind(Y)
    # Generating W from conditional Ising model
    gybar <- matrix(colMeans(fy), nrow = 1)
    r2 <- ncol(C2)
    q <- ncol(beta)
    kk <- q * (q + 1) / 2 # q diagonal terms and q(q-1)/2 off-diagonal
    if (nrow(A) != kk) {
      stop("A has the wrong number of rows")
    }
    tau <- A %*% t(C2) # dimension kk \times r
    ACf <- gybar %*% t(tau)
    thresholds <- ACf[1:q]
    graphs <- matrix(0, nrow = q, ncol = q)
    # arr.ind for upper triangular elements
    upper_tri_index <- which(upper.tri(graphs), arr.ind = T)
    lower_tri_index <- upper_tri_index[, c(2, 1)]
    graphs[upper_tri_index] <- ACf[-c(1:q)]
    graphs[lower_tri_index] <- ACf[-c(1:q)]
    W <- IsingSampler(n = 1, graphs, thresholds)
    
    # Generating (time-variant) X from X|W, Y
    
    if (d == 1) {
      X <- fy %*% C1 %*% t(Gammai) +
        pracma::repmat((W - muW) %*% t(beta), m, 1) + mvtnorm::rmvnorm(m, rep(0, p), Delta)
      X <- sweep(X, 2, Gammai, "+")
    } else if (d == 2) {
      fy <- cbind(1 + Y, 1 / 2 + 1 / 2 * Y ^ 2)
      X <- fy %*% t(Gammai) +
        pracma::repmat((W - muW) %*% t(beta), m, 1) + mvtnorm::rmvnorm(m, rep(0, p), Delta)
    }
    
    return(list(
      X = X,
      Y = as.matrix(Y),
      W = W,
      beta = beta
    ))
}

gen_allclusters <- function(n, p, d, q, 
                            C1 = matrix(c(3), ncol = 1), 
                            C2 = matrix(c(6), ncol = 1), 
                            beta, Delta, A, SigmaTilde){
  # if (d==1){
  #   Gamma0 = matrix(rep(1, p)/sqrt(p), ncol = d)
  # }
  v <- gen_gamma(n, p, d, SigmaTilde)
  Gamma0 <- v[[1]]; Gamma <- v[[2]]; Sigma <- v[[3]]
  clusterSize <- sample(15:20, size = n, replace = TRUE)
  gen_all_data <- mapply(FUN = function(oneclusterSize, oneGamma){
    gen_onecluster(oneclusterSize, oneGamma, beta, Delta, A, C1, C2)
  }, clusterSize, Gamma, SIMPLIFY = FALSE)
  
  X <- Reduce("rbind", gen_all_data %>% map(1))
  Y <- Reduce("rbind", gen_all_data %>% map(2)) 
  W <- Reduce("rbind", gen_all_data %>% map(3)) 
  Y <- data.frame(Y = Y, cluster = rep(1:n, clusterSize))
  return(list(Gamma0 = Gamma0, C1 = C1, C2 = C2,
              X = X, Y = Y, W = W,
              Gamma = Gamma,  
              Sigma = Sigma))
}

sim1 <- function(n, d = 1, p = 7, q = 4, 
                 kk = q*(q+1)/2, kq = q*(q-1)/2,
                 Delta = 0.5^abs(outer(1:p, 1:p,  "-")),
                 SigmaTilde,
                 A = matrix(c(rep(1, kk-2), 10, 10))){
  #browser()
  A <- A/sqrt(sum(A^2))                
  beta = matrix(rbinom(p*q, size = 1, prob = 0.2) * sample(c(-1, 1), size = p*q, replace = TRUE), p, q)
  dat <- gen_allclusters(n, p, d, q, beta = beta, Delta = Delta, A = A, SigmaTilde = SigmaTilde)
  X <- dat$X; Y = dat$Y; W = dat$W; 
  clusterSize = as.vector(table(Y$cluster))
  # true theta_i continuous part
  trueThetai <- lapply(dat$Gamma, function(Gamma){
    rbind(solve(Delta, Gamma),
          -t(beta) %*% solve(Delta, Gamma)
          )
  })
  trueTheta0ContinuousPart <- solve(Delta, projection(dat$Gamma0))
  
  trueTheta0BinaryPart <- 
    rbind(
      cbind(crossprod(beta, solve(Delta, dat$Gamma0)),  -A[1:q,, drop = F]),
      cbind(matrix(0, kq, ncol(dat$Gamma0)),  A[-c(1:q),, drop = F] )      
    )
  
  # Fit a global Ising model on W|Ybar
  FY <- data.frame(Y$Y)
  FYbar <- Reduce("rbind", tapply(FY, Y$cluster, colMeans, simplify = TRUE))
  fitWandY <- fit_Ising(W, FYbar, control = list(tol = 1e-5, maxIter = 20))
  estA <- fitWandY$estGamma
  estC2 <- fitWandY$estC
  norm(projection(estA) - projection(A))
  
  # Fit a RMIR model on X|Y, W
  # Global fit
  # fit a multivariate reduce rank model on X using MLE
  Wrep <- W[rep(1:n, clusterSize),]
  
  GMIR <- list()
  fy = ldr::bf(Y$Y, degree = 3)
  fitXandWY <- RRRR::RRR(y = X, x = fy, z = Wrep, 
                         mu = TRUE,
                         r = d)
  estbeta <- fitXandWY$D 
  estGamma0 <- qr.Q(qr(fitXandWY$A))
  estC1 <- qr.R(qr(fitXandWY$A)) %*% t(fitXandWY$B)
  estSigma <- fitXandWY$Sigma
  Matrix::norm(projection(estGamma0) - projection(dat$Gamma0), "F")
  
  GMIR$estDelta <- estSigma
  GMIR$estbeta <- fitXandWY$D 
  GMIR$estGamma0 <- estGamma0
  # Need to update this
  GMIR$Sigma <- NA 
  GMIR$estGamma <- NA
  GMIR$estBasisContinuousPart <- solve(GMIR$estDelta) %*% projection(estGamma0)
  GMIR$estBasisBinaryPart <- 
    rbind(
      cbind(crossprod(GMIR$estbeta, solve(GMIR$estDelta, estGamma0)),  - estA[1:q,, drop = F]),
      cbind(matrix(0, kq, ncol(estGamma0)),  estA[-c(1:q),, drop = F])      
    )
  
  ### separate fit fto get initial estimates for Gamma and SigmaStar
  SMIR <- list()
  n <- max(Y$cluster)
  sep_fit <- lapply(1:n, function(ii){
    try({
      Xi <- X[Y$cluster == ii, ]
      Yi <- Y[Y$cluster == ii,]$Y
      fyi <- fy[Y$cluster == ii,, drop = FALSE]
      ni <- nrow(Xi)
      #Wi <- Wrep[Y$cluster == ii,, drop = FALSE] # This is not possible, since everything is a constant
      fit <- RRRR::RRR(y = Xi, x = fyi, z = NULL, r = d)
        sep_fit_Gamma <-  qr.Q(qr(fit$A))
      sep_fit_C1 <- qr.R(qr(fit$A)) %*% t(fit$B)
      sep_fit_beta <- matrix(0, p, q)
      sep_fit_Delta <- fit$Sigma
      # separate fit basis
      sep_fit_basis <- rbind(solve(sep_fit_Delta, sep_fit_Gamma),
              -t(sep_fit_beta) %*% solve(sep_fit_Delta, sep_fit_Gamma)
        )
      return(list(sep_fit_Gamma, sep_fit_Delta, sep_fit_basis))  
    })
  })
  
  # sep_fit for Gamma0
  sep_fit_Gamma <- lapply(sep_fit, "[[", 1)
  SMIR$estGamma0 <- frechetMeanGr(sep_fit_Gamma)
  SMIR$estV <- lapply(sep_fit_Gamma, invexpmap, Gamma0 = SMIR$estGamma0)
  SMIR$Sigma  <- Reduce("+", lapply(SMIR$estV, tcrossprod))/length(SMIR$estV)
  
  SMIR$estDelta <- Reduce("+", lapply(sep_fit, "[[", 2))/n
  SMIR$estMeanCS <- solve(SMIR$estDelta) %*% projection(SMIR$estGamma0)
  SMIR$estThetai <- lapply(sep_fit, "[[", 3)
  
  # RMIR fit
  # Initialization
  InitEst <- list()
  InitEst$Gamma0 <- GMIR$estGamma0
  InitEst$Delta <- fitXandWY$Sigma
  InitEst$beta <- fitXandWY$D
  InitEst$muX <- t(fitXandWY$mu)
  InitEst$muW <- t(colMeans(W))
  InitEst$Gamma <- sep_fit_Gamma
  InitEst$estV <- lapply(sep_fit_Gamma, invexpmap, Gamma0 = GMIR$estGamma0)
  InitEst$Sigma  <- Reduce("+", lapply(InitEst$estV, tcrossprod))/length(SMIR$estV)
  
  
  fy <- cbind(1, fy)
  InitEst$C <- matrix(rnorm(d * ncol(fy)), nrow = d)
  RMIRfit <- MCEM_update_noGamma0(X, Y, fy, Wrep, InitEst, niter = 40, burnin = 10, 
                                  tol = 0.05, mc.cores = detectCores()/2, B = 800)
  RMIRfit$estBasisContinuousPart <-
    solve(RMIRfit$estDelta, projection(RMIRfit$estGamma0))
  RMIRfit$estBasisBinaryPart <-
    rbind(
      cbind(crossprod(RMIRfit$estbeta, solve(RMIRfit$estDelta, RMIRfit$estGamma0)),  -estA[1:q,, drop = F]),
      cbind(matrix(0, kq, ncol(RMIRfit$estGamma0)),  estA[-c(1:q),, drop = F])      
    )
  
  
  #### Prediction of cluster-specific spaces
  # Final E-step
  doEstepfinal <- Estep_allclusters(X, Y, fy, Wrep, 
                                    RMIRfit$estmuX, RMIRfit$estGamma0, RMIRfit$estC, 
                                    RMIRfit$estmuW,
                                    RMIRfit$estbeta, RMIRfit$estDelta, RMIRfit$estSigma, B = 2000) 
  
  wb <- doEstepfinal$wb
  Gamma_array <- doEstepfinal$Gamma
  Theta_array <- doEstepfinal$Theta # on tangent space at estGamma0
  
  posteriorGamma <- lapply(1:n, function(ii){
    tmp <- Reduce("+", Map("*",  Theta_array, wb[ii, ]))
    gi <- expmap(tmp, RMIRfit$estGamma0)
  })
  
  RMIRfit$estThetai <- lapply(posteriorGamma, function(Gamma){
    rbind(solve(RMIRfit$estDelta, Gamma),
          -t(RMIRfit$estbeta) %*% solve(RMIRfit$estDelta, Gamma)
    )
  })
  ##### Evaluation of estimation
  result <- list()
  ### Mean (Fixed) spaces
  result$estMeanSpaceContinuous = c(
    RMIR = norm(RMIRfit$estBasisContinuousPart - trueTheta0ContinuousPart, "F"),
    GMIR = norm(GMIR$estBasisContinuousPart - trueTheta0ContinuousPart, "F"),
    SMIR = norm(SMIR$estMeanCS - trueTheta0ContinuousPart, "F"))
  result$estSpaceBinary = c(
    RMIR = norm(projection(RMIRfit$estBasisBinaryPart) - projection(trueTheta0BinaryPart), "F"), 
    SMIR = norm(projection(GMIR$estBasisBinaryPart) - projection(trueTheta0BinaryPart), "F"))
  
  #### Random effect covariance
  result$SigmaTilde = c(RMIR = norm(dat$Sigma - RMIRfit$estSigma, "F"),
                        SMIR = norm(dat$Sigma - SMIR$Sigma, "F"))
  
  ### prediction of cluster-specific CS
  RMIR_MSPE <- sapply(1:n, function(ii) {
    norm(projection(RMIRfit$estThetai[[ii]]) - projection(trueThetai[[ii]]),
         "F")
  })
  SMIR_MSPE <- sapply(1:n, function(ii) {
    norm(projection(SMIR$estThetai[[ii]]) - projection(trueThetai[[ii]]), "F")
  })
  result$MSPE <- c(RMIR = mean(RMIR_MSPE), SMIR = mean(SMIR_MSPE))
  return(list(
    result = result,
    #RMIRfit = RMIRfit,
    dat = dat
  ))
}

Settings <- expand.grid(n = c(100, 500, 1000),
                        SigmaTildeType = 1:3)

# Divide into 20 batches, each having 10 runs
id_artemis <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX")) # from 1-80
if(is.na(id_artemis)) id_artemis <- 2
index_table <- (id_artemis -1) %/% 20 + 1


p <- 7
d <- 1
nsim <- 10

n <- Settings$n[index_table]
SigmaTildetype <- Settings$SigmaTildeType[index_table]
if (SigmaTildetype == 2){
  SigmaTilde <- 0.5*0.5^abs(outer(1:p, 1:p, "-"))
}else if (SigmaTildetype == 1){
  SigmaTilde <- 0.5*diag(p)
}else if (SigmaTildetype == 3){
  SigmaTilde <- matrix(0.1, p, p)
  diag(SigmaTilde) <- 0.5 
}

Sim_onesigma <- lapply(1:nsim, function(seed){
  print(seed)
  try(sim1(n = n, SigmaTilde = SigmaTilde))
}) 
lapply(Sim_onesigma, "[[", 1)


save(Sim_onesigma, file = sprintf("./simulation-results/RMIR/timeinvariant/results_sim_index%.2d.Rdata", id_artemis))

