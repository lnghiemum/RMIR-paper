# Sufficient reductions with mixed types
# continuous predictors X and categorical predictor W

utis <- new.env()
library(purrr)
library(RRRR)
library(lme4)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("RMIR/globalIsing.cpp")
sourceCpp("RMIR/RMIR_MCEM_update_allcontinuous+qbinary.cpp")
source("RMIR/global_Ising_with_covariates.R")
source("RMIR/ancillary_functions.R")
source("RMIR/IntrinsicMeanGrassMan.R")
source("RMIR/RMIR_MCEM_code_allcontinuous+qbinary.R")

# Examples follow Bura et al. (2015) Appendix

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

gen_onecluster <- function(m, Gammai, Delta, C1, 
                           b2, Ai, sigmay = 1,
                           mu_Y = 0){
  #browser()
  # n i
  # A: 1 \times q matrix of random effects on the W|Y part, each column corresponding to one binary, each row corresponds to one cluster
  Ai <- matrix(Ai, nrow = 1)
  q <- ncol(Ai)
  p <- nrow(Gammai);
  d <- ncol(Gammai)
  Y <- rnorm(m, mean = 0, sd = sqrt(sigmay))
  fy <- as.matrix(Y)
  # Generating W from q logistic model of W | Y
  probs <- exp(Y %*% Ai)/(1+exp(Y %*% Ai))
  W <- apply(probs, 2, rbinom, n = m, size = 1)
  W <- as.matrix(W)
  # Generating X from X|W, Y
  if (d==1){
    X <- fy %*% C1 %*% t(Gammai) + 
      W %*% t(b2) + mvtnorm::rmvnorm(m, rep(0, p), Delta)
  } else if (d==2){
    fy <- cbind(Y, Y^2)
    X <- fy %*% C1 %*% t(Gammai) + 
      W %*% t(b2) + mvtnorm::rmvnorm(m, rep(0, p), Delta)
  }
  
  return(list(X = X, Y = as.matrix(Y), W = as.matrix(W), b2 = b2))
}
function(){
  n = 100; p = 7; d = 1; q = 4; A0 = matrix(0.5, n, q);
  SigmaTilde = 0.3*diag(p)
  A <- A0 + matrix(mvtnorm::rmvnorm(n, sigma = 0.3*diag(q)), n, q)
  b2 <- matrix(rnorm(p*q), q, p)
}

gen_allclusters <- function(n, p, d, A0, 
                            C1 = matrix(3, ncol = 1),
                            Delta, b2, sigma_a = 0.3, SigmaTilde){
  if (d==2){
    C1 = matrix(c(3, 0, 0, 1/2), ncol = d)
  }
  q <- ncol(A0)
  v <- gen_gamma(n, p, d, SigmaTilde)
  A <- A0 + matrix(mvtnorm::rmvnorm(n, sigma = sigma_a*diag(q)), n, q)
  Gamma0 <- v[[1]]; Gamma <- v[[2]]; Sigma <- v[[3]]
  clusterSize <- sample(15:20, size = n, replace = TRUE)
  gen_all_data <- mapply(FUN = function(oneclusterSize, oneGamma, Ai){
    gen_onecluster(oneclusterSize, oneGamma, Delta, C1, b2, Ai)
  }, clusterSize, Gamma, data.frame(t(A)), SIMPLIFY = FALSE)
  X <- Reduce("rbind", gen_all_data %>% map(1))
  Y <- Reduce("rbind", gen_all_data %>% map(2)) 
  W <- Reduce("rbind", gen_all_data %>% map(3)) 
  Y <- data.frame(Y = Y, cluster = rep(1:n, clusterSize))
  return(list(Gamma0 = Gamma0, X = X, Y = Y, W = W,
              Gamma = Gamma,  A = A,
              Sigma = Sigma))
}


##### Testing a function
sim1 <- function(n, SigmaTilde, 
                 q = 4, d = 1, p = 7, A0 = matrix(0.5, n, q), 
                 sigma_a = 0.6,
                 C1 = matrix(3, ncol = 1), Delta = 0.5^abs(outer(1:p, 1:p, "-"))){
  
  b2 = matrix(rbinom(p*q, size = 1, prob = 0.2) * sample(c(-1, 1), size = p*q, replace = TRUE), p, q)
  dat <- gen_allclusters(n, p, d, A0, C1,
                         Delta, b2, sigma_a, SigmaTilde)
  X <- dat$X; Y = dat$Y; W = dat$W; 
  A = dat$A
  # true dimension reduction
  trueTheta0ContinuousPart <- solve(Delta, projection(dat$Gamma0))
  trueThetai <- lapply(1:n, function(i){
    rbind(solve(Delta, dat$Gamma[[i]])  %*% C1,
          -crossprod(b2, solve(Delta, dat$Gamma[[i]])) %*% C1 + A[i,]
    )
  })
  
  #### Fit a mixed logistic model on W|Y
  df = data.frame(Y, W)
  estAi <- apply(W, 2, function(Wk){
    fitWandY <- lme4::glmer(Wk ~ Y + (-1+Y|cluster), data = df, family = binomial)
    est_ai <- coef(fitWandY)$cluster[, 2]
    #est_a0 <- summary(fitWandY)$coefficients[2, 1]
  })
  
  ##### X | W, Y
  #### Global fit
  # fit a multivariate reduce rank model on X using MLE
  fy = cbind(Y$Y)
  fitXandWY <- RRRR::RRR(y = X, x = fy, z = W, 
                         mu = TRUE,
                         r = d)
  estb2 <- fitXandWY$D # same direction as in b2, but not 
  estA <- fitXandWY$A
  estGamma0 <- qr.Q(qr(estA))
  estC1 <- qr.R(qr(fitXandWY$A)) %*% t(fitXandWY$B)
  estSigma <- fitXandWY$Sigma
  
  GMIR <- list()
  GMIR$estbeta <- fitXandWY$D 
  GMIR$estGamma0 <- estGamma0
  GMIR$Sigma <- NA 
  GMIR$estGamma <- NA
  GMIR$estDelta <- estSigma
  GMIR$estBasisContinuousPart <- solve(GMIR$estDelta, projection(estGamma0))
  ### SMIR FIT
  SMIR <- list()
  n <- max(Y$cluster)
  sep_fit <- lapply(1:n, function(ii){
    Xi <- X[Y$cluster == ii, ]
    fyi <- fy[Y$cluster == ii,, drop = FALSE]
    Wi <- W[Y$cluster == ii,, drop = FALSE]
    fit <- RRRR::RRR(y = Xi, x = fyi, z = Wi, r = d)
    sep_fit_Gamma <-  qr.Q(qr(fit$A))
    sep_fit_Delta <- fit$Sigma
    sep_fit_estC1 <- qr.R(qr(fit$A)) %*% t(fit$B) 
    sep_fit_basis <- rbind(
      solve(sep_fit_Delta, sep_fit_Gamma) %*% sep_fit_estC1,
      -crossprod(fit$D, solve(sep_fit_Delta, sep_fit_Gamma)) %*% sep_fit_estC1  + estAi[ii, ])
    return(list(sep_fit_Gamma, sep_fit_basis, sep_fit_Delta))
  })
  # sep_fit for Gamma0
  sep_fit_Gamma <- lapply(sep_fit, "[[", 1)
  SMIR$estGamma0 <- frechetMeanGr(sep_fit_Gamma)
  SMIR$estV <- lapply(sep_fit_Gamma, invexpmap, Gamma0 = SMIR$estGamma0)
  SMIR$estSigma <- Reduce("+", lapply(SMIR$estV, tcrossprod))/length(SMIR$estV)
  SMIR$estMeanCS <- frechetMeanGr(lapply(sep_fit, "[[", 1))
  
  SMIR$estDelta <- Reduce("+", lapply(sep_fit, "[[", 3))/n
  SMIR$estBasisContinuousPart <- solve(SMIR$estDelta) %*% projection(SMIR$estGamma0)
  ### MCEM fit
  
  ################
  InitEst <- list()
  InitEst$Gamma0 <- qr.Q(qr(estA))
  InitEst$Delta <- fitXandWY$Sigma
  InitEst$beta <- fitXandWY$D
  InitEst$muX <- t(fitXandWY$mu)
  InitEst$muW <- t(colMeans(W))
  InitEst$Gamma <- sep_fit_Gamma
  estV <- lapply(sep_fit_Gamma, invexpmap, Gamma0 = InitEst$Gamma0)
  InitEst$Sigma <- Reduce("+", lapply(estV, tcrossprod))/length(estV)
  
  fy <- cbind(1, Y$Y)
  InitEst$C <- matrix(rnorm(d * ncol(fy)), nrow = d)
  RMIRfit <- MCEM_update_noGamma0(X, Y, fy, W, InitEst, niter = 80, burnin = 10, 
                                  tol = 0.05, B = 800)
  RMIRfit$estBasisContinuousPart <-
    solve(RMIRfit$estDelta, projection(RMIRfit$estGamma0))
  

  # Estimation of some parameters
  
  # Final E-step
  doEstepfinal <- Estep_allclusters(X, Y, fy, W, 
                                    RMIRfit$estmuX, RMIRfit$estGamma0, RMIRfit$estC, 
                                    RMIRfit$estmuW,
                                    RMIRfit$estbeta, RMIRfit$estDelta, RMIRfit$estSigma, B = 2000) 
  
  wb <- doEstepfinal$wb
  Gamma_array <- doEstepfinal$Gamma
  Theta_array <- doEstepfinal$Theta # on tangent space at estGamma0
  
  posteriorGamma <- lapply(1:n, function(ii){
    tmp <- Reduce("+", Map("*",  Theta_array, wb[ii, ]))
    gi <- expmap(tmp, RMIRfit$estGamma0)
    tmpi <- rbind(
      solve(RMIRfit$estDelta, gi) %*% RMIRfit$estC[, -1],
      -crossprod(RMIRfit$estbeta, solve(RMIRfit$estDelta, gi)) %*% RMIRfit$estC[, -1] + estAi[ii,])
  })
  
  ##### Evaluation of estimation
  result <- list()
  ### Mean (Fixed) spaces
  result$estMeanSpaceContinuous = c(
    GMIR = norm(GMIR$estBasisContinuousPart - trueTheta0ContinuousPart, "F"),
    RMIR = norm(RMIRfit$estBasisContinuousPart - trueTheta0ContinuousPart, "F"),
    SMIR = norm(SMIR$estBasisContinuousPart - trueTheta0ContinuousPart, "F"))
  #### Random effect covariance
  result$SigmaTilde = c(RMIR = norm(dat$Sigma - RMIRfit$estSigma, "F"),
                        SMIR = norm(dat$Sigma - SMIR$estSigma, "F"))
  
  ### prediction of cluster-specific CS
  RMIR_MSPE <- sapply(1:n, function(ii){
    norm(projection(posteriorGamma[[ii]]) - projection(trueThetai[[ii]]), "F")
  }) 
  SMIR_MSPE <- sapply(1:n, function(ii){
    norm(projection(sep_fit[[ii]][[2]]) - projection(trueThetai[[ii]]), "F") 
  })
  result$MSPE <- c(RMIR = mean(RMIR_MSPE), SMIR = mean(SMIR_MSPE))
  return(list(result = result))
}
                                                                      
Settings <- expand.grid(n = c(100, 500, 1000),
                        SigmaTildeType = 1:3)

# Divide into 20 batches, each having 10 runs
id_artemis <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX")) # from 1-80
if(is.na(id_artemis)) id_artemis <- 21
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
save(Sim_onesigma, file = sprintf("simulation-results/RMIR/timevariant/results_sim_index%.2d.Rdata", id_artemis))


