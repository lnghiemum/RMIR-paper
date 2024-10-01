### Inverse regression random effect_v2


marginal_llh <- function(X, y, fy, W, 
                         muX, Gamma0, C, muW,
                         beta, Delta, Sigma, B=400, ...){
  
  n <- max(y$cluster)
  p <- ncol(X)
  Gamma0 <- as.matrix(Gamma0)
  d <- ncol(Gamma0)
  K <- diag(p) - tcrossprod(Gamma0)
  ### Generate data from truncated MVN
  theta <- lapply(
    1:B, function(b) rMVNormC_eigen(d, rep(0, p), Sigma))
  GeneratedGamma <- lapply(theta, expmap, Gamma0)
  marginal_llh_onecluster <- computeMarginalLLH(X, y, fy, W, muX, Gamma0, C, muW, beta, Delta, GeneratedGamma, B)
  mean(marginal_llh_onecluster[is.finite(marginal_llh_onecluster)], na.rm=TRUE)
}

Estep_allclusters <- function(X, y, fy, W, 
                              muX, Gamma0, C, muW,
                              beta, Delta, Sigma, B = 400, ...){
  #browser()
  n <- max(y$cluster)
  ### sample B Gamma from the prior distribution
  p <- ncol(X)
  d <- ncol(Gamma0)
  pass = FALSE
  K <- diag(p) - tcrossprod(Gamma0)
  
  theta <- lapply(
    1:B, function(b) rMVNormC_eigen(d, rep(0, p), Sigma))
  GeneratedGamma <- lapply(theta, expmap, Gamma0 = Gamma0)
  wb <- computeWeightsEstep(X, y, fy, W, muX, 
                            Gamma0, C, muW, beta, Delta, GeneratedGamma, B)  
  
  return(list(wb = wb, Gamma = GeneratedGamma, Theta = theta))
}
updateSigma <- function(Theta_array, p, d, wb){
  # Theta_array: data from E-step
  tmp <- sapply(Theta_array, tcrossprod, simplify = FALSE)
  rowSumsTheta <- colSums(wb)
  
  S <- Reduce("+", Map("*", tmp, rowSumsTheta)) /(sum(wb)*d)

  return(S)
}  

MCEM_update_noGamma0 <- function(X, y, fy, W, InitEst, niter = 20, burnin = niter/2, 
                                   tol = 0.005, mc.cores = detectCores()/2, B = 400){
  #browser()
  estGamma0 <- InitEst$Gamma0
  estDelta <- InitEst$Delta
  estbeta <- InitEst$beta
  estmuX <- InitEst$muX
  estmuW <- InitEst$muW
  estC <- InitEst$C
  estSigma <- InitEst$Sigma
  estGamma  <- InitEst$Gamma
  p <- ncol(X)
  d = ncol(estGamma0)
  
  logllk_cw <- rep(NA, niter)
  # compute loglikelihood at the initial value
  
  n <- max(y$cluster)
  logllk_cw[1] <-
   marginal_llh(X, y, fy, W,
                estmuX, estGamma0, estC, estmuW,
                estbeta, estDelta, estSigma)
  iter = 1; conv = FALSE
  
  while(iter < niter & conv == FALSE){
    # E-step
    message("Iter ", iter)
    tmp <- Estep_allclusters(X, y, fy, W, 
                             estmuX, estGamma0, estC, estmuW,
                             estbeta, estDelta, estSigma, B = B, mc.cores = mc.cores) 
    print("Estep done")
    wb <- tmp$wb
    Gamma_array <- tmp$Gamma
    Theta_array <- tmp$Theta
    # browser()
    estmuX <- updatemuX(X, y, fy, W, estbeta, estC, estmuW, wb, Gamma_array)
    print("Update muX done")
    estmuW <- updatemuW(X, y, fy, W, estmuX, estbeta, estC, wb, Gamma_array)
    print("Update muW done")
  
    estbeta <- updatebeta(X, y, fy, W, estmuX, estC, estmuW,  wb, Gamma_array)
    print("Update beta done")

    
    estDelta <- updateDelta(X, y, fy, W, estmuX, estC, estmuW, estbeta, wb, Gamma_array) 
    print("Update Delta done")
    estC <- updateC(X, y, fy, W, estmuX, estbeta, estmuW, estDelta, wb, Gamma_array)
    print("Update C done")
    
    estSigma <- updateSigma(Theta_array, p, d, wb)
    print("Update Sigma done")
    
    ### monitoring the marginal log likelihood for convergence
    # Since the log likelihood is the function of Gamma_i, here average the theta and project it to get Gamma_i
    logllk_cw[iter + 1] <-  marginal_llh(X, y, fy, W, 
                                         estmuX, estGamma0, estC, estmuW,
                                         estbeta, estDelta, estSigma, mc.cores = mc.cores)
    # Monitoring convergence through the change in logllk
    diff_logllk <- logllk_cw[iter + 1] - logllk_cw[iter]  
    message("diff in log llh: ", diff_logllk)
    
    if (abs(diff_logllk) < tol & iter > burnin){
      conv = TRUE
    }
    iter = iter  + 1
  }
  d <- ncol(estGamma0)
  log_llh <- logllk_cw[iter]
  num_par <- 
    n*p + # mu_i for each cluster
    (p-d)*d + # Gamma0  
    d*ncol(fy) + # beta
    p*(p+1)/2 +  # Delta 
    (p-d)*(p-d+1)/2 
  AIC <- -2*n*log_llh + 2*num_par
  BIC <- -2*n*log_llh + log(nrow(y))*num_par
  
  return(list(estSigma = estSigma,
              estmuX = estmuX,
              estmuW = estmuW,
              estC = estC, 
              estGamma0 = estGamma0, 
              estDelta = estDelta, 
              estbeta = estbeta,
              iter = iter, llh = logllk_cw, AIC = AIC, BIC = BIC))
}
  


