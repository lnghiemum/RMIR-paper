### Inverse regression random effect_v2
### The code estimates the parameters using MCEM algorithm for isotropic Random effect 
require(Rfast)
require(pracma)
require(matrixNormal)

marginal_llh <- function(X, y, fy, C, Delta, Gamma0, sigma2, B=400, ...){
  n <- max(y$cluster)
  p <- ncol(X)
  Gamma0 <- as.matrix(Gamma0)
  d <- ncol(Gamma0)
  K <- diag(p) - tcrossprod(Gamma0)
  
  v <-  array(rnorm(B*p*d, sd = sqrt(sigma2)), dim = c(p, d, B))
  theta <- apply(v, 3, function(onev) K %*% onev, simplify = FALSE) 
  
  GeneratedGamma <- lapply(theta, expmap, Gamma0)
  marginal_llh_onecluster <- computeMarginalLLH(X, y, fy, Gamma0, C, Delta, GeneratedGamma, B)
  mean(marginal_llh_onecluster[is.finite(marginal_llh_onecluster)], na.rm=TRUE)
}

Estep_allclusters <- function(X, y, fy, Gamma0, C, Delta, sigma2, B = 400, ...){
  #browser()
  n <- max(y$cluster)
  ### sample B Gamma from the prior distribution
  p <- ncol(X)
  d <- ncol(Gamma0)
  pass = FALSE
  K <- diag(p) - tcrossprod(Gamma0)
  
  v <-  array(rnorm(B*p*d, sd = sqrt(sigma2)), dim = c(p, d, B))
  theta <- apply(v, 3, function(onev) K %*% onev, simplify = FALSE) 
  
  GeneratedGamma <- lapply(theta, expmap, Gamma0 = Gamma0)
  wb <- computeWeightsEstep(X, y, fy, Gamma0, C, Delta, GeneratedGamma, B)  
  return(list(wb = wb, Gamma = GeneratedGamma, Theta = theta))
}


updateSigma <- function(Theta_array, p, d, wb){
  # Theta_array: data from E-step
  tmp <- sapply(Theta_array, tcrossprod, simplify = FALSE)
  rowSumsTheta <- colSums(wb)
  S <- Reduce("+", Map("*", tmp, rowSumsTheta)) /(sum(wb)*d)
  estSigmastar <- S# MLESigmaSingular(S, p-d)
  return(estSigmastar)
}  


updatesigma2 <- function(Theta_array, Gamma0, wb){
  
  n <- nrow(wb)
  p <- nrow(Gamma0)
  d <- ncol(Gamma0)
  K <- diag(p) - tcrossprod(Gamma0)
  
  tmp <- sapply(Theta_array, tcrossprod, simplify = FALSE)
  tmp2 <- sapply(1:n, function(i) {
    S <- Reduce("+", Map("*", tmp, t(wb[i, ])))
    sum(diag(K %*% S))
  })
  hatsigma2 <- sum(tmp2)/(sum(wb)*(p-d)*d) 
  return(hatsigma2)
}  


MCEM_update_fixedGamma0 <- function(X, y, fy, InitEst, niter = 20, burnin = niter/2, 
                                 tol = 0.005, mc.cores = detectCores()/2, B = 400){
  estGamma0 <- InitEst$Gamma0
  estDelta <- InitEst$Delta
  estC <- InitEst$C
  estsigma2 <- InitEst$sigma2
  estGamma  <- InitEst$Gamma
  p <- ncol(X)
  d = ncol(estGamma0)
  
  logllk_cw <- rep(NA, niter)
  # compute loglikelihood at the initial value
  
  n <- max(y$cluster)
  logllk_cw[1] <-
    marginal_llh(X, y, fy, estC, estDelta,  estGamma0, estsigma2)
  iter = 1; conv = FALSE
  
  while(iter < niter & conv == FALSE){
    # E-step
    message("Iter ", iter)
    tmp <- Estep_allclusters(X, y, fy, estGamma0, estC, estDelta, estsigma2, B = B, mc.cores = mc.cores) 
    print("Estep done")
    wb <- tmp$wb
    Gamma_array <- tmp$Gamma
    Theta_array <- tmp$Theta
    #browser()
    estDelta <- updateDelta(X, y, fy, estC, wb, Gamma_array) 
    print("Update Delta done")
    estC <- updateC(X, y, fy,  estDelta, wb, Gamma_array)
    print("Update C done")

    estsigma2 <- updatesigma2(Theta_array, estGamma0, wb)
    print("Update Sigma done")
    
    ### monitoring the marginal log likelihood for convergence
    # Since the log likelihood is the function of Gamma_i, here average the theta and project it to get Gamma_i
    logllk_cw[iter + 1] <- marginal_llh(X, y, fy , estC, estDelta, estGamma0, estsigma2)
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
  K = diag(p) - tcrossprod(estGamma0)
  estSigma = estsigma2 * K
  return(list(estsigma2 = estsigma2,
              estSigma = estSigma, 
              estC = estC, 
              estGamma0 = estGamma0, 
              estDelta = estDelta, 
              iter = iter, llh = logllk_cw))
}






