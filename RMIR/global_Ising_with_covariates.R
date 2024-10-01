require(huge)
require(dplyr)
require(IsingSampler)
require(pbmcapply) 
require(Rcpp)
require(RcppArmadillo)
### Global Ising model with covariate

fit_Ising <- function(X, fy, d = 1, 
                      control = list(tol = 1e-4, maxIter = 20)){
  # Initialization
  q <- ncol(X)
  kk <- q*(q+1)/2
  r <- ncol(fy)
  # id_grid to keep track of index 
  id_grid <- rbind(
    cbind(Var1 = 1:q, Var2 = 1:q),
    expand.grid(1:q, 1:q) %>% filter(Var1 < Var2) %>%
      arrange(Var1)
  )
  grid = as.matrix(id_grid - 1)
  
  estC <- matrix(0, nrow = r, ncol = d)
  esttau0 <- matrix(0, nrow = kk , ncol = 1) # each dimension 
  estGamma <- matrix(rnorm(kk*d), ncol = d, nrow = kk)
  estGamma <- qr.Q(qr(estGamma))
  logllh <- rep(NA, control$maxIter)
  logllh[1] <- pseudologLLH(X, fy, esttau0, estGamma, estC, grid)
  # Update parameters until convergence
  conv = F; iter = 1
  currentTau0 <- esttau0
  currentGamma <- estGamma
  currentC <- estC
  
  while(!conv & iter < control$maxIter){
    #print(iter)
    iter = iter + 1
    currentTau0 <- globalupdateTau0(X, fy, currentTau0, currentGamma, currentC, grid)
    currentC <- globalupdateC(X, fy, currentTau0, currentGamma, currentC, grid)
    currentGamma <- globalupdateGamma(X, fy, currentTau0, currentGamma, currentC, grid)
    currentGamma <- qr.Q(qr(currentGamma))
    logllh[iter] <- pseudologLLH(X, fy, currentTau0, currentGamma, currentC, grid)
    diff = logllh[iter] - logllh[iter-1]
    if (abs(diff) < control$tol){
      conv = TRUE
    }
  }
  return(list(estGamma = currentGamma, 
              estTau0 = currentTau0,
              estC = currentC,
              logllh = logllh[iter]))
}
# test function
function(){
  n <- 1000; q <- 4; d <- 1
  # for each observation, both the threshold and the graph are functions of y
  Y <- rnorm(n, mean = 0.3)
  kk <- q*(q+1)/2 # q diagonal terms and q(q-1)/2 off-diagonal
  fy <- as.matrix(Y) # n \times r
  r <- ncol(fy)
  A <- matrix(c(rep(1, kk-2), 10, 10), ncol = 1)
  #  matrix(rnorm(kk), ncol = 1)
  A <- A/sqrt(sum(A^2))
  C <- matrix(2)
  ACt <- A %*% t(C) # dimension kk \times r
  ACf <- fy %*% t(ACt)
  X <- matrix(nrow = n, ncol = q)
  for(ii in 1:n){
    thresholds <- ACf[ii, 1:q]
    graphs <- matrix(0, nrow = q, ncol = q)
    # arr.ind for upper triangular elements
    upper_tri_index <- which(upper.tri(graphs), arr.ind = T)
    lower_tri_index <- upper_tri_index[, c(2, 1)]
    graphs[upper_tri_index] <- ACf[ii, -c(1:q)]
    graphs[lower_tri_index] <- ACf[ii, -c(1:q)]
    X[ii, ] <- IsingSampler(n=1, graphs, thresholds)
  }
  fy = cbind(Y, Y^2) 
  globalIsing_fit = fit_Ising(X, fy)
  estGamma = globalIsing_fit[[1]]
  proj1 <- projection(estGamma)
  trueproj <- projection(A)
  norm(trueproj - proj1, "F")
}
  



