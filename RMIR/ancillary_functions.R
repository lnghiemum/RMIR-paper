### some supporting functions for singular normal
projection <- function(alpha){
  alpha <- as.matrix(alpha)
  A <- t(alpha) %*% alpha
  if (rcond(A)<1e-3){
    return(NA)
  }else{
    alpha %*% solve(A, t(alpha))  
  }
}

trunc_eigen_sqrt <- function(Sigma, inv){
  es <- eigen(Sigma)
  
  # Small negative eigenvalues can occur just due to numerical error and should
  # be set back to zero. 
  es$values[abs(es$values) < 1e-12] <- 0
  
  # If any eigen values are large negative throw error (something wrong)
  if (any(es$values < -1e-12)) stop("Non-trivial negative eigenvalues present")
  
  # calculate square root and reveal rank (k)
  k <- sum(es$values > 0)
  if (!inv){
    L <- es$vectors %*% diag(sqrt(es$values))  
  } else if (inv) {
    L <- es$vectors %*% diag(1/sqrt(es$values))  
  }
  return(L[,1:k, drop=F])
}

#' Covariance parameterization (Eigen decomposition)
#' @param n number of samples to draw
#' @param mu p-vector mean
#' @param Sigma covariance matrix (p x p)
#' @return matrix of dimension p x n of samples
rMVNormC_eigen <- function(n, mu, Sigma){
  p <- length(mu)
  L <- trunc_eigen_sqrt(Sigma, inv=FALSE)
  k <- ncol(L)
  Z <- matrix(rnorm(k*n), k, n)
  X <- L%*%Z
  X <- sweep(X, 1, mu, FUN=`+`)
}


#### solve function to ensure the Delta mat is always invertible in the MCEM algorithm
solveIf <- function(A, b = diag(ncol(A))){
  Ainv <- try(solve(A, b))
  if(class(Ainv)[1] == "try-error"){
    p <- ncol(A)
    eigenA <- eigen(A, only.values = TRUE)$values
    A <- A + 0.001*max(eigenA)*diag(p)
    Ainv <- solve(A, b)
  }
  Ainv
}

function(){
  d <- 2; m = 8; n = 300; p = 4
  A <- matrix(rnorm(p^2), p, p, byrow=TRUE)
  gs <- pracma::gramSchmidt(A)
  Gamma0 <- gs$Q[, 1:d]  
  K <- diag(p) - tcrossprod(Gamma0)
  
  #v <- replicate(n, matrixNormal::rmatnorm(M = pracma::zeros(p, d), 
  #                                                 U= sigma2*diag(p), V = diag(d)))
  v <- array(rnorm(n * p *d, sd = sqrt(sigma2)), dim = c(p, d, n)) 
  theta <- apply(v, 3, function(onev) K %*% onev, simplify = F)
  Gamma <- sapply(theta, expmap, Gamma0, simplify = F)
  
  datcluster <- sapply(1:n, function(i){
    X <- Rfast::rmvnorm(n = m, mu = rep(0, p), sigma = diag(p))
    y <- (X %*% Gamma[[i]][, 1])^3/2 + X %*% Gamma[[i]][, 2] +  0.1 * rnorm(m) 
    y <- data.frame(cluster = i, time = 1:m, yresp=y)
    return(list(X=X, Y=y))
  })
  X <- Reduce(rbind, datcluster[1,])
  y <- Reduce(rbind, datcluster[2,])
  fy <- cbind(y$yresp, y$yresp^2, y$yresp^3)
  beta <- matrix(0, nrow = d, ncol = ncol(fy))
  #Delta <- diag(p)
}


### Exponential map function from Theta_i (on tangent space) to Gamma_i (on manifold)
expmap <- function(theta, Gamma0){
  # both theta and gamma0 should be a p \times d matrix
  p <- nrow(Gamma0); d = ncol(Gamma0)
  K <- diag(p) - tcrossprod(Gamma0)
  qrdecom <- qr(K %*% theta) 
  Q <- qr.Q(qrdecom)
  R <- qr.R(qrdecom)
  svd.R <- svd(R)
  if (d > 1){
    M <- svd.R$v %*% diag(cos(svd.R$d)) %*% t(svd.R$v)
    NE <- svd.R$u %*% diag(sin(svd.R$d)) %*% t(svd.R$v)
  } else{
    M <- svd.R$v %*% cos(svd.R$d) %*% t(svd.R$v)
    NE <- svd.R$u %*% sin(svd.R$d) %*% t(svd.R$v)
  }
  Utilde <- Gamma0 %*% M + Q %*% NE
  return(Utilde)
}
### inverse exponential map function from gamma (a point on the manifold) to theta (a point on tangent space)
invexpmap <- function(Gamma, Gamma0){
  Gamma0 <- as.matrix(Gamma0)
  d <- ncol(Gamma0)
  Gamma <- as.matrix(Gamma)
  p <- nrow(Gamma0)
  K <- diag(p) - tcrossprod(Gamma0)
  M <- crossprod(Gamma0, Gamma)
  P <- K %*% Gamma %*% solve(M)
  svd.P <- svd(P)
  if (d > 1){
    svd.P$u %*% diag(atan(svd.P$d)) %*% t(svd.P$v)
  } else{
    svd.P$u %*% atan(svd.P$d) %*% t(svd.P$v)
  }
}

MLESigmaSingular <- function(S, r){
  ### Calculate the MLE of singular covariance matrix that has rank r
  eigenS <- eigen(S)
  Shat  <- eigenS$vectors[, 1:r] %*% diag(eigenS$values[1:r]) %*% t(eigenS$vectors[, 1:r])
  return(Shat) 
}  
