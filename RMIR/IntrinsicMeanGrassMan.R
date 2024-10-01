### function to compute intrinsic mean of points on Grassmann manifold

### The algorithm is based on Chakraborty (2015) 

# Geodesic function Gamma (p.4231)

computeGeodesic <- function(X, Y, t){
  # X and Y are two semi-orthogonal bases of the same dimension
  # t is a scalar
  n <- nrow(X); p <- ncol(X)
  projX <- X %*% solve(crossprod(X), t(X))
  M <- (diag(n) - projX) %*% Y %*% solve(crossprod(X, Y))
  svd.M <- svd(M)
  
  if (p>1){
    Theta = atan(diag(svd.M$d))
  }else{
    Theta = atan(svd.M$d)
  }
  A <- X %*% svd.M$v %*% cos(Theta*t) + svd.M$u %*% sin(Theta*t)   
  # finding orthogonal basis for A
  pracma::gramSchmidt(A)$Q
}
# test data
function(){
  n <- 5; p <- 2
  X <- matrix(rnorm(n*p), n, p)
  Y <- matrix(rnorm(n*p), n, p)
  computeGeodesic(X, Y, t=1)  
}

# compute the IGA

frechetMeanGr <- function(X){
  # X is a list of n matrices, each having dimension n \times p
  M <- X[[1]]
  n <- length(X)
  for (k in 2:n){
    M <- computeGeodesic(M, X[[k]], 1/(k))
  }
  return(M)
}

function(){
  n <- 200; p <- 4; d <- 2
  Xlist <- replicate(n, matrix(0.5*rnorm(p*d), p, d), simplify = F)
  frechetMeanGr(Xlist)    
}