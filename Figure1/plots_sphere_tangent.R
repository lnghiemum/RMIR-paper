library(RiemBase)
library(rgl)

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

frechetVar <- function(mydata){
  rgrassmann <- RiemBase::riemfactory(mydata, "grassmann")
  frechetmean <- RiemBase::rbase.mean(rgrassmann)
  
  crossprod(frechetmean$x)
  riemmean <- RiemBase::riemfactory(list(frechetmean$x), "grassmann")
  frechetVar <- mean(RiemBase::rbase.pdist2(rgrassmann, riemmean))
  frechetVar
}


library(Rfast)
library(rgl)
# Gamma0 = (0, 0, 1), so any point on tangent plane is (a, b, 0) where (a,b) from N_2
Gamma0 <- matrix(c(0, 0, 1), nrow = 3, ncol = 1)
K <- diag(3) - tcrossprod(Gamma0)
dat <- Rfast::rmvnorm(100, rep(0, 3), 0.3*diag(3), seed = 200)
theta <- dat %*% K
# Exponential mapping points to the sphere
Gamma <- apply(theta, 1, expmap, Gamma0 = Gamma0, simplify = TRUE)

### Visualization
open3d()
shapelist3d(octahedron3d(), x = theta[, 1], y = theta[, 2], z = theta[, 3]+1, color = "green", size = 0.03)
planes3d(a = 0, b = 0, c = 1, d = -1, color = "lightgray", alpha = 0.2)

rgl.spheres(x = 0, y=0, z = 0, radius = 1, phi = 90, color = "lightblue")
points3d(0, 0, 1, col = "red", size = 15)
rgl.spheres(x = Gamma[1,], y = Gamma[2,], z = Gamma[3, ], r = 0.03, color = "darkgreen")
close3d()

#### Computation to include Frechet variance
#### Create Points on the tangent plane
set.seed(20)
### Generating a big sample size
n <- 10^6
V <- Rfast::rmvnorm(n, rep(0, 3), 0.3*diag(3), seed = 1000)
theta <- V %*% K
# Exponential mapping points to the sphere
Gamma <- apply(theta, 1, expmap, Gamma0 = Gamma0, simplify = FALSE)
case1frechetVar <- frechetVar(Gamma)


# Separate two pictures separately
open3d()
bg3d(color = "white")
shapelist3d(octahedron3d(), x = theta[, 1], y = theta[, 2], z = theta[, 3]+1, color = "green", size = 0.06)
planes3d(a = 0, b = 0, c = 1, d = -1, color = "lightgray", alpha = 0.2)
snapshot3d("TangentSpace1.png", webshot = FALSE,
           width = 500, height = 500)
close3d()

open3d()
bg3d(color = "white")
rgl.spheres(x = 0, y=0, z = 0, radius = 1, phi = 90, color = "lightblue")
points3d(0, 0, 1, col = "red", size = 15)
rgl.spheres(x = Gamma[1,], y = Gamma[2,], z = Gamma[3, ], r = 0.03, color = "darkgreen")
um <- par3d()$userMatrix
view3d(userMatrix = um)
snapshot3d("Sphere1.png", webshot = FALSE, 
           width = 500, height = 500)
close3d()

### Second situation when the matrix Sigma is non-diagonal
###

p <- 3
set.seed(20)
#### Visualization
Sigma2 <- matrix(0.2, 3, 3)
diag(Sigma2) <- 0.5
Gamma0 <- matrix(c(0, 0, 1), nrow = 3, ncol = 1)
theta <- dat %*% K
# Exponential mapping points to the sphere
Gamma <- apply(theta, 1, expmap, Gamma0 = Gamma0, simplify = TRUE)
open3d()
bg3d(color = "white")
shapelist3d(octahedron3d(), x = theta[, 1], y = theta[, 2], z = theta[, 3] + 1, color = "green", size = 0.03)
planes3d(a = 0, b = 0, c = 1, d = -1, color = "lightgray", alpha = 0.2)
snapshot3d("TangentSpace2.png", webshot = FALSE,
           width = 500, height = 500)
close3d()

open3d()
rgl.spheres(x = 0, y=0, z = 0, radius = 1, phi = 90, color = "lightblue")
points3d(0, 0, 1, col = "red", size = 15)
rgl.spheres(x = Gamma[1,], y = Gamma[2,], z = Gamma[3, ], r = 0.03, color = "darkgreen")

um <- par3d()$userMatrix

view3d(userMatrix = um)
snapshot3d("Sphere2.png", webshot = FALSE, 
           width = 500, height = 500)
close3d()


n <- 10^6
V <- Rfast::rmvnorm(n, rep(0, 3), Sigma2, seed = 1000)
K <- diag(3) - tcrossprod(Gamma0)
# Create Points on the tangent plane
theta <- V %*% K
# Exponential mapping points to the sphere
Gamma <- apply(theta, 1, expmap, Gamma0 = Gamma0, simplify = FALSE)
case2frechetvar = frechetVar(Gamma)

# Separate two pictures separately
open3d()
bg3d(color = "white")
shapelist3d(octahedron3d(), x = x, y = y, z = z+1, color = "green", size = 0.06)
planes3d(a = 0, b = 0, c = 1, d = -1, color = "lightgray", alpha = 0.2)
snapshot3d("TangentSpace2.png", webshot = FALSE,
           width = 500, height = 500)
close3d()

open3d()
bg3d(color = "white")
rgl.spheres(x = 0, y=0, z = 0, radius = 1, phi = 90, color = "lightblue")
points3d(0, 0, 1, col = "red", size = 15)
rgl.spheres(x = Gamma[1,], y = Gamma[2,], z = Gamma[3, ], r = 0.03, color = "darkgreen")
um <- par3d()$userMatrix

view3d(userMatrix = um)
snapshot3d("Sphere2.png", webshot = FALSE, 
           width = 500, height = 500)
close3d()


