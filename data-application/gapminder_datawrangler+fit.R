library(tidyverse)
library(lme4)
library(gapminder)
library(ldr)
library(Rcpp)
library(RcppArmadillo)

source("RMIR/global_Ising_with_covariates.R")
source("RMIR/ancillary_functions.R")
source("RMIR/RMIR_MCEM_code_allcontinuous+qbinary.R")
sourceCpp("RMIR/RMIR_MCEM_update_allcontinuous+qbinary.cpp")
sourceCpp("RMIR/globalIsing.cpp")

##############################
### Data wrangling ###########
##############################
dat_income <- read_delim("data/ddf--datapoints--income_per_person_long_series--by--geo--time.csv",
                         col_types = c("c", "c", "d"), delim = ",") %>%
  mutate(logincome = log(income_per_person_long_series))

dat_infant_mortality <- read_delim("data/ddf--datapoints--infant_mortality_rate_per_1000_births--by--geo--time.csv",
                                   col_types = c("c", "c", "d"), delim = ",") 

dat_sexratio <- read_delim("data/ddf--datapoints--sex_ratio_all_age_groups--by--geo--time.csv",
                           col_types = c("c", "c", "d"), delim = ",") 

dat_lifeExpectancy_female <- read_delim("data/ddf--datapoints--life_expectancy_female--by--geo--time.csv",
                                        col_types = c("c", "c", "d"), delim = ",")

dat_co2_emission <- read_delim("data/ddf--datapoints--consumption_emissions_tonnes_per_person--by--geo--time.csv",
                               col_types = c("c", "c", "d"), delim = ",")

dat_childperwoman <- read_delim("data/ddf--datapoints--children_per_woman_total_fertility--by--geo--time.csv",
                                col_types = c("c", "c", "d"), delim = ",")

dat_gini <- read_delim("data/ddf--datapoints--gapminder_gini--by--geo--time.csv",
                       col_types = c("c", "c", "d"), delim = ",")

# Merging dataset
df_merge <- left_join(dat_lifeExpectancy_female, dat_income, by = c("geo", "time")) 
df_merge <- left_join(df_merge, dat_sexratio, by = c("geo", "time")) %>% 
  left_join(., dat_infant_mortality, by = c("geo", "time")) %>%
  left_join(., dat_co2_emission, by = c("geo", "time")) %>%
  left_join(., dat_childperwoman, by = c("geo", "time")) %>%
  left_join(., dat_gini, by = c("geo", "time")) %>%
  filter(!if_any(everything(), is.na))

### time-invariant binary covariates
country_data <- read_csv("data/ddf--entities--geo--country.csv") %>% dplyr::select(country, un_sdg_ldc, west_and_rest) %>%
  rename(geo = country)

df_final <- left_join(df_merge, country_data, by = "geo")
#############################################
# summary
df_final %>% group_by(geo) %>% summarize(n())

#############################################
###### Extracting X, Y, W components #######
#############################################

X <- df_final[, c(5:10)] %>% mutate(across(1:6, ~ scale(.x, center = FALSE))) # column 4 is income, column 5 is log income, use log income
p <- ncol(X)
Wdf <- df_final[, c(1, 11:12)] %>% unique() %>%
  mutate(
    # Convert character to numeric
    west_and_rest = ifelse(west_and_rest == "rest", 0, 1),
    un_sdg_ldc = ifelse(un_sdg_ldc == "un_not_least_developed", 0, 1))
W <- as.matrix(Wdf[, -1]) 
q <- ncol(W); kq <- q*(q-1)/2; kk <- q*(q+1)/2
y_tmp <- df_final[, 1:3] 
# create a cluster variable for y
levels_country <- unique(y_tmp$geo)
n <- length(levels_country)
cluster_level <- 1:n
df_cluster_geo <- data.frame(geo = levels_country, cluster = cluster_level)
Y <- left_join(y_tmp, df_cluster_geo, by = "geo")
names(Y)[3] <- "Y"

countryinfo <- data.frame(geo = levels_country) %>% mutate(iso_alpha = toupper(geo)) %>% 
  left_join(., gapminder::country_codes) %>% select(-iso_num) 
countryinfo[countryinfo$iso_alpha == "LAO", ]$country <- "Laos"
countryinfo[countryinfo$iso_alpha == "KGZ", ]$country <- "Kyrgyzstan"
country = unique(countryinfo$country)

##############################
##### Model fits #############
##############################

## Binary parts Fit a global Ising model on W|Ybar
FY <- data.frame(Y$Y)
FYbar <- Reduce("rbind", tapply(FY, Y$cluster, colMeans, simplify = TRUE))
fitWandY <- fit_Ising(W, FYbar, control = list(tol = 1e-5, maxIter = 20))
estA <- fitWandY$estGamma
estC2 <- fitWandY$estC


clusterSize <- c(table(Y$geo))
Wrep <- W[rep(1:n, clusterSize),]
          
# Fitting models on continuous part X | W, Y
# SAIC method for choosing d

d_cand = 0:3
separatePFCfit <- lapply(1:n, function(jj){
  try({
    Xi  <- X[Y$cluster == jj, ]
    yi <- Y$Y[Y$cluster==jj]
    fyi <- fy[Y$cluster == jj, ]
    pfc.fit <- pfc(Xi, yi, fyi, numdir = max(d_cand),
                   structure = "unstr", numdir.test = TRUE)
    # Getting the full log lik (for d > 0)
    pfc.fit$loglik2 <-  mapply(function(Gamma, Delta, beta){
      Xic <- scale(Xi, center = TRUE, scale = FALSE)
      lkh <- Rfast::dmvnorm(Xic - fyi %*% t(beta) %*% t(Gamma), mu = rep(0, p), sigma = Delta, logged = TRUE)
    }, pfc.fit$Gammahat, pfc.fit$Deltahat, pfc.fit$Betahat) %>% colSums()
    pfc.fit$loglik2 <- c(pfc.fit$loglik[1], pfc.fit$loglik2)
    return(pfc.fit)
  }, silent = T)
})
aic_sepPFC <- sapply(separatePFCfit, function(j){
  if (class(j) == "try-error") return(rep(NA, 4))
  else return(-2*j$loglik2 + 2*j$numpar)
})
SAIC <- rowSums(aic_sepPFC, na.rm = TRUE) # select d = 2

d <- d_cand[which.min(SAIC)]

GMIR <- list()
fy = ldr::bf(Y$Y, case = "poly", degree = 3, scale = TRUE)
fitXandWY <- RRRR::RRR(y = X, x = fy, z = Wrep, 
                       mu = TRUE,
                       r = d)
estbeta <- fitXandWY$D 
estGamma0 <- qr.Q(qr(fitXandWY$A))
estC1 <- qr.R(qr(fitXandWY$A)) %*% t(fitXandWY$B)
estSigma <- fitXandWY$Sigma

GMIR$estDelta <- estSigma
GMIR$estbeta <- fitXandWY$D 
GMIR$estGamma0 <- estGamma0
GMIR$Sigmastar <- NA 
GMIR$estGamma <- NA
GMIR$estBasisContinuousPart <- solve(GMIR$estDelta) %*% projection(estGamma0)

### SMIR fit
SMIR <- list()
n <- max(Y$cluster)
X <- as.matrix(X)
sep_fit <- lapply(1:n, function(ii){
  Xi <- X[Y$cluster == ii, ]
  Yi <- Y[Y$cluster == ii,]$Y
  fyi <- fy[Y$cluster == ii,, drop = FALSE]
  ni <- nrow(Xi)
  fit <- try( RRRR::RRR(y = Xi, x = fyi, z = NULL, r = d))
  if(class(fit) == "try-error"){
    fit <- fitXandWY
  }
  sep_fit_Gamma <-  qr.Q(qr(fit$A))
  sep_fit_C1 <- qr.R(qr(fit$A)) %*% t(fit$B)
  sep_fit_Delta <- fit$Sigma
  sep_fit_basis <- solve(sep_fit_Delta, sep_fit_Gamma)
  return(list(sep_fit_Gamma, sep_fit_Delta, sep_fit_basis))  
})

# sep_fit for Gamma0
sep_fit_Gamma <- lapply(sep_fit, "[[", 1)
SMIR$estGamma0 <- frechetMeanGr(sep_fit_Gamma)
SMIR$estV <- lapply(sep_fit_Gamma, invexpmap, Gamma0 = SMIR$estGamma0)
SMIR$Sigmastar <- Reduce("+", lapply(SMIR$estV, tcrossprod))/length(SMIR$estV)
SMIR$estDelta <- Reduce("+", lapply(sep_fit, "[[", 2))/n
SMIR$estMeanCS <- solve(SMIR$estDelta) %*% SMIR$estGamma0
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
estV <- lapply(sep_fit_Gamma, invexpmap, Gamma0 = InitEst$Gamma0)
InitEst$Sigmastar <- Reduce("+", lapply(estV, tcrossprod))/length(estV)

fy <- cbind(1, fy)
InitEst$C <- matrix(rnorm(d * ncol(fy)), nrow = d)
RMIRfit <- MCEM_update_noGamma0(X, Y, fy, Wrep, InitEst, niter = 80, burnin = 10, 
                                tol = 0.01, mc.cores = detectCores()/2, B = 800)
r <- ncol(RMIRfit$estC)
RMIRfit$estBasisContinuousPart <-
  qr.Q(qr(solve(RMIRfit$estDelta, RMIRfit$estGamma0)))
RMIRfit$meanCSbasis <- 
  rbind(solve(RMIRfit$estDelta, RMIRfit$estGamma0),
        -t(RMIRfit$estbeta) %*% solve(RMIRfit$estDelta, RMIRfit$estGamma0))

#### Prediction of cluster-specific spaces
# Final E-step
doEstepfinal <- Estep_allclusters(X, Y, fy, Wrep, 
                                  RMIRfit$estmuX, RMIRfit$estGamma0, RMIRfit$estC, 
                                  RMIRfit$estmuW,
                                  RMIRfit$estbeta, RMIRfit$estDelta, RMIRfit$estSigmastar, B = 2000) 

wb <- doEstepfinal$wb
Gamma_array <- doEstepfinal$Gamma
Theta_array <- doEstepfinal$Theta # on tangent space at estGamma0

posteriorGamma <- lapply(1:n, function(ii){
  tmp <- Reduce("+", Map("*",  Theta_array, wb[ii, ]))
  gi <- expmap(tmp, RMIRfit$estGamma0)
})
### Cluster specific dimension components
# all interactions of Wrep
vechWW <- data.frame(W) %>%
  rename(W1 = un_sdg_ldc, W2 = west_and_rest) %>%
  mutate(W12 = W1 * W2) %>% as.matrix()

RMIRfit$estThetai <- mapply(function(Gammai, ni){
  col1 <- solve(RMIRfit$estDelta, Gammai)
  col2 <- -t(RMIRfit$estbeta) %*% solve(RMIRfit$estDelta, Gammai)
  rbind(col1, col2)
}, posteriorGamma, clusterSize, SIMPLIFY = F)

clusterSpecificDR <- lapply(1:n, function(j){
  Xi <- X[Y$cluster == j, ]
  estGamma <- qr.Q(qr(RMIRfit$estThetai[[j]]))
  # Make the first components of each cluster-specific basis to be positive
  if (estGamma[1, 1] > 0)  estGamma[, 1] <- -1 *  estGamma[, 1] 
  if (estGamma[1, 2] < 0)  estGamma[, 2] <- -1 *  estGamma[, 2]
  sweep(Xi %*% estGamma[1:p, ], 2, as.matrix(vechWW[j, 1:2, drop = F]) %*% t(RMIRfit$estbeta) %*% estGamma[1:p,]) 
})

meanDR <- X %*% RMIRfit$estBasisContinuousPart +
  Wrep %*% t(RMIRfit$estbeta) %*% RMIRfit$estBasisContinuousPart

data_application <- list(data  = df_final,
                         X = X, 
                         Y = Y, 
                         W = W,
                         vechWW = vechWW,
                         Wrep = Wrep,
                         RMIR_fit = RMIRfit,
                         RMIR_pred_ranef = RMIRfit$estThetai,
                         RMIR_DR = clusterSpecificDR,
                         RMIR_meanDR = meanDR,
                         A = estA)
data_application$SMIR <- SMIR 
data_application$GMIR <- GMIR
data_application$country_info <- countryinfo
data_application$finalEstep = doEstepfinal

save(data_application, file = "gapminder_applications_withbinarycovariates.Rdata")

