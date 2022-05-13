################################################
library(parallel)
library(Rcpp)
library(tidyverse)
library(fastGHQuad)
source("data_generating_functions.R")
sourceCpp("fast_estimation_functions.cpp")
################################################
# Factors #
n <- 500
mnratio <- 3
################################################
# Data Generation #
Mu_YX <- c(2,1,1,1)
SigMat_YX <- ar1_cor(4, 0.9)
SigMat_YX[-1,-1] <- ar1_cor(3, 0.3)
SigMat_YX <- SigMat_YX+0.1*diag(4)
SigMat_YX[1,1] <- 1.44
Mu_Y_T <- 1.5
Sig_Y_T <- 1.5
Mu_Y_S <- Mu_YX[1]
Sig_Y_S <- sqrt(SigMat_YX[1,1])

trueBetaRho <- Compute_Rho_Parameters(Mu_Y_T, Sig_Y_T, Mu_YX, SigMat_YX)
yx_dist <- Compute_Y_Given_X(Mu_YX, SigMat_YX)
coef_y_x_s_true <- c(yx_dist$Beta0, yx_dist$Beta1)
var_y_x_s_true <- yx_dist$VarYX
sigma_y_x_s_true <- sqrt(var_y_x_s_true)
################################################
m <- mnratio*n
beta_rho <- trueBetaRho
B1 <- 10 # Monte-Carlo Sample Size
B2 <- 10 # Bootstrap (Perturbation Size)

gh_num <- 10
ghxw <- gaussHermiteData(gh_num)
xList <- ghxw$x
wList <- ghxw$w
################################################
set.seed(888)
data_list_mc <- mclapply(1:B1, function(x) {
  Generate_Dat(n, m, Mu_YX, SigMat_YX, Mu_Y_T, Sig_Y_T)
}, mc.cores = detectCores())
################################################