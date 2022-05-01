################################################
library(parallel)
library(Rcpp)
library(tidyverse)
library(fastGHQuad)
source("data_generating_functions.R")
source("estimation_functions.R")
sourceCpp("fast_estimation_functions.cpp")
################################################
# Factors
n <- 3000
mnratio <- 1
################################################
# Data Generating Parameters
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

# fityx <- lm(Y~., data = as.data.frame(sData))
# coef_y_x_s_hat <- coef(fityx)
# sigma_y_x_s_hat <- sigma(fityx)
# 
# Mu_Y_S_hat <- mean(sData[,"Y"])
# Sig_Y_S_hat <- sd(sData[,"Y"])
################################################
m <- mnratio*n
beta_rho <- trueBetaRho
B1 <- 720 # Monte-Carlo Sample Size
# B2 <- 10 # Bootstrap (Perturbation Size)

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
results_output <- mclapply(data_list_mc, function(dataList) {
  sData <- dataList$sDat
  tData <- dataList$tDat
  piVal <- n/(n+m)
  ispar <- T
  parameters <- list(y_vec = sData[,"Y"], mu = Mu_Y_S, sigma = Sig_Y_S, xList = xList, wList= wList)
  
  betaHat <- optim(beta_rho, EstimateBetaFunc_CPP, sData = sData, tData = tData, piVal = piVal,
                   tDat_ext = tData, coef_y_x_s = coef_y_x_s_true, sigma_y_x_s = sigma_y_x_s_true, ispar = ispar,
                   parameters = parameters, xList = xList, wList = wList)
  betaHat <- betaHat$par
  betaSd1 <- EstimateBetaVarCenterFunc(betaHat, sData, tData, piVal, tData, coef_y_x_s_true, sigma_y_x_s_true, ispar, parameters, xList = xList, wList = wList)
  # betaSd1 <- EstimateBetaVarFunc_CPP(betaHat, sData, tData, piVal, tData, coef_y_x_s_true, sigma_y_x_s_true, ispar, parameters, xList, wList)
  # betaBoot <- ComputeRandomizedWeightBootstrap(beta_rho, sData, tData, piVal, tData, coef_y_x_s_true, sigma_y_x_s_true, ispar, parameters, xList, wList, B2 = B2)
  
  CI1 <- matrix(betaHat, nrow = 2, ncol = 2, byrow = T)
  Sd1 <- matrix(betaSd1, nrow = 2, ncol = 2, byrow = T)
  Sd1[1,] <- -1.96*Sd1[1,]
  Sd1[2,] <- 1.96*Sd1[2,]
  CI1 <- CI1+Sd1
  
  # CI2 <- matrix(nrow = 2, ncol = 2)
  # CI2 <- apply(betaBoot, MARGIN = 2, quantile, probs = c(0.025, 0.975))
  
  CP1 <- c(trueBetaRho[1] >= CI1[1,1] & trueBetaRho[1] <= CI1[2,1],
           trueBetaRho[2] >= CI1[1,2] & trueBetaRho[2] <= CI1[2,2])
  # CP2 <- c(trueBetaRho[1] >= CI2[1,1] & trueBetaRho[1] <= CI2[2,1],
  #          trueBetaRho[2] >= CI2[1,2] & trueBetaRho[2] <= CI2[2,2])
  
  return(list(
    TrueBeta = trueBetaRho,
    BetaHat = betaHat,
    Sd = betaSd1,
    CI1 = CI1,
    CP1 = CP1
    # CI2 = CI2,
    # CP2 = CP2
  ))
}, mc.cores = detectCores())
################################################
sapply(results_output, function(out) {
  out$CP1
}) %>% rowMeans()

sapply(results_output, function(out) {
  out$Sd
}) %>% rowMeans()

sapply(results_output, function(out) {
  out$BetaHat
}) %>% apply(MARGIN = 1, sd)

sapply(results_output, function(out) {
  out$BetaHat-out$TrueBeta
}) %>% rowMeans()
