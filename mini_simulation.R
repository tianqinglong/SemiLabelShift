################################################
library(parallel)
source("data_generating_functions.R")
source("estimation_functions.R")
################################################

################################################
# Hyper-parameters
################################################

# Parameters for data generating
Mu_YX <- c(2,1,1,1)
SigMat_YX <- ar1_cor(4, 0.9)
SigMat_YX[-1,-1] <- ar1_cor(3, 0.3)
SigMat_YX <- SigMat_YX+0.1*diag(4)
SigMat_YX[1,1] <- 1.44
Mu_Y_T <- 1.5
Sig_Y_T <- 1.5

Mu_Y_S <- Mu_YX[1]
Sig_Y_S <- sqrt(SigMat_YX[1,1])

## True Beta
trueBetaRho <- Compute_Rho_Parameters(Mu_Y_T, Sig_Y_T, Mu_YX, SigMat_YX)

# Sample size

n <- 500
m <- 500

# Extra Data (or CV)

n_ext <- 10000
m_ext <- 10000

data_list_extra <- Generate_Dat(n_ext, m_ext, Mu_YX, SigMat_YX, Mu_Y_T, Sig_Y_T)

# True Distribution
yx_dist <- Compute_Y_Given_X(Mu_YX, SigMat_YX)
coef_y_x_s_true <- c(yx_dist$Beta0, yx_dist$Beta1)
var_y_x_s_true <- yx_dist$VarYX

################################################
# Data Generation
################################################

dataList <- Generate_Dat(n, m, Mu_YX, SigMat_YX, Mu_Y_T, Sig_Y_T)
sData <- dataList$sDat
tData <- dataList$tDat
piVal <- n/(n+m)

# External Data for experimental use only
sData_ext <- data_list_extra$sDat
tDat_ext <- data_list_extra$tDat

################################################
# Estimation and Inference (single dataset)
################################################

beta_rho <- trueBetaRho

# Use True for both E_s(\cdot) and E_s(\cdot|x)
## For E_s(\cdot): use empirical distribution
ispar <- T
parameters <- list(mu = Mu_Y_S, sigma = Sig_Y_S, num_of_repl = 18)
## For E_s(\cdot|x)
coef_y_x_s <- coef_y_x_s_true
sigma_y_x_s <- sqrt(c(var_y_x_s_true))

estBeta <- optim(beta_rho, EstimateBetaFunc,
                 sData = sData, tData = tData, piVal = piVal, tDat_ext = tData,
                 coef_y_x_s = coef_y_x_s, sigma_y_x_s = sigma_y_x_s, ispar = ispar, parameters = parameters)
(betaHat <- estBeta$par)
EstimateBetaVarFunc(betaHat, sData, tData, piVal, tData, coef_y_x_s, sigma_y_x_s, ispar, parameters)

# Use estimated version of E_s(\cdot) and E_s(\cdot|x)
## For E_s(\cdot)
ispar <- F
parameters <- sData[,"Y"]
## For E_s(\cdot|x)
fityx <- lm(Y~., data = as.data.frame(sData))
coef_y_x_s <- coef(fityx)
sigma_y_x_s <- sigma(fityx)

estBeta <- optim(beta_rho, EstimateBetaFunc,
                 sData = sData, tData = tData, piVal = piVal, tDat_ext = tData,
                 coef_y_x_s = coef_y_x_s, sigma_y_x_s = sigma_y_x_s, ispar = ispar, parameters = parameters)
(betaHat <- estBeta$par)
EstimateBetaVarFunc(betaHat, sData, tData, piVal, tData, coef_y_x_s, sigma_y_x_s, ispar, parameters)

################################################
# Estimation and Inference: Monte Carlo
################################################
B <- 500
data_list_mc <- mclapply(1:B, function(x) {
  Generate_Dat(n, m, Mu_YX, SigMat_YX, Mu_Y_T, Sig_Y_T)
},
mc.cores = 12)

# Use True for both E_s(\cdot) and E_s(\cdot|x)
## For E_s(\cdot): use empirical distribution
ispar <- T
parameters <- list(mu = Mu_Y_S, sigma = Sig_Y_S, num_of_repl = 18)
## For E_s(\cdot|x)
coef_y_x_s <- coef_y_x_s_true
sigma_y_x_s <- sqrt(c(var_y_x_s_true))

results_output <- mclapply(data_list_mc, function(dataList) {
  sData <- dataList$sDat
  tData <- dataList$tDat
  piVal <- n/(n+m)
  
  estBeta <- optim(beta_rho, EstimateBetaFunc,
                   sData = sData, tData = tData, piVal = piVal, tDat_ext = tData,
                   coef_y_x_s = coef_y_x_s, sigma_y_x_s = sigma_y_x_s, ispar = ispar, parameters = parameters)
  betaHat <- estBeta$par
  sdVec <- EstimateBetaVarFunc(betaHat, sData, tData, piVal, tData, coef_y_x_s, sigma_y_x_s, ispar, parameters)
  cp1 <- (trueBetaRho[1] >= betaHat[1]-1.96*sdVec[1]) && (trueBetaRho[1] <= betaHat[1]+1.96*sdVec[1])
  cp2 <- (trueBetaRho[2] >= betaHat[2]-1.96*sdVec[2]) && (trueBetaRho[2] <= betaHat[2]+1.96*sdVec[2])
  
  return(list(Estimated = fopt$par,
              Sd = sdVec,
              Cp1 = cp1,
              Cp2 = cp2))
},
mc.cores = 12)
mean(sapply(results_output, function(out) {out$Cp1}))
mean(sapply(results_output, function(out) {out$Cp2}))

# Use estimated version of E_s(\cdot) and E_s(\cdot|x)
## For E_s(\cdot)
ispar <- F
parameters <- sData[,"Y"]
## For E_s(\cdot|x)
fityx <- lm(Y~., data = as.data.frame(sData))
coef_y_x_s <- coef(fityx)
sigma_y_x_s <- sigma(fityx)

results_output <- mclapply(data_list_mc, function(dataList) {
  sData <- dataList$sDat
  tData <- dataList$tDat
  piVal <- n/(n+m)
  
  estBeta <- optim(beta_rho, EstimateBetaFunc,
                   sData = sData, tData = tData, piVal = piVal, tDat_ext = tData,
                   coef_y_x_s = coef_y_x_s, sigma_y_x_s = sigma_y_x_s, ispar = ispar, parameters = parameters)
  betaHat <- estBeta$par
  sdVec <- EstimateBetaVarFunc(betaHat, sData, tData, piVal, tData, coef_y_x_s, sigma_y_x_s, ispar, parameters)
  cp1 <- (trueBetaRho[1] >= betaHat[1]-1.96*sdVec[1]) && (trueBetaRho[1] <= betaHat[1]+1.96*sdVec[1])
  cp2 <- (trueBetaRho[2] >= betaHat[2]-1.96*sdVec[2]) && (trueBetaRho[2] <= betaHat[2]+1.96*sdVec[2])
  
  return(list(Estimated = fopt$par,
              Sd = sdVec,
              Cp1 = cp1,
              Cp2 = cp2))
},
mc.cores = 12)
mean(sapply(results_output, function(out) {out$Cp1}))
mean(sapply(results_output, function(out) {out$Cp2}))

