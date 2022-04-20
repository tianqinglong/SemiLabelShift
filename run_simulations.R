#################################
library(parallel)
source("data_generating_functions.R")
source("estimation_functions.R")
#################################

# There are len(n_vec)*len(m_div_n_vec) combinations; each combination should run in a separate task

n_vec <- c(500, 600, 700, 800, 900, 1000)
m_div_n_vec <- c(1, 2, 3, 4, 5)

#################################

# Data generating setting (fixed for the setting)

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
#################################

# Monte Carlo Sample Size

B <- 1000

#################################

n <- n_vec[1]
m_div_n <- m_div_n_vec[1]
m <- m_div_n*n
beta_rho <- trueBetaRho

#################################

data_list_mc <- mclapply(1:B, function(x) {
  Generate_Dat(n, m, Mu_YX, SigMat_YX, Mu_Y_T, Sig_Y_T)
},
mc.cores = 12)

results_output <- mclapply(data_list_mc, function(dataList) {
  sData <- dataList$sDat
  tData <- dataList$tDat
  piVal <- n/(n+m)
  
  ## For E_s(\cdot): use empirical distribution of Y on the source
  ispar <- F
  parameters <- sData[,"Y"]
  ## For E_s(\cdot|x): use correctly specified model of f_s(y|x)
  fityx <- lm(Y~., data = as.data.frame(sData))
  coef_y_x_s <- coef(fityx)
  sigma_y_x_s <- sigma(fityx)
  
  estBeta <- optim(beta_rho, EstimateBetaFunc,
                   sData = sData, tData = tData, piVal = piVal, tDat_ext = tData,
                   coef_y_x_s = coef_y_x_s, sigma_y_x_s = sigma_y_x_s, ispar = ispar, parameters = parameters)
  betaHat <- estBeta$par
  sdVec <- EstimateBetaVarFunc(betaHat, sData, tData, piVal, tData, coef_y_x_s, sigma_y_x_s, ispar, parameters)
  cp1 <- (trueBetaRho[1] >= betaHat[1]-1.96*sdVec[1]) && (trueBetaRho[1] <= betaHat[1]+1.96*sdVec[1])
  cp2 <- (trueBetaRho[2] >= betaHat[2]-1.96*sdVec[2]) && (trueBetaRho[2] <= betaHat[2]+1.96*sdVec[2])
  
  return(list(Estimated = estBeta$par,
              Sd = sdVec,
              Cp1 = cp1,
              Cp2 = cp2))
},
mc.cores = 12)
# mean(sapply(results_output, function(out) {out$Cp1}))
# mean(sapply(results_output, function(out) {out$Cp2}))

saveRDS(results_output, file = paste("n",n,"ratio",m_div_n,"_correct.rds", sep = ""))
#################################

results_output <- mclapply(data_list_mc, function(dataList) {
  sData <- dataList$sDat
  tData <- dataList$tDat
  piVal <- n/(n+m)
  
  ## For E_s(\cdot): use empirical distribution of Y on the source
  ispar <- T
  parameters <- list(mu = Mu_Y_S, sigma = Sig_Y_S, num_of_repl = 18)
  ## For E_s(\cdot|x): use correctly specified model of f_s(y|x)
  coef_y_x_s <- coef_y_x_s_true
  sigma_y_x_s <- sqrt(var_y_x_s_true)
  
  estBeta <- optim(beta_rho, EstimateBetaFunc,
                   sData = sData, tData = tData, piVal = piVal, tDat_ext = tData,
                   coef_y_x_s = coef_y_x_s, sigma_y_x_s = sigma_y_x_s, ispar = ispar, parameters = parameters)
  betaHat <- estBeta$par
  sdVec <- EstimateBetaVarFunc(betaHat, sData, tData, piVal, tData, coef_y_x_s, sigma_y_x_s, ispar, parameters)
  cp1 <- (trueBetaRho[1] >= betaHat[1]-1.96*sdVec[1]) && (trueBetaRho[1] <= betaHat[1]+1.96*sdVec[1])
  cp2 <- (trueBetaRho[2] >= betaHat[2]-1.96*sdVec[2]) && (trueBetaRho[2] <= betaHat[2]+1.96*sdVec[2])
  
  return(list(Estimated = estBeta$par,
              Sd = sdVec,
              Cp1 = cp1,
              Cp2 = cp2))
},
mc.cores = 12)
#################################
saveRDS(results_output, file = paste("n",n,"_ratio",m_div_n,"_true.rds", sep = ""))