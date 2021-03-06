################################################
# Testing Rcpp functions
# (make sure they agree with the R versions)
################################################

library(Rcpp)
library(tidyverse)
library(fastGHQuad)
library(microbenchmark)
source("data_generating_functions.R")
source("estimation_functions.R")
sourceCpp("fast_estimation_functions.cpp")

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
m <- 2000

# Extra Data (or CV)

n_ext <- 1000
m_ext <- 1000

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
# Prepare Objects
################################################

theta = Mu_Y_T
beta_rho <- trueBetaRho
rho_pwr <- 2
x_mat_no_intercept <- tData
coef_y_x_s <- coef_y_x_s_true
sigma_y_x_s <- sqrt(var_y_x_s_true)
ghObj <- gaussHermiteData(10)
xList <- ghObj$x
wList <- ghObj$w

parameters <- list(y_vec = sData[,"Y"], mu = Mu_Y_S, sigma = Sig_Y_S, num_of_repl = 10,
                   xList = xList, wList= wList)

xMatAll <- rbind(sData[,-1], tData)
################################################
# Testing
################################################

E_S_RHO_X_CPP(beta_rho, rho_pwr, x_mat_no_intercept,
              coef_y_x_s, sigma_y_x_s, xList, wList) %>% sd
E_S_RHO_X(beta_rho, rho_pwr, x_mat_no_intercept,
          coef_y_x_s, sigma_y_x_s, 10) %>% sd
# microbenchmark(E_S_RHO_X_CPP(beta_rho, rho_pwr, x_mat_no_intercept,
#                              coef_y_x_s, sigma_y_x_s, xList, wList),
#                E_S_RHO_X(beta_rho, rho_pwr, x_mat_no_intercept,
#                          coef_y_x_s, sigma_y_x_s, 10))
################################################

ispar <- T
E_S_RHO_CPP(beta_rho, ispar, parameters)
E_S_RHO(beta_rho, ispar, parameters)

# microbenchmark(E_S_RHO_CPP(beta_rho, ispar, parameters),
#                E_S_RHO(beta_rho, ispar, parameters))

ispar <- F
E_S_RHO_CPP(beta_rho, ispar, parameters)
E_S_RHO(beta_rho, ispar, sData[,"Y"])

# microbenchmark(E_S_RHO_CPP(beta_rho, ispar, parameters),
#                E_S_RHO(beta_rho, ispar, sData[,"Y"]))
################################################

e_s_rho_x <- E_S_RHO_X_CPP(beta_rho, 1, x_mat_no_intercept,
                           coef_y_x_s, sigma_y_x_s, xList, wList)
e_s_rho2_x <- E_S_RHO_X_CPP(beta_rho, 2, x_mat_no_intercept,
                           coef_y_x_s, sigma_y_x_s, xList, wList)
ispar <- T
c_ps <- E_S_RHO_CPP(beta_rho, ispar, parameters)

COMPUTE_TAU_CPP(e_s_rho_x, e_s_rho2_x, c_ps, piVal) %>% sd
COMPUTE_TAU(e_s_rho_x, e_s_rho2_x, c_ps, piVal) %>% sd

# microbenchmark(COMPUTE_TAU_CPP(e_s_rho_x, e_s_rho2_x, c_ps, piVal),
#                COMPUTE_TAU(e_s_rho_x, e_s_rho2_x, c_ps, piVal))
################################################

E_T_TAU_CPP(e_s_rho_x, e_s_rho2_x, c_ps, piVal)
E_T_TAU(e_s_rho_x, e_s_rho2_x, c_ps, piVal)

# microbenchmark(E_T_TAU_CPP(e_s_rho_x, e_s_rho2_x, c_ps, piVal),
#                E_T_TAU(e_s_rho_x, e_s_rho2_x, c_ps, piVal))

################################################
ispar <- T
E_T_D_LOG_RHO_DIV_D_BETA_CPP(beta_rho, ispar, parameters, c_ps)
E_T_D_LOG_RHO_DIV_D_BETA(beta_rho, ispar, parameters, c_ps)
microbenchmark(E_T_D_LOG_RHO_DIV_D_BETA_CPP(beta_rho, ispar, parameters, c_ps),
               E_T_D_LOG_RHO_DIV_D_BETA(beta_rho, ispar, parameters, c_ps))

ispar <- F
E_T_D_LOG_RHO_DIV_D_BETA_CPP(beta_rho, ispar, parameters, c_ps)
E_T_D_LOG_RHO_DIV_D_BETA(beta_rho, ispar, sData[,"Y"], c_ps)
microbenchmark(E_T_D_LOG_RHO_DIV_D_BETA_CPP(beta_rho, ispar, parameters, c_ps),
               E_T_D_LOG_RHO_DIV_D_BETA(beta_rho, ispar, sData[,"Y"], c_ps))
################################################
E_T_D_LOG_RHO_DIV_D_BETA_X_CPP(beta_rho, x_mat_no_intercept,
                               coef_y_x_s, sigma_y_x_s, e_s_rho_x,
                               xList, wList) %>% sd
E_T_D_LOG_RHO_DIV_D_BETA_X(beta_rho, x_mat_no_intercept,
                           coef_y_x_s, sigma_y_x_s, e_s_rho_x, 10) %>% sd
microbenchmark(E_T_D_LOG_RHO_DIV_D_BETA_X_CPP(beta_rho, x_mat_no_intercept,
                                              coef_y_x_s, sigma_y_x_s, e_s_rho_x,
                                              xList, wList),
               E_T_D_LOG_RHO_DIV_D_BETA_X(beta_rho, x_mat_no_intercept,
                                          coef_y_x_s, sigma_y_x_s, e_s_rho_x, 10))
################################################
ispar <- T
Compute_S_CPP(beta_rho, x_mat_no_intercept, coef_y_x_s,
              sigma_y_x_s, e_s_rho_x, xList, wList, ispar, parameters, c_ps) %>% sd
Compute_S(beta_rho, x_mat_no_intercept, coef_y_x_s,
          sigma_y_x_s, e_s_rho_x, ispar, parameters, c_ps) %>% sd
microbenchmark(Compute_S_CPP(beta_rho, x_mat_no_intercept, coef_y_x_s,
                             sigma_y_x_s, e_s_rho_x, xList, wList, ispar, parameters, c_ps),
               Compute_S(beta_rho, x_mat_no_intercept, coef_y_x_s,
                         sigma_y_x_s, e_s_rho_x, ispar, parameters, c_ps))

ispar <- F
Compute_S_CPP(beta_rho, x_mat_no_intercept, coef_y_x_s,
              sigma_y_x_s, e_s_rho_x, xList, wList, ispar, parameters, c_ps) %>% sd
Compute_S(beta_rho, x_mat_no_intercept, coef_y_x_s,
          sigma_y_x_s, e_s_rho_x, ispar, sData[,"Y"], c_ps) %>% sd
microbenchmark(Compute_S_CPP(beta_rho, x_mat_no_intercept, coef_y_x_s,
                             sigma_y_x_s, e_s_rho_x, xList, wList, ispar, parameters, c_ps),
               Compute_S(beta_rho, x_mat_no_intercept, coef_y_x_s,
                         sigma_y_x_s, e_s_rho_x, ispar, sData[,"Y"], c_ps))
################################################
ispar <- T
ComputeEfficientScore_CPP(beta_rho, sData, tData, piVal, tData, coef_y_x_s,
                          sigma_y_x_s, ispar, parameters, xList, wList) %>% sd
ComputeEfficientScore(beta_rho, sData, tData, piVal, tData, coef_y_x_s,
                      sigma_y_x_s, ispar, parameters) %>% sd
microbenchmark(ComputeEfficientScore_CPP(beta_rho, sData, tData, piVal, tData, coef_y_x_s,
                                         sigma_y_x_s, ispar, parameters, xList, wList),
               ComputeEfficientScore(beta_rho, sData, tData, piVal, tData, coef_y_x_s,
                                     sigma_y_x_s, ispar, parameters))

ispar <- F
ComputeEfficientScore_CPP(beta_rho, sData, tData, piVal, tData, coef_y_x_s,
                          sigma_y_x_s, ispar, parameters, xList, wList) %>% sd
ComputeEfficientScore(beta_rho, sData, tData, piVal, tData, coef_y_x_s,
                      sigma_y_x_s, ispar, sData[,"Y"]) %>% sd
################################################
ispar <- T
EstimateBetaFunc_CPP(beta_rho, sData, tData, piVal, tData, coef_y_x_s,
                     sigma_y_x_s, ispar, parameters, xList, wList)
EstimateBetaFunc(beta_rho, sData, tData, piVal, tData, coef_y_x_s,
                 sigma_y_x_s, ispar, parameters)
microbenchmark(EstimateBetaFunc_CPP(beta_rho, sData, tData, piVal, tData, coef_y_x_s,
                                    sigma_y_x_s, ispar, parameters, xList, wList),
               EstimateBetaFunc(beta_rho, sData, tData, piVal, tData, coef_y_x_s,
                                sigma_y_x_s, ispar, parameters))

ispar <- F
EstimateBetaFunc_CPP(beta_rho, sData, tData, piVal, tData, coef_y_x_s,
                     sigma_y_x_s, ispar, parameters, xList, wList)
EstimateBetaFunc(beta_rho, sData, tData, piVal, tData, coef_y_x_s,
                 sigma_y_x_s, ispar, sData[,"Y"])
microbenchmark(EstimateBetaFunc_CPP(beta_rho, sData, tData, piVal, tData, coef_y_x_s,
                                    sigma_y_x_s, ispar, parameters, xList, wList),
               EstimateBetaFunc(beta_rho, sData, tData, piVal, tData, coef_y_x_s,
                                sigma_y_x_s, ispar, sData[,"Y"]))
################################################
ispar <- T
optim(beta_rho, EstimateBetaFunc_CPP, sData = sData, tData = tData, piVal = piVal,
      tDat_ext = tData, coef_y_x_s = coef_y_x_s, sigma_y_x_s = sigma_y_x_s, ispar = ispar,
      parameters = parameters, xList = xList, wList = wList) -> fop
optim(beta_rho, EstimateBetaFunc, sData = sData, tData = tData, piVal = piVal,
      tDat_ext = tData, coef_y_x_s = coef_y_x_s, sigma_y_x_s = sigma_y_x_s, ispar = ispar,
      parameters = parameters)
EstimateBetaVarFunc(fop$par, sData, tData, piVal, tData, coef_y_x_s, sigma_y_x_s, ispar, parameters)
EstimateBetaVarFunc_CPP(fop$par, sData, tData, piVal, tData, coef_y_x_s, sigma_y_x_s, ispar, parameters, xList, wList)
microbenchmark(EstimateBetaVarFunc(fop$par, sData, tData, piVal, tData, coef_y_x_s, sigma_y_x_s, ispar, parameters),
               EstimateBetaVarFunc_CPP(fop$par, sData, tData, piVal, tData, coef_y_x_s, sigma_y_x_s, ispar, parameters, xList, wList))

optim(beta_rho, EstimateBetaFunc_CPP, sData = sData, tData = tData, piVal = piVal,
      tDat_ext = tData, coef_y_x_s = coef_y_x_s, sigma_y_x_s = sigma_y_x_s, ispar = ispar,
      parameters = parameters, xList = xList, wList = wList, weights = TRUE) -> fop
fop$par

ComputeRandomizedWeightBootstrap(beta_rho, sData, tData, piVal, tData, coef_y_x_s, sigma_y_x_s, ispar,
                                 parameters, xList, wList)
################################################
e_t_tau <- E_T_TAU_CPP(e_s_rho_x, e_s_rho2_x, c_ps, piVal)
tau_x <- COMPUTE_TAU_CPP(e_s_rho_x, e_s_rho2_x, c_ps, piVal)
COMPUTE_B(tau_x, e_t_tau, c_ps, piVal) %>% sd
COMPUTE_B_CPP(tau_x, e_t_tau, c_ps, piVal) %>% sd
################################################
E_S_RHO2_PSI_X(beta_rho, x_mat_no_intercept, coef_y_x_s, sigma_y_x_s,
               xList, wList) %>% sd
E_S_RHO2_PSI_X_CPP(beta_rho, x_mat_no_intercept, coef_y_x_s, sigma_y_x_s,
                   xList, wList) %>% sd
################################################
e_s_rho_x_all <- E_S_RHO_X_CPP(beta_rho, 1, xMatAll, coef_y_x_s, sigma_y_x_s, xList, wList)
e_s_rho2_x_all <- E_S_RHO_X_CPP(beta_rho, 2, xMatAll, coef_y_x_s, sigma_y_x_s, xList, wList)

e_s_rho_x_t <- E_S_RHO_X_CPP(beta_rho, 1, tData, coef_y_x_s, sigma_y_x_s, xList, wList)
e_s_rho2_x_t <- E_S_RHO_X_CPP(beta_rho, 2, tData, coef_y_x_s, sigma_y_x_s, xList, wList)

tau_x_all <- COMPUTE_TAU_CPP(e_s_rho_x_all, e_s_rho2_x_all, c_ps, piVal)
tau_x_ext <- COMPUTE_TAU_CPP(e_s_rho_x_t, e_s_rho2_x_t, c_ps, piVal)

e_s_rho2_psi_x_all <- E_S_RHO2_PSI_X_CPP(beta_rho, xMatAll, coef_y_x_s, sigma_y_x_s, xList, wList)
e_s_rho2_psi_x_ext <- E_S_RHO2_PSI_X_CPP(beta_rho, tData, coef_y_x_s, sigma_y_x_s, xList, wList)

COMPUTE_A_CPP(beta_rho = beta_rho, e_t_tau = e_t_tau,
              tau_x_internal = tau_x_all, tau_x_external = tau_x_ext,
              e_s_rho2_psi_x_internal = e_s_rho2_psi_x_all, e_s_rho2_x_internal = e_s_rho2_x_all,
              e_s_rho2_psi_x_external = e_s_rho2_psi_x_ext, e_s_rho2_x_external = e_s_rho2_x_t) %>% sd
COMPUTE_A(beta_rho = beta_rho, e_t_tau = e_t_tau,
              tau_x_internal = tau_x_all, tau_x_external = tau_x_ext,
              e_s_rho2_psi_x_internal = e_s_rho2_psi_x_all, e_s_rho2_x_internal = e_s_rho2_x_all,
              e_s_rho2_psi_x_external = e_s_rho2_psi_x_ext, e_s_rho2_x_external = e_s_rho2_x_t) %>% sd
################################################

e_s_rho2_psi_x_source <- E_S_RHO2_PSI_X_CPP(beta_rho, sData[,-1], coef_y_x_s,
                                            sigma_y_x_s, xList, wList)
e_s_rho2_x_source <- E_S_RHO_X_CPP(beta_rho, 2, sData[,-1], coef_y_x_s,
                                   sigma_y_x_s, xList, wList)
e_s_rho1_x_source <- E_S_RHO_X_CPP(beta_rho, 1, sData[,-1], coef_y_x_s,
                                   sigma_y_x_s, xList, wList)
tau_for_x_source <- COMPUTE_TAU_CPP(e_s_rho1_x_source, e_s_rho2_x_source, c_ps, piVal)

COMPUTE_H_CPP(beta_rho, theta, c_ps, sData, e_t_tau, tau_x_ext, piVal,
              coef_y_x_s, sigma_y_x_s, e_s_rho2_psi_x_ext, e_s_rho2_x_t,
              e_s_rho2_psi_x_source, e_s_rho2_x_source, tau_for_x_source) %>% sd
COMPUTE_H(beta_rho, theta, c_ps, sData, e_t_tau, tau_x_ext, piVal, coef_y_x_s, sigma_y_x_s, e_s_rho2_psi_x_ext, e_s_rho2_x_external = e_s_rho2_x_t) %>% sd
################################################

GH_Materials <- list(mu = Mu_Y_S, sigma = Sig_Y_S,
                     xList = xList, wList = wList)
E_S_RHO_Y_CPP(beta_rho, T, GH_Materials, sData)
E_S_RHO_Y_CPP(beta_rho, F, GH_Materials, sData)

E_S_RHO_Y(beta_rho, 1, T, GH_Materials, sData)
E_S_RHO_Y(beta_rho, 2, T, GH_Materials, sData)
E_S_RHO_Y(beta_rho, 3, T, GH_Materials, sData)

################################################

ispar <- T
COMPUTE_1ST_PART_D_WRAPPER(theta, beta_rho, ispar, GH_Materials,
                           sData, c_ps, tau_x_ext, piVal, coef_y_x_s, sigma_y_x_s,
                           e_s_rho2_psi_x_ext, e_s_rho2_x_t)

e_s_rho2_psi_x_internal_source <- E_S_RHO2_PSI_X_CPP(beta_rho, sData[,-1], coef_y_x_s, sigma_y_x_s, xList, wList)
e_s_rho2_x_internal_source <- E_S_RHO_X_CPP(beta_rho, 2, sData[,-1], coef_y_x_s, sigma_y_x_s, xList, wList)
e_s_rho1_x_internal_source <- E_S_RHO_X_CPP(beta_rho, 1, sData[,-1], coef_y_x_s, sigma_y_x_s, xList, wList)
tau_for_x_source <- COMPUTE_TAU_CPP(e_s_rho1_x_internal_source, e_s_rho2_x_internal_source, c_ps, piVal)

COMPUTE_D_MULTIPLIER(beta_rho, ispar, GH_Materials, sData,
                     theta, piVal, c_ps, e_t_tau,
                     tau_x_ext, coef_y_x_s, sigma_y_x_s,
                     e_s_rho2_psi_x_ext, e_s_rho2_x_t,
                     e_s_rho2_psi_x_internal_source, e_s_rho2_x_internal_source, tau_for_x_source,
                     e_s_rho1_x_internal_source, xList, wList)

################################################

num_of_target <- nrow(tData)
yVec <- sData[,"Y"]
rhoValSource <- exp(c(cbind(yVec, yVec^2) %*% matrix(beta_rho, ncol = 1)))
c_ps <- E_S_RHO_CPP(beta_rho, ispar, GH_Materials)
e_s_rho_x_ext <- E_S_RHO_X_CPP(beta_rho, 1, tData, coef_y_x_s,
                               sigma_y_x_s, xList, wList)
e_s_rho2_x_ext <- E_S_RHO_X_CPP(beta_rho, 2, tData, coef_y_x_s,
                                sigma_y_x_s, xList, wList)
e_t_tau <- E_T_TAU_CPP(e_s_rho_x = e_s_rho_x_ext, e_s_rho2_x = e_s_rho2_x_ext, c_ps = c_ps, piVal = piVal)
tau_x_external <- COMPUTE_TAU_CPP(e_s_rho_x_ext, e_s_rho2_x_ext, c_ps, piVal)

xMatAll <- rbind(sData[,-1], tData)
e_s_rho_x_all <- E_S_RHO_X_CPP(beta_rho, 1, xMatAll, coef_y_x_s,
                               sigma_y_x_s, xList, wList)
e_s_rho2_x_all <- E_S_RHO_X_CPP(beta_rho, 2, xMatAll, coef_y_x_s,
                                sigma_y_x_s, xList, wList)
tau_x_internal_all <- COMPUTE_TAU_CPP(e_s_rho_x_all, e_s_rho2_x_all, c_ps, piVal)
e_s_rho2_psi_x_internal_all <- E_S_RHO2_PSI_X_CPP(beta_rho, xMatAll, coef_y_x_s, sigma_y_x_s, xList, wList)
e_s_rho2_psi_x_external <- E_S_RHO2_PSI_X_CPP(beta_rho, tData, coef_y_x_s, sigma_y_x_s, xList, wList)

MatInv <- EstimateBetaCovMat_CPP(beta_rho, sData, tData, piVal, tData,
                                 coef_y_x_s, sigma_y_x_s, ispar, GH_Materials,
                                 xList, wList)
e_s_rho2_psi_x_internal_source <- E_S_RHO2_PSI_X_CPP(beta_rho, sData[,-1], coef_y_x_s, sigma_y_x_s, xList, wList)
e_s_rho2_x_internal_source <- E_S_RHO_X_CPP(beta_rho, 2, sData[,-1], coef_y_x_s, sigma_y_x_s, xList, wList)
e_s_rho_x_internal_source <- E_S_RHO_X_CPP(beta_rho, 1, sData[,-1], coef_y_x_s, sigma_y_x_s, xList, wList)
tau_for_x_source <- COMPUTE_TAU_CPP(e_s_rho_x = e_s_rho_x_internal_source, e_s_rho2_x = e_s_rho2_x_internal_source, c_ps = c_ps, piVal = piVal)

SEff <- ComputeEfficientScore_CPP(beta_rho, sData, tData, piVal, tData,
                                  coef_y_x_s, sigma_y_x_s, ispar, GH_Materials,
                                  xList, wList)
COMPUTE_EFFICIENT_IF_FOR_THETA_CPP(theta, num_of_target, piVal,
                                   rhoValSource, c_ps, yVec,
                                   beta_rho, e_t_tau,
                                   tau_x_internal_all, tau_x_external,
                                   e_s_rho2_psi_x_internal_all,
                                   e_s_rho2_x_all, e_s_rho2_psi_x_external,
                                   e_s_rho2_x_ext, MatInv, ispar, GH_Materials,
                                   sData, coef_y_x_s, sigma_y_x_s,
                                   e_s_rho2_psi_x_internal_source,
                                   e_s_rho2_x_internal_source,
                                   tau_for_x_source, e_s_rho_x_internal_source,
                                   xList, wList, SEff) %>% sd
