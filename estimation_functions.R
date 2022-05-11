library(fastGHQuad)
# This function computes $E_s{\rho^{1 or 2}(y)|x}$ for the internal x, coef_y_x_s and sigma_y must be estimated using external data
E_S_RHO_X <- function(beta_rho, rho_pwr, x_mat_no_intercept, coef_y_x_s, sigma_y_x_s,
                      num_of_replications = 10) {
  coef_y_x_s <- matrix(coef_y_x_s, ncol = 1)
  beta_rho <- matrix(beta_rho, ncol = 1)
  
  x_mat <- cbind(1, x_mat_no_intercept)
  mean_y_x_s <- c(x_mat %*% coef_y_x_s)
  coef_y_x_s <- matrix(coef_y_x_s, ncol = 1)
  
  x_mat <- cbind(1, x_mat_no_intercept)
  mean_y_x_s <- c(x_mat %*% coef_y_x_s)
  
  num_of_x <- nrow(x_mat_no_intercept)
  e_s_rho_x <- numeric(num_of_x)
  
  ghList <- gaussHermiteData(num_of_replications)
  xList <- ghList$x
  wList <- ghList$w
  
  for (i in 1:num_of_x) {
    y_plus_error <- sqrt(2)*c(sigma_y_x_s)*xList+c(mean_y_x_s[i])
    
    yMat <- cbind(y_plus_error, y_plus_error^2)
    rho_values <- c(yMat %*% beta_rho)
    
    e_s_rho_x[i] <- sum(wList*exp(rho_pwr*rho_values))/sqrt(pi)
  }
  
  return(e_s_rho_x)
}

# This function computes $E_s(\rho(Y))$ using different approaches
# ispar = T: we are going to compute using parametric model (either true or estimated); parameters contains mu and sigma and number of points is GH
# ispar = F: we use the empirical distribution; parameters are sample of y_s
E_S_RHO <- function(beta_rho, ispar, parameters) {
  if (ispar) {
    Mu_Y_S <- parameters$mu
    Sig_Y_S <- parameters$sigma
    num_of_repl <- parameters$num_of_repl
    
    ghList <- gaussHermiteData(num_of_repl)
    xList <- ghList$x
    wList <- ghList$w
    
    beta_rho <- matrix(beta_rho, ncol = 1)
    y_s_vector <- sqrt(2)*Sig_Y_S*xList+Mu_Y_S
    y_mat <- cbind(y_s_vector, y_s_vector^2)
    
    rho_values <- c(y_mat %*% beta_rho)
    e_s_rho <- sum(wList*exp(rho_values))/sqrt(pi)
  }
  else {
    y_vec <- matrix(parameters, ncol = 1)
    y_mat <- cbind(y_vec, y_vec^2)
    
    rho_values <- c(y_mat %*% beta_rho)
    e_s_rho <- mean(exp(rho_values))
  }
  return(e_s_rho)
}

# This function computes $\tau(\x)$ for the internal/external x, depending on the choice of e_s_rho_x and e_s_rho2_x
COMPUTE_TAU <- function(e_s_rho_x, e_s_rho2_x, c_ps, piVal) {
  e_t_rho_x <- e_s_rho2_x/e_s_rho_x
  tmp <- e_t_rho_x/c_ps/piVal
  tauX <- tmp/(tmp+1/(1-piVal))
  
  return(tauX)
}

# This function computes $E_t(\tau)$, x_t_mat_external.
# Here e_s_rho_x and e_s_rho2_x are related to x_t_mat_external, not the target sample
E_T_TAU <- function(e_s_rho_x, e_s_rho2_x, c_ps, piVal) {
  tau_vec <- COMPUTE_TAU(e_s_rho_x, e_s_rho2_x, c_ps, piVal)
  out <- mean(tau_vec)
  return(out)
}

# This function computes $E_t(d\log\Rho/d\beta)$, y_s_external is external
E_T_D_LOG_RHO_DIV_D_BETA <- function(beta_rho, ispar, parameters, c_ps) {
  beta_rho <- matrix(beta_rho, ncol = 1)
  if (ispar) {
    Mu_Y_S <- parameters$mu
    Sig_Y_S <- parameters$sigma
    num_of_rep <- parameters$num_of_repl
    
    ghList <- gaussHermiteData(num_of_rep)
    xList <- ghList$x
    wList <- ghList$w
    y_s_external <- sqrt(2)*Sig_Y_S*xList+Mu_Y_S
    
    yMat <- cbind(y_s_external, y_s_external^2)
    rhoVecSrc <- c(exp(yMat %*% beta_rho))
    # \beta_1: that is y
    outBeta1 <- sum(wList*rhoVecSrc*y_s_external)/c_ps/sqrt(pi)
    # \beta_2: that is y^2
    outBeta2 <- sum(wList*rhoVecSrc*(y_s_external^2))/c_ps/sqrt(pi)
  }
  else
  {
    y_s_external <- matrix(parameters, ncol = 1)
    yMat <- cbind(y_s_external, y_s_external^2)
    rhoVecSrc <- c(exp(yMat %*% beta_rho))
    # \beta_1: that is y
    outBeta1 <- mean(rhoVecSrc*y_s_external)/c_ps
    # \beta_2: that is y^2
    outBeta2 <- mean(rhoVecSrc*(y_s_external^2))/c_ps
  }
  
  return(c(outBeta1, outBeta2))
}

# This function computes $\E_t(d\log\Rho/d\beta|x)$, x might be internal/external
E_T_D_LOG_RHO_DIV_D_BETA_X <- function(beta_rho, x_mat_no_intercept, coef_y_x_s, sigma_y_x_s, e_s_rho_x, num_of_rep = 10) {
  x_mat <- cbind(1, x_mat_no_intercept)
  y_x_mean <- c(x_mat %*% coef_y_x_s)
  beta_rho <- matrix(beta_rho, ncol = 1)
  
  num_of_x <- nrow(x_mat)
  out_mat <- matrix(nrow = num_of_x, ncol = 2)

  ghList <- gaussHermiteData(num_of_rep)
  xList <- ghList$x
  wList <- ghList$w
    
  for (i in 1:num_of_x) {
    y_plus_error <- sqrt(2)*c(sigma_y_x_s)*xList+c(y_x_mean[i])
    y_tmp_mat <- cbind(y_plus_error, y_plus_error^2)
    rhoVecSrc <- exp(y_tmp_mat %*% beta_rho)
    
    tmp1 <- sum(wList*rhoVecSrc*y_plus_error)/sqrt(pi)
    tmp2 <- sum(wList*rhoVecSrc*y_plus_error^2)/sqrt(pi)
    
    out_mat[i,] <- c(tmp1, tmp2)
  }
  e_s_rho_x_mat <- matrix(e_s_rho_x, nrow = num_of_x, ncol = length(c(beta_rho)))
  out_mat <- out_mat/e_s_rho_x_mat
  
  return(out_mat)
}

# This function computes S(\x), for both internal/external x
Compute_S <- function(beta_rho, x_mat_no_intercept, coef_y_x_s, sigma_y_x_s, e_s_rho_x, ispar, parameters, c_ps) {
  e_t_d_log_rho_div_d_beta_x <- E_T_D_LOG_RHO_DIV_D_BETA_X(beta_rho, x_mat_no_intercept, coef_y_x_s, sigma_y_x_s, e_s_rho_x)
  e_t_d_log_rho_div_d_beta <- E_T_D_LOG_RHO_DIV_D_BETA(beta_rho, ispar, parameters, c_ps)
  e_t_d_log_rho_div_d_beta <- matrix(e_t_d_log_rho_div_d_beta, ncol = 2, nrow = nrow(x_mat_no_intercept), byrow = T)
  
  return(e_t_d_log_rho_div_d_beta_x-e_t_d_log_rho_div_d_beta)
}

# This function computes the efficient score function for the internal data
ComputeEfficientScore <- function(beta_rho, sData, tData, piVal, tDat_ext, coef_y_x_s, sigma_y_x_s, ispar, parameters) {
  num_of_source <- nrow(sData)
  num_of_target <- nrow(tData)
  num_of_total <- num_of_source+num_of_target
  
  beta_rho <- matrix(beta_rho, ncol = 1)
  y_internal <- sData[,"Y"]
  y_internal <- matrix(y_internal, ncol = 1)
  
  y_internal_mat <- cbind(y_internal, y_internal^2)
  rho_internal <- c(exp(y_internal_mat %*% beta_rho))
  
  c_ps <- E_S_RHO(beta_rho, ispar, parameters)
  
  xMatAll <- rbind(sData[,-1], tData)
  
  # Find the multiplier for the S_eff
  ## For the source sample
  multiplier1 <- rho_internal/c_ps/piVal
  ## For the target sample
  multiplier2 <- rep(-1/(1-piVal), num_of_target)
  ## Concatenate
  multiplier <- matrix(c(multiplier1, multiplier2), ncol = 2, nrow = num_of_total)
  
  # Find the first part of b_1: -(1-\pi)*(1-\tau(x)) for all the internal x
  e_s_rho2_x <- E_S_RHO_X(beta_rho, 2, xMatAll, coef_y_x_s, sigma_y_x_s)
  e_s_rho_x <- E_S_RHO_X(beta_rho, 1, xMatAll, coef_y_x_s, sigma_y_x_s)
  
  tau_x_internal <- COMPUTE_TAU(e_s_rho_x, e_s_rho2_x, c_ps, piVal)
  tau_x_internal <- matrix(c(tau_x_internal), ncol = 2, nrow = length(tau_x_internal), byrow = F)
  
  # Find the second part of b_1
  ## First find S(x), where x is internal
  s_x_internal <- Compute_S(beta_rho, xMatAll, coef_y_x_s, sigma_y_x_s, e_s_rho_x, ispar, parameters, c_ps)
  ## For E_t(\tau\S) and E_t(\tau) (using tDat_ext)
  e_s_rho_x_ext <- E_S_RHO_X(beta_rho, 1, tDat_ext, coef_y_x_s, sigma_y_x_s)
  e_s_rho2_x_ext <- E_S_RHO_X(beta_rho, 2, tDat_ext, coef_y_x_s, sigma_y_x_s)
  tau_x_external <- COMPUTE_TAU(e_s_rho_x_ext, e_s_rho2_x_ext, c_ps, piVal)
  s_x_external <- Compute_S(beta_rho, tDat_ext, coef_y_x_s, sigma_y_x_s, e_s_rho_x_ext, ispar, parameters, c_ps)
  
  e_t_tau <- mean(tau_x_external)
  e_t_tau_s <- colMeans(matrix(tau_x_external, ncol = 2, nrow = length(tau_x_external), byrow = F)*s_x_external)
  e_t_tau_s <- matrix(e_t_tau_s, ncol = 2, nrow = num_of_total, byrow = T)
  
  b_2nd_part <- s_x_internal-e_t_tau_s/(e_t_tau-1)
  
  b1_x <- -(1-piVal)*(1-tau_x_internal)*b_2nd_part
  SEffMat <- multiplier*b1_x
  
  return(SEffMat)
}

# This function is used for estimating Beta
EstimateBetaFunc <- function(beta_rho, sData, tData, piVal, tDat_ext, coef_y_x_s, sigma_y_x_s, ispar, parameters) {
  SEffMat <- ComputeEfficientScore(beta_rho, sData, tData, piVal, tDat_ext, coef_y_x_s, sigma_y_x_s, ispar, parameters)
  return(sum(colMeans(SEffMat)^2))
}

# This function is for the inference of Beta
EstimateBetaVarFunc <- function(beta_rho, sData, tData, piVal, tDat_ext, coef_y_x_s, sigma_y_x_s, ispar, parameters) {
  SEffMat <- ComputeEfficientScore(beta_rho, sData, tData, piVal, tDat_ext, coef_y_x_s, sigma_y_x_s, ispar, parameters)
  
  num_of_total <- nrow(sData)+nrow(tData)
  out <- matrix(0, ncol = 2, nrow = 2)
  for (i in 1:num_of_total) {
    tmp <- matrix(SEffMat[i,], ncol = 1)
    out <- out + tmp %*% t(tmp)
  }
  out <- out/num_of_total
  out <- solve(out)/num_of_total
  
  return(sqrt(diag(out)))
}

# This function repeat the pertubations
ComputeRandomizedWeightBootstrap <- function(beta_rho, sData, tData, piVal, tDat_ext,
                                             coef_y_x_s, sigma_y_x_s, ispar, parameters, xList, wList, B2 = 1000) {
  i <- 1
  outMat <- matrix(nrow = B2, ncol = 2)
  
  while(i <= B2)
  {
    optim(beta_rho, EstimateBetaFunc_CPP, sData = sData, tData = tData, piVal = piVal,
          tDat_ext = tDat_ext, coef_y_x_s = coef_y_x_s, sigma_y_x_s = sigma_y_x_s, ispar = ispar,
          parameters = parameters, xList = xList, wList = wList, weights = TRUE) -> fop
    fop$par -> outMat[i,]
    i <- i+1
  }
  
  return(outMat)
}

#################################
# Functions for estimating \psi #
#################################

COMPUTE_B <- function(tau_all, e_t_tau, c_ps, piVal) {
  out <- (e_t_tau-tau_all)/(1-e_t_tau)
  
  return(out)
}

E_S_RHO_Y <- function(beta_rho, pwr, ispar, parameters, sData)
{
  beta_rho <- matrix(beta_rho, ncol = 1)
  if(ispar)
  {
    xList <- parameters$xList
    wList <- parameters$wList
    
    mu_y_s <- parameters$Mu
    sigma_y_s <- parameters$Sigma
    
    yVec <- sqrt(2)*sigma_y_s*xList+mu_y_s
    yMat <- cbind(yVec, yVec^2)
    
    e_s_rho_y <- sum(wList*exp(c(yMat%*%beta_rho))*(yVec^pwr))/sqrt(pi)
    
    return(mean(e_s_rho_y))
  }
  
  yVec <- sData[,"Y"]
  yMat <- cbind(yVec, yVec^2)
  e_s_rho_y <- exp(c(yMat%*%beta_rho))*(yVec^pwr)
  
  return(mean(e_s_rho_y))
}

E_S_RHO2_PSI_X <- function(beta_rho, xmat_no_intercept, coef_y_x_s, sigma_y_x_s, xList, wList) {
  coef_y_x_s <- matrix(coef_y_x_s, ncol = 1)
  x_mat <- cbind(1, xmat_no_intercept)
  mean_y_x_s <- c(x_mat %*%  coef_y_x_s)
  beta_rho <- matrix(beta_rho, ncol = 1)
  
  num_of_x <- nrow(xmat_no_intercept)
  outMat <- numeric(num_of_x)
  for (i in 1:num_of_x) {
    y_plus_error <- sqrt(2)*c(sigma_y_x_s)*xList+c(mean_y_x_s)[i]
    yMat <- cbind(y_plus_error, y_plus_error^2)
    rho_values <- c(yMat %*% beta_rho)
    outMat[i] <- sum(wList*exp(2*rho_values)*y_plus_error)/sqrt(pi)
  }
  
  return(outMat)
}

COMPUTE_A <- function(beta_rho, e_t_tau,
                      tau_x_internal, tau_x_external,
                      e_s_rho2_psi_x_internal, e_s_rho2_x_internal,
                      e_s_rho2_psi_x_external, e_s_rho2_x_external) {
  tmp <- mean(tau_x_external*e_s_rho2_psi_x_external/e_s_rho2_x_external)
  aVec <- 1/(1-e_t_tau)*tmp-tau_x_internal*(1/(1-e_t_tau)*tmp+e_s_rho2_psi_x_internal/e_s_rho2_x_internal)
  
  return(aVec)
}

COMPUTE_H <- function(beta_rho, theta, c_ps, sData, e_t_tau, tau_x_external, piVal,
                      coef_y_x_s, sigma_y_x_s, e_s_rho2_psi_x_external, e_s_rho2_x_external)
{
  beta_rho <- matrix(beta_rho, ncol = 1)
  yVec <- sData[,"Y"]
  xMat <- sData[,-1]
  
  yMat <- cbind(yVec, yVec^2)
  rho_val <- exp(c(yMat%*%beta_rho))
  
  e_s_rho2_x_internal <- E_S_RHO_X(beta_rho, 2, xMat, coef_y_x_s, sigma_y_x_s)
  e_s_rho_x_internal <- E_S_RHO_X(beta_rho, 1, xMat, coef_y_x_s, sigma_y_x_s)
  tau_x_internal <- COMPUTE_TAU(e_s_rho_x_internal, e_s_rho2_x_internal, c_ps, piVal)
  
  e_s_rho2_psi_x_internal <- E_S_RHO2_PSI_X(beta_rho, xMat, coef_y_x_s, sigma_y_x_s, xList, wList)
  COMPUTE_A(beta_rho, e_t_tau, tau_x_internal, tau_x_external, e_s_rho2_psi_x_internal, e_s_rho2_x_internal,
            e_s_rho2_psi_x_external, e_s_rho2_x_external) -> aVec
  COMPUTE_B(tau_x_internal, e_t_tau, c_ps, piVal) -> bVec
  outVec <- 1/piVal*rho_val/c_ps*(yVec-theta+aVec-bVec*theta)
  
  return(outVec)
}

COMPUTE_1ST_PART_D <- function(theta, e_s_rho_y, e_s_rho_y2, e_s_rho_y3, c_ps, rho_val_s, h_val_s, s_val_s, piVal)
{
  first_part <- matrix(c(e_s_rho_y2, e_s_rho_y3), nrow = 1)/c_ps
  second_part <- -matrix(c(e_s_rho_y, e_s_rho_y2), nrow = 1)/c_ps*theta
  
  third_part <- matrix(rho_val_s*h_val_s, ncol = 2, nrow = length(rho_val_s), byrow = F)*
    s_val_s
  third_part <- -(1-piVal)*colMeans(third_part)/c_ps
  
  return(first_part+second_part+third_part)
}

COMPUTE_1ST_PART_D_WRAPPER <- function(theta, beta_rho, ispar, parameters, sData,
                                       c_ps, tau_x_external, piVal, coef_y_x_s, sigma_y_x_s,
                                       e_s_rho2_psi_x_external, e_s_rho2_x_external) {
  e_s_rho_y <- E_S_RHO_Y(beta_rho, 1, ispar, parameters, sData)
  e_s_rho_y2 <- E_S_RHO_Y(beta_rho, 2, ispar, parameters, sData)
  e_s_rho_y3 <- E_S_RHO_Y(beta_rho, 3, ispar, parameters, sData)
  
  yVec <- sData[,"Y"]
  yMat <- cbind(yVec, yVec^2)
  
  beta_rho <- matrix(beta_rho, ncol = 1)
  rho_val_s <- exp(c(yMat%*%beta_rho))
  h_val_s <- COMPUTE_H(beta_rho, theta, c_ps, sData, e_t_tau,
                       tau_x_external, piVal, coef_y_x_s, sigma_y_x_s,
                       e_s_rho2_psi_x_external, e_s_rho2_x_external)
  e_s_rho_x_s <- E_S_RHO_X(beta_rho, 1, sData[,-1], coef_y_x_s, sigma_y_x_s)
  s_val_s <- Compute_S(beta_rho, sData[,-1], coef_y_x_s, sigma_y_x_s, e_s_rho_x_s, ispar, parameters, c_ps)
  dVec_first <- COMPUTE_1ST_PART_D(theta, e_s_rho_y, e_s_rho_y2, e_s_rho_y3, c_ps, rho_val_s, h_val_s, s_val_s, piVal)
  
  return(dVec_first)
}

COMPUTE_EFFICIENT_IF_THETA <- function(theta, beta_rho, sData, tData,
                                       ispar, parameters, piVal,
                                       tDat_ext, coef_y_x_s, sigma_y_x_s,
                                       xList, wList) {
  num_of_source <- nrow(sData)
  num_of_target <- nrow(tData)
  
  # First Part
  yVec <- sData[,"Y"]
  yMat <- cbind(yVec, yVec^2)
  rho_val_s <- exp(c(yMat%*%beta_rho))
  c_ps <- E_S_RHO(beta_rho, ispar, parameters)
  first_part_0 <- 1/piVal*rho_val_s/c_ps*(yVec-theta)
  first_part <- c(first_part_0, rep(0,num_of_target))
  
  # Second Part
  second_part_0 <- c(1/piVal*rho_val_s/c_ps, -rep(1/(1-piVal),num_of_target))
  e_s_rho_x_ext <- E_S_RHO_X(beta_rho, 1, tDat_ext, coef_y_x_s, sigma_y_x_s)
  e_s_rho2_x_ext <- E_S_RHO_X(beta_rho, 2, tDat_ext, coef_y_x_s, sigma_y_x_s)
  tau_x_external <- COMPUTE_TAU(e_s_rho_x_ext, e_s_rho2_x_ext, c_ps, piVal)
  e_t_tau <- mean(tau_x_external)
  
  xMatAll <- rbind(sData[,-1], tData)
  e_s_rho2_x_all <- E_S_RHO_X(beta_rho, 2, xMatAll, coef_y_x_s, sigma_y_x_s)
  e_s_rho_x_all <- E_S_RHO_X(beta_rho, 1, xMatAll, coef_y_x_s, sigma_y_x_s)
  tau_x_all <- COMPUTE_TAU(e_s_rho_x, e_s_rho2_x, c_ps, piVal)
  
  e_s_rho2_psi_x_all <- E_S_RHO2_PSI_X(beta_rho, xMatAll, coef_y_x_s, sigma_y_x_s, xList, wList)
  e_s_rho2_psi_x_ext <- E_S_RHO2_PSI_X(beta_rho, tDat_ext, coef_y_x_s, sigma_y_x_s, xList, wList)
  
  aVec <- COMPUTE_A(beta_rho, e_t_tau, tau_x_all, tau_x_external,
                    e_s_rho2_psi_x_all, e_s_rho2_x_all, e_s_rho2_psi_x_ext, e_s_rho2_x_ext)
  bVec <- COMPUTE_B(tau_x_all, e_t_tau, c_ps)
  second_part <- second_part_0*(aVec-bVec)
  
  d1st <- COMPUTE_1ST_PART_D_WRAPPER(theta, beta_rho, ispar, parameters,
                                     sData, c_ps, tau_x_external, piVal,
                                     coef_y_x_s, sigma_y_x_s, e_s_rho2_psi_x_ext, e_s_rho2_x_ext)
  covMat <- EstimateBetaCovMat_CPP(beta_rho, sData, tData, piVal, tDat_ext, coef_y_x_s, sigma_y_x_s, ispar, parameters, xList, wList)
  dMat <- matrix(d1st, nrow = 1) %*% covMat
  sEff <- ComputeEfficientScore_CPP(beta_rho, sData, tData, piVal, tDat_ext, coef_y_x_s, sigma_y_x_s, ispar, parameters, xList, wList)
  third_part <- t(dMat %*% t(sEff))
  
  phi <- c(first_part)+c(second_part)+c(third_part)
  
  return(phi)
}
