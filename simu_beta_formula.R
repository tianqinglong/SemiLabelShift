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
# n <- 2000
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
beta_rho <- trueBetaRho
B1 <- 2000 # Monte-Carlo Sample Size
# B2 <- 10 # Bootstrap (Perturbation Size)

gh_num <- 10
ghxw <- gaussHermiteData(gh_num)
xList <- ghxw$x
wList <- ghxw$w
################################################
set.seed(888)
for (n in c(200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000)){
  for (mnratio in c(0.2, 0.6, 1, 1.4, 1.8, 2)) {
    m <- mnratio*n
    data_list_mc <- mclapply(1:B1, function(x) {
      Generate_Dat(n, m, Mu_YX, SigMat_YX, Mu_Y_T, Sig_Y_T)
    }, mc.cores = detectCores())
    
    ################################################
    results_output <- mclapply(data_list_mc, function(dataList) {
      sData <- dataList$sDat
      tData <- dataList$tDat
      piVal <- n/(n+m)
      ispar <- T
      
      fityx <- lm(Y~., data = as.data.frame(sData))
      coef_y_x_s_hat <- coef(fityx)
      sigma_y_x_s_hat <- sigma(fityx)

      Mu_Y_S_hat <- mean(sData[,"Y"])
      Sig_Y_S_hat <- sd(sData[,"Y"])
      
      parameters <- list(y_vec = sData[,"Y"], mu = Mu_Y_S_hat, sigma = Sig_Y_S_hat, xList = xList, wList= wList)
      
      betaHat <- optim(beta_rho, EstimateBetaFunc_CPP, sData = sData, tData = tData, piVal = piVal,
                       tDat_ext = tData, coef_y_x_s = coef_y_x_s_hat, sigma_y_x_s = sigma_y_x_s_hat, ispar = ispar,
                       parameters = parameters, xList = xList, wList = wList)
      betaHat <- betaHat$par
      # betaSd1 <- EstimateBetaVarCenterFunc(betaHat, sData, tData, piVal, tData, coef_y_x_s_true, sigma_y_x_s_true, ispar, parameters, xList = xList, wList = wList)
      betaSd1 <- EstimateBetaVarFunc_CPP(betaHat, sData, tData, piVal, tData, coef_y_x_s_hat, sigma_y_x_s_hat, ispar, parameters, xList, wList)
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
    saveRDS(results_output, file = paste("dat1/output_ratio", mnratio, "_n", n,"_.RDS", sep = ""))
  }
}
################################################
read_in_data <- function(filename, dir) {
  rt <- str_match(filename, "ratio(.*?)_n")[2] %>% as.numeric()
  n <- str_match(filename, "_n(.*?)_.RDS")[2] %>% as.numeric()
  dat <- readRDS(paste(dir, filename, sep = ""))
  
  return(list(rt = rt, dat = dat, n = n))
}

dat_dir <- "dat1/"
filenames <- list.files(dat_dir)
out <- NULL

for (filename in filenames) {
  rtDat <- read_in_data(filename, dat_dir)
  rt <- rtDat[[1]]
  dat <- rtDat[[2]]
  n <- rtDat[[3]]
  
  sapply(dat, function(out) {
    out$CP1
  }) %>% rowMeans() -> Cp
  
  sapply(dat, function(out) {
    out$Sd
  }) %>% rowMeans() -> Sd
  
  sapply(dat, function(out) {
    out$BetaHat
  }) %>% apply(MARGIN = 1, sd) -> Se
  
  sapply(dat, function(out) {
    out$BetaHat-out$TrueBeta
  }) %>% rowMeans() -> Bias
  
  out <- rbind(out, c(rt, Cp, Sd, Se, Bias, n))
}
colnames(out) <- c("ratio", "Cp_Beta_0", "Cp_Beta_1", "Sd_Beta_0", "Sd_Beta_1", "Se_Beta_0", "Se_Beta_1", "Bias_0", "Bias_1", "n")
out <- as.data.frame(out)

# Coverage Probability
outCP <- out[,c(1, 2, 3, 10)]
outCP <- reshape2::melt(outCP, value.name = "CP", measure.vars = c("Cp_Beta_0", "Cp_Beta_1"), variable.name = "Parameter")
ggplot(aes(x = n, y = CP), data = outCP)+geom_line(aes(col = Parameter, linetype = Parameter))+geom_point(aes(col = Parameter, shape = Parameter))+
  geom_hline(aes(yintercept = 0.95),linetype="dashed")+
  scale_color_discrete(labels = c(expression(paste(beta[0])), expression(paste(beta[1]))))+
  scale_shape_discrete(labels = c(expression(paste(beta[0])), expression(paste(beta[1]))))+
  scale_linetype_discrete(labels = c(expression(paste(beta[0])), expression(paste(beta[1]))))+
  facet_wrap( ~ ratio, nrow=2, labeller = labeller(ratio = c("0.25" = "m/n=0.25", "0.5" = "m/n=0.5", "1" = "m/n=1",
                                                             "1.5" = "m/n=1.5", "2" = "m/n=2")))+ylim(c(0.9, 1))+
  ylab("Coverage Probability")

# SE and SD
outSeSd <- out[,c(1, 4, 6, 5, 7, 10)]
outSeSd <- reshape2::melt(outSeSd, value.name = "SE_or_SD", measure.vars = c("Sd_Beta_0", "Se_Beta_0", "Sd_Beta_1", "Se_Beta_1"))
outSeSd$Beta <- ifelse(outSeSd$variable %in% c("Sd_Beta_0", "Se_Beta_0"), "beta_0", "beta_1")
outSeSd$IsSe <- ifelse(outSeSd$variable %in% c("Sd_Beta_0", "Sd_Beta_1"), "Estimated", "Empirical")
outSeSd$ratio <- factor(outSeSd$ratio, labels = c("0.25" = 1, "0.5" = 2, "1" = 3,"1.5" = 4, "2" = 5))
outSeSd$FacetBeta <- ifelse(outSeSd$Beta == "beta_0", "beta[0]", "beta[1]")
outSeSd$Facetratio <- c('m/n*"=0.25"', 'm/n*"=0.5"', 'm/n*"=1"', 'm/n*"=1.5"', 'm/n*"=2"')[outSeSd$ratio]

ggplot(aes(x = n, y = SE_or_SD), data = outSeSd)+
  geom_line(aes(col = IsSe, linetype = IsSe))+
  geom_point(aes(col = IsSe, shape = IsSe))+
  scale_color_discrete(name = "Standard Deviation")+
  scale_shape_discrete(name = "Standard Deviation")+
  scale_linetype_discrete(name = "Standard Deviation")+
  facet_grid(FacetBeta ~ Facetratio, scales = "free_y", labeller = label_parsed)+
  ylab("Empirical v.s. Estimated Standard Deviation")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Bias
outBias <- out[,c(1, 8, 9, 10)]
outBias <- reshape2::melt(outBias, value.name = "Bias", measure.vars = c("Bias_0", "Bias_1"))
outBias$ratio <- factor(outBias$ratio, labels = c("0.25" = 1, "0.5" = 2, "1" = 3,"1.5" = 4, "2" = 5))
outBias$Facetratio <- c('m/n*"=0.25"', 'm/n*"=0.5"', 'm/n*"=1"', 'm/n*"=1.5"', 'm/n*"=2"')[outBias$ratio]
outBias$FacetBias <- ifelse(outBias$variable == "Bias_0", "beta[0]", "beta[1]")
ggplot(aes(x = n, y = Bias), data = outBias)+geom_line(aes(col = variable))+
  geom_point(aes(col = variable))+
  facet_grid(FacetBias~Facetratio, scales = "free_y", labeller = label_parsed)+
  scale_color_discrete(name = "Parameter",
                       labels = c(expression(paste(beta[0])), expression(paste(beta[1]))))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
################################################
