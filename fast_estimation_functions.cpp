#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector E_S_RHO_X_CPP(NumericVector beta_rho, double rho_pwr,
                            NumericMatrix x_mat_no_intercept, NumericVector coef_y_x_s, double sigma_y_x_s,
                            NumericVector xList, NumericVector wList) {
  int i, j;
  int num_row = x_mat_no_intercept.nrow(), num_col = x_mat_no_intercept.ncol(), ghpoints = xList.length();
  NumericVector mean_y_x_s(num_row), e_s_rho_x(num_row);
  
  for (i=0; i<num_row; i++) {
    mean_y_x_s(i) = coef_y_x_s(0);
    for (j=0; j<num_col; j++) {
      mean_y_x_s(i) += coef_y_x_s(j+1)*x_mat_no_intercept(i, j);
    }//
  }
  
  double yTemp;
  for (i=0; i<num_row; i++) {
    yTemp = 0;
    for (j=0; j<ghpoints; j++) {
      yTemp = sqrt(2)*sigma_y_x_s*xList(j)+mean_y_x_s(i);
      e_s_rho_x(i) += exp(rho_pwr*yTemp*beta_rho(0)+rho_pwr*yTemp*yTemp*beta_rho(1))*wList(j)/sqrt(3.1415926);
    }
  }
  
  return e_s_rho_x;
}

// [[Rcpp::export]]
double E_S_RHO_CPP(NumericVector beta_rho, bool ispar, List parameters) {
  double e_s_rho=0;
  
  if (ispar) {
    double Mu_Y_S = parameters["mu"], Sig_Y_S = parameters["sigma"];
    NumericVector xList = parameters["xList"], wList = parameters["wList"];
    int ghpoint = xList.length();
    
    double yTemp;
    for (int i=0; i<ghpoint; i++) {
      yTemp = sqrt(2)*Sig_Y_S*xList(i)+Mu_Y_S;
      e_s_rho += exp(yTemp*beta_rho(0)+yTemp*yTemp*beta_rho(1))*wList(i);
    }
    e_s_rho /= sqrt(3.1415926);
  }
  else {
    NumericVector y_vec = parameters[0];
    int yLen = y_vec.length();
    
    for (int i=0; i<yLen; i++) {
      e_s_rho += exp(beta_rho(0)*y_vec(i)+beta_rho(1)*y_vec(i)*y_vec(i));
    }
    e_s_rho /= yLen;
  }
  
  return e_s_rho;
}

// [[Rcpp::export]]
NumericVector COMPUTE_TAU_CPP(NumericVector e_s_rho_x, NumericVector e_s_rho2_x, double c_ps, double piVal) {
  int num_of_x = e_s_rho_x.length();
  double tmp;
  NumericVector tauX(num_of_x);
  for (int i=0; i<num_of_x; i++) {
    tmp = e_s_rho2_x(i)/e_s_rho_x(i)/c_ps/piVal;
    tauX(i) = tmp/(tmp+1/(1-piVal));
  }
  
  return tauX;
}

// [[Rcpp::export]]
double E_T_TAU_CPP(NumericVector e_s_rho_x, NumericVector e_s_rho2_x, double c_ps, double piVal) {
  NumericVector tau_vec = COMPUTE_TAU_CPP(e_s_rho_x, e_s_rho2_x, c_ps, piVal);
  double mean_tau = mean(tau_vec);
  
  return mean_tau;
}

// [[Rcpp::export]]
NumericVector E_T_D_LOG_RHO_DIV_D_BETA_CPP(NumericVector beta_rho, bool ispar, List parameters, double c_ps) {
  
  NumericVector outBeta(2);
  
  if (ispar) {
    double Mu_Y_S = parameters["mu"], Sig_Y_S = parameters["sigma"];
    NumericVector xList = parameters["xList"], wList = parameters["wList"];
    int ghnum = xList.length();
    
    double yTemp, rhoTmp;
    for (int i=0; i<ghnum; i++) {
      yTemp = sqrt(2)*Sig_Y_S*xList(i)+Mu_Y_S;
      rhoTmp = exp(yTemp*beta_rho(0)+yTemp*yTemp*beta_rho(1));
      outBeta(0) += wList(i)*rhoTmp*yTemp;
      outBeta(1) += wList(i)*rhoTmp*yTemp*yTemp;
    }
    outBeta(0) /= (c_ps*sqrt(3.1415926));
    outBeta(1) /= (c_ps*sqrt(3.1415926));
  }
  else
  {
    NumericVector y_s_external = parameters[0];
    int y_num = y_s_external.length();
    double yCurrent, rhoCurrent;
    for (int i=0; i<y_num; i++) {
      
      yCurrent = y_s_external(i);
      rhoCurrent = exp(yCurrent*beta_rho(0)+yCurrent*yCurrent*beta_rho(1));
      
      outBeta(0) += rhoCurrent*yCurrent;
      outBeta(1) += rhoCurrent*yCurrent*yCurrent;
      
    }
    outBeta(0) /= (y_num*c_ps);
    outBeta(1) /= (y_num*c_ps);
  }
  
  return outBeta;
}

// [[Rcpp::export]]
NumericMatrix E_T_D_LOG_RHO_DIV_D_BETA_X_CPP(NumericVector beta_rho, NumericMatrix x_mat_no_intercept,
                                             NumericVector coef_y_x_s, double sigma_y_x_s, NumericVector e_s_rho_x,
                                             NumericVector xList, NumericVector wList) {
  int num_of_x = x_mat_no_intercept.nrow(), num_of_covariate = x_mat_no_intercept.ncol(), ghnum = xList.length();
  int i, j;
  NumericMatrix out_mat(num_of_x, 2);
  NumericVector y_x_mean(num_of_x);
  
  for (i=0; i<num_of_x; i++) {
    y_x_mean(i) = coef_y_x_s(0);
    for (j=0; j<num_of_covariate; j++) {
      y_x_mean(i) += coef_y_x_s(j+1)*x_mat_no_intercept(i,j);
    }
  }
  
  double yTmp, rhoSrc, tmp1, tmp2;
  for (i=0; i<num_of_x; i++) {
    tmp1 = 0;
    tmp2 = 0;
    for (j=0; j<ghnum; j++) {
      yTmp = sqrt(2)*sigma_y_x_s*xList(j)+y_x_mean(i);
      rhoSrc = exp(yTmp*beta_rho(0)+yTmp*yTmp*beta_rho(1));
      tmp1 += wList(j)*rhoSrc*yTmp;
      tmp2 += wList(j)*rhoSrc*yTmp*yTmp;
    }
    tmp1 /= (sqrt(3.1415926)*e_s_rho_x(i));
    tmp2 /= (sqrt(3.1415926)*e_s_rho_x(i));
    
    out_mat(i,0) = tmp1;
    out_mat(i,1) = tmp2;
  }
  
  return out_mat;
}

// [[Rcpp::export]]
NumericMatrix Compute_S_CPP(NumericVector beta_rho, NumericMatrix x_mat_no_intercept,
                            NumericVector coef_y_x_s, double sigma_y_x_s, NumericVector e_s_rho_x,
                            NumericVector xList, NumericVector wList,
                            bool ispar, List parameters, double c_ps) {
  NumericMatrix e_t_d_log_rho_div_d_beta_x = E_T_D_LOG_RHO_DIV_D_BETA_X_CPP(beta_rho, x_mat_no_intercept, coef_y_x_s, sigma_y_x_s,
                                                                            e_s_rho_x, xList, wList);
  NumericVector e_t_d_log_rho_div_d_beta = E_T_D_LOG_RHO_DIV_D_BETA_CPP(beta_rho, ispar, parameters, c_ps);
  int num_of_x = e_t_d_log_rho_div_d_beta_x.nrow();
  for (int i=0; i< num_of_x; i++) {
    for (int j=0; j<2; j++) {
      e_t_d_log_rho_div_d_beta_x(i,j) -= e_t_d_log_rho_div_d_beta(j);
    }
  }
  
  return e_t_d_log_rho_div_d_beta_x;
}

// [[Rcpp::export]]
NumericMatrix ComputeEfficientScore_CPP(NumericVector beta_rho, NumericMatrix sData, NumericMatrix tData,
                                        double piVal, NumericMatrix tDat_ext, NumericVector coef_y_x_s, double sigma_y_x_s,
                                        bool ispar, List parameters, NumericVector xList, NumericVector wList) {
  int num_of_source = sData.nrow(), num_of_target = tData.nrow(), num_of_cov = tData.ncol();
  int num_of_total = num_of_source+num_of_target;
  
  NumericVector y_internal = sData( _ , 0);
  NumericVector rho_internal(num_of_source);
  for (int i=0; i<num_of_source; i++) {
    rho_internal(i) = exp(y_internal(i)*beta_rho(0)+y_internal(i)*y_internal(i)*beta_rho(1));
  }
  double c_ps = E_S_RHO_CPP(beta_rho, ispar, parameters);
  
  NumericMatrix xMatAll(num_of_total, num_of_cov);
  for (int i=0; i<num_of_total; i++) {
    for (int j=0; j<num_of_cov; j++) {
      if (i < num_of_source)
      {
        xMatAll(i,j) = sData(i, j+1);
      }
      else
      {
        xMatAll(i,j) = tData(i-num_of_source, j);
      }
    }
  }
  
  NumericVector multiplier(num_of_total);
  for(int i=0; i<num_of_total; i++) {
    if(i<num_of_source)
    {
      multiplier(i) = rho_internal(i)/c_ps/piVal;
    }
    else
    {
      multiplier(i) = -1/(1-piVal);
    }
  }
  
  NumericVector e_s_rho2_x = E_S_RHO_X_CPP(beta_rho, 2, xMatAll, coef_y_x_s, sigma_y_x_s, xList, wList);
  NumericVector e_s_rho_x = E_S_RHO_X_CPP(beta_rho, 1, xMatAll, coef_y_x_s, sigma_y_x_s, xList, wList);
  NumericVector tau_x_internal = COMPUTE_TAU_CPP(e_s_rho_x, e_s_rho2_x, c_ps, piVal);
  NumericMatrix s_x_internal = Compute_S_CPP(beta_rho, xMatAll, coef_y_x_s, sigma_y_x_s, e_s_rho_x, xList, wList,
                                             ispar, parameters, c_ps);

  NumericVector e_s_rho_x_ext = E_S_RHO_X_CPP(beta_rho, 1, tDat_ext, coef_y_x_s, sigma_y_x_s, xList, wList);
  NumericVector e_s_rho2_x_ext = E_S_RHO_X_CPP(beta_rho, 2, tDat_ext, coef_y_x_s, sigma_y_x_s, xList, wList);
  NumericVector tau_x_external = COMPUTE_TAU_CPP(e_s_rho_x_ext, e_s_rho2_x_ext, c_ps, piVal);
  NumericMatrix s_x_external = Compute_S_CPP(beta_rho, tDat_ext, coef_y_x_s, sigma_y_x_s, e_s_rho_x_ext, xList, wList,
                                             ispar, parameters, c_ps);

  double e_t_tau = mean(tau_x_external);
  int num_of_ext = tDat_ext.nrow();
  NumericVector e_t_tau_s(2);
  for(int i=0; i<num_of_ext; i++) {
    e_t_tau_s(0) += tau_x_external(i)*s_x_external(i,0);
    e_t_tau_s(1) += tau_x_external(i)*s_x_external(i,1);
  }
  e_t_tau_s(0) /= num_of_ext;
  e_t_tau_s(1) /= num_of_ext;

  double b2ndpart1, b2ndpart2;
  NumericMatrix b1_x(num_of_total, 2);
  for(int i=0; i<num_of_total; i++) {
    b2ndpart1 = s_x_internal(i,0)-e_t_tau_s(0)/(e_t_tau-1);
    b2ndpart2 = s_x_internal(i,1)-e_t_tau_s(1)/(e_t_tau-1);

    b1_x(i,0) = -multiplier(i)*(1-piVal)*(1-tau_x_internal(i))*b2ndpart1;
    b1_x(i,1) = -multiplier(i)*(1-piVal)*(1-tau_x_internal(i))*b2ndpart2;
  }
  
  return b1_x;
}

// [[Rcpp::export]]
double EstimateBetaFunc_CPP(NumericVector beta_rho, NumericMatrix sData, NumericMatrix tData,
                            double piVal, NumericMatrix tDat_ext, NumericVector coef_y_x_s, double sigma_y_x_s,
                            bool ispar, List parameters, NumericVector xList, NumericVector wList, bool weights = false) {
  double out;
  NumericMatrix Seff = ComputeEfficientScore_CPP(beta_rho, sData, tData, piVal, tDat_ext, coef_y_x_s, sigma_y_x_s,
                                                 ispar, parameters, xList, wList);
  int num_of_row = Seff.nrow();
  
  NumericVector rexpVec(num_of_row);
  if (!weights) {
    for (int i=0; i<num_of_row; i++) {
      rexpVec(i) = 1;
    }
  }
  else {
    for (int i=0; i<num_of_row; i++) {
      rexpVec(i) = R::rexp(1);
    }
  }
  
  double sumTmp1 = 0, sumTmp2 = 0;
  for(int i = 0; i<num_of_row; i++) {
    sumTmp1 += rexpVec(i)*Seff(i,0);
    sumTmp2 += rexpVec(i)*Seff(i,1);
  }
  sumTmp1 /= num_of_row;
  sumTmp2 /= num_of_row;
  
  out = sumTmp1*sumTmp1+sumTmp2*sumTmp2;

  return out;
}

// [[Rcpp::export]]
double EstimateBetaFuncSum_CPP(NumericVector beta_rho, NumericMatrix sData, NumericMatrix tData,
                               double piVal, NumericMatrix tDat_ext, NumericVector coef_y_x_s, double sigma_y_x_s,
                               bool ispar, List parameters, NumericVector xList, NumericVector wList, bool weights = false) {
  double out;
  NumericMatrix Seff = ComputeEfficientScore_CPP(beta_rho, sData, tData, piVal, tDat_ext, coef_y_x_s, sigma_y_x_s,
                                                 ispar, parameters, xList, wList);
  int num_of_row = Seff.nrow();
  
  NumericVector rexpVec(num_of_row);
  if (!weights) {
    for (int i=0; i<num_of_row; i++) {
      rexpVec(i) = 1;
    }
  }
  else {
    for (int i=0; i<num_of_row; i++) {
      rexpVec(i) = R::rexp(1);
    }
  }
  
  double sumTmp1 = 0, sumTmp2 = 0;
  for(int i = 0; i<num_of_row; i++) {
    sumTmp1 += rexpVec(i)*Seff(i,0);
    sumTmp2 += rexpVec(i)*Seff(i,1);
  }
  
  out = sumTmp1*sumTmp1+sumTmp2*sumTmp2;
  
  return out;
}

// [[Rcpp::export]]
NumericVector EstimateBetaVarFunc_CPP(NumericVector beta_rho, NumericMatrix sData, NumericMatrix tData,
                                  double piVal, NumericMatrix tDat_ext, NumericVector coef_y_x_s, double sigma_y_x_s,
                                  bool ispar, List parameters, NumericVector xList, NumericVector wList) {
  NumericMatrix Seff = ComputeEfficientScore_CPP(beta_rho, sData, tData, piVal, tDat_ext, coef_y_x_s, sigma_y_x_s,
                                                 ispar, parameters, xList, wList);
  NumericMatrix tmpMat(2,2);
  NumericVector outVec(2);
  
  int num_of_row = Seff.nrow();
  for(int i=0; i<num_of_row; i++) {
    tmpMat(0,0) += Seff(i,0)*Seff(i,0);
    tmpMat(0,1) += Seff(i,0)*Seff(i,1);
    tmpMat(1,0) += Seff(i,0)*Seff(i,1);
    tmpMat(1,1) += Seff(i,1)*Seff(i,1);
  }
  
  for(int i=0; i<2; i++){
    for(int j=0; j<2; j++){
      tmpMat(i,j) /= num_of_row;
    }
  }
  
  double det = tmpMat(0,0)*tmpMat(1,1)-tmpMat(1,0)*tmpMat(0,1);
  outVec(0) = sqrt(tmpMat(1,1)/det/num_of_row);
  outVec(1) = sqrt(tmpMat(0,0)/det/num_of_row);

  return outVec;
}

