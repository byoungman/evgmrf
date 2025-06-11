// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double xieps = 0.0001;

// [[Rcpp::export(.tppu0)]]
double tppu0(arma::vec pars, double u, double w)
{
  
  double nllh = 0.0;
  
  double mu = pars(0);
  double lpsi = pars(1);
  double txi = pars(2);
  double xi = 1.5 / (1.0 + exp(-txi)) - 1.0;
  
  double ee1;
  
  ee1 = xi * (u - mu) / exp(lpsi);
  nllh += w * R_pow(1.0 + ee1, -1/xi);
  
  return(nllh);
  
}

// [[Rcpp::export(.tppu1)]]
arma::vec tppu1(arma::vec pars, double u, double w)
{
  
  arma::vec g(3);
  
  double mu = pars(0);
  double lpsi = pars(1);
  double txi = pars(2);
  
  double ee2, ee3, ee5, ee6, ee7, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
  double ee20, ee21, ee24, ee25, ee26, ee29;
  double ee30, ee31, ee32, ee33;
  
  ee2 = exp(-txi);
  ee3 = 1 + ee2;
  ee5 = 1.5/ee3 - 1;
  ee6 = exp(lpsi);
  ee7 = u - mu;
  ee9 = ee5 * ee7/ee6;
  ee10 = 1/ee5;
  ee11 = ee9 + 1;
  ee12 = 1 + ee10;
  ee13 = R_pow(ee11, ee12);
  ee14 = R_pow(ee3, 2);
  ee15 = log1p(ee9);
  ee16 = R_pow(ee11, ee10 + 2);
  ee17 = ee16 * ee6;
  ee18 = ee13 * ee6;
  ee20 = R_pow(ee11, ee10);
  ee21 = ee12 * ee5;
  ee24 = ee14 * ee5;
  ee25 = (1.5 * (ee15/(ee13 * R_pow(ee5, 2))) - 1.5 * (ee12 * ee7/ee17)) * ee2;
  ee26 = ee7/ee18;
  ee29 = ee21 * ee7/ee17;
  ee30 = ee3 * ee5;
  ee31 = ee14 * ee6;
  ee32 = (1.5 * (ee15/(ee20 * ee5)) - 1.5 * ee26) * ee2;
  ee33 = 1/ee13;
  
  g(0) = 1/ee18;
  g(1) = ee26;
  g(2) = ee32/ee24;
  // out(j, 3) = ee21/(ee16 * R_pow(ee6, 2));
  // out(j, 4) = (ee29 - ee33)/ee6;
  // out(j, 5) = ee25/ee31;
  // out(j, 6) = -((ee33 - ee29) * ee7/ee6);
  // out(j, 7) = ee25 * ee7/ee31;
  // out(j, 8) = ((((2.25/ee30 - 3) * ee2/ee3 + 1.5)/ee13 - 1.5 * (ee25/ee14)) * ee7/ee6 +
  //   ((2.25 * (ee2 * ee7/(ee11 * ee14 * ee6)) -
  //   ((4.5/ee30 - 3) * ee2/ee3 + 1.5) * ee15)/ee20 + 1.5 * (ee32 * ee15/ee24))/ee5) * ee2/ee24;

  return w * g;
  
}

// [[Rcpp::export(.tppu2)]]
arma::mat tppu2(arma::vec pars, double u, double w)
{
  
  arma::mat h(3, 3);
  
  double mu = pars(0);
  double lpsi = pars(1);
  double txi = pars(2);
  
  double ee2, ee3, ee5, ee6, ee7, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
  double ee20, ee21, ee24, ee25, ee26, ee29;
  double ee30, ee31, ee32, ee33;
  
  ee2 = exp(-txi);
  ee3 = 1 + ee2;
  ee5 = 1.5/ee3 - 1;
  ee6 = exp(lpsi);
  ee7 = u - mu;
  ee9 = ee5 * ee7/ee6;
  ee10 = 1/ee5;
  ee11 = ee9 + 1;
  ee12 = 1 + ee10;
  ee13 = R_pow(ee11, ee12);
  ee14 = R_pow(ee3, 2);
  ee15 = log1p(ee9);
  ee16 = R_pow(ee11, ee10 + 2);
  ee17 = ee16 * ee6;
  ee18 = ee13 * ee6;
  ee20 = R_pow(ee11, ee10);
  ee21 = ee12 * ee5;
  ee24 = ee14 * ee5;
  ee25 = (1.5 * (ee15/(ee13 * R_pow(ee5, 2))) - 1.5 * (ee12 * ee7/ee17)) * ee2;
  ee26 = ee7/ee18;
  ee29 = ee21 * ee7/ee17;
  ee30 = ee3 * ee5;
  ee31 = ee14 * ee6;
  ee32 = (1.5 * (ee15/(ee20 * ee5)) - 1.5 * ee26) * ee2;
  ee33 = 1/ee13;
  
  // g(0) = 1/ee18;
  // g(1) = ee26;
  // g(2) = ee32/ee24;
  
  h(0, 0) = ee21/(ee16 * R_pow(ee6, 2));
  h(1, 0) = (ee29 - ee33)/ee6;
  h(2, 0) = ee25/ee31;
  h(1, 1) = -((ee33 - ee29) * ee7/ee6);
  h(2, 1) = ee25 * ee7/ee31;
  h(2, 2) = ((((2.25/ee30 - 3) * ee2/ee3 + 1.5)/ee13 - 1.5 * (ee25/ee14)) * ee7/ee6 +
    ((2.25 * (ee2 * ee7/(ee11 * ee14 * ee6)) -
    ((4.5/ee30 - 3) * ee2/ee3 + 1.5) * ee15)/ee20 + 1.5 * (ee32 * ee15/ee24))/ee5) * ee2/ee24;

  h(0, 1) = h(1, 0);
  h(0, 2) = h(2, 0);
  h(1, 2) = h(2, 1);
  
  return w * h;
  
}

// [[Rcpp::export(.tppz0)]]
double tppz0(arma::vec pars, arma::vec zv)
{
  
  int n = zv.size(); // number of locations

  double nllh = 0.0;
  
  double mu = pars(0);
  double lpsi = pars(1);
  double txi = pars(2);
  double xi = 1.5 / (1.0 + exp(-txi)) - 1.0;
  
  double w = 1.0;
  double z;
  
  double ee1, ee2;
  
  for (int i=0; i < n; i++) {
    
    z = zv(i);
    
    ee1 = xi * (z - mu) / exp(lpsi);
    
    if (ee1 <= -1.0) {
      
      nllh = 1e20;
      break;
      
    } else {
      
      nllh += w * (lpsi + (1 / xi + 1) * log1p(ee1));
      
    }
    
  }
  
  return(nllh);
  
}

// [[Rcpp::export(.tppz1)]]
arma::vec tppz1(arma::vec pars, arma::vec zv)
{
  
  int n = zv.size(); // number of locations
  
  arma::vec g(3, arma::fill::zeros);
  
  double mu = pars(0);
  double lpsi = pars(1);
  double txi = pars(2);
  double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee17, ee18, ee19;
  double ee20;
  double z;
  double w = 1.0;
  
  for (int i=0; i < n; i++) {
    
  z = zv(i);
  
  ee2 = exp(-txi);
  ee3 = 1 + ee2;
  ee5 = 1.5/ee3 - 1;
  ee6 = exp(lpsi);
  ee7 = z - mu;
  ee8 = ee5 * ee7;
  ee9 = ee8/ee6;
  ee10 = ee9 + 1;
  ee11 = ee10 * ee6;
  ee12 = 1 + 1/ee5;
  ee13 = R_pow(ee3, 2);
  ee14 = ee8/ee11;
  ee15 = R_pow(ee5, 2);
  ee17 = ee10 * ee13 * ee6;
  ee18 = (ee12 * (1.5 - 1.5 * ee14) - 1.5/ee5) * ee2;
  ee19 = ee12 * ee5;
  ee20 = log1p(ee9);
  
  g(0) += -(ee19/ee11);
  g(1) += 1 - ee19 * ee7/ee11;
  g(2) += (1.5 * (ee12 * ee7/ee11) - 1.5 * (ee20/ee15)) * ee2/ee13;

  // h(0, 0) =  - (ee12 * ee15/(R_pow(ee10, 2) * R_pow(ee6, 2)));
  // h(1, 0) = (1 - ee14) * ee12 * ee5/ee11;
  // h(2, 0) = -(ee18/ee17);
  // h(1, 1) = -((ee14 - 1) * ee12 * ee5 * ee7/ee11);
  // h(2, 1) = -(ee18 * ee7/ee17);
  // h(2, 2) = -(((((2.25 * (ee7/(ee10 * ee3 * ee6)) - 3) * ee2/ee3 +
  //   1.5) * ee12 + 2.25 * (ee2/(ee13 * ee15))) * ee7/ee11 +
  //   (2.25 * (ee2 * ee7/ee17) - ((4.5/(ee3 * ee5) - 3) * ee2/ee3 +
  //   1.5) * ee20)/ee15) * ee2/ee13);
  
  // h(0, 1) = h(1, 0);
  // h(0, 2) = h(2, 0);
  // h(2, 1) = h(1, 2);
  
  }
  
  return g;
  
}

// [[Rcpp::export(.tppz2)]]
arma::mat tppz2(arma::vec pars, arma::vec zv)
{
  
  int n = zv.size(); // number of locations
  
  arma::mat h(3, 3, arma::fill::zeros);
  
  double mu = pars(0);
  double lpsi = pars(1);
  double txi = pars(2);
  
  double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee17, ee18, ee19;
  double ee20;
  double z;
  double w = 1.0;
  
  for (int i=0; i < n; i++) {
    
    z = zv(i);
    
    ee2 = exp(-txi);
    ee3 = 1 + ee2;
    ee5 = 1.5/ee3 - 1;
    ee6 = exp(lpsi);
    ee7 = z - mu;
    ee8 = ee5 * ee7;
    ee9 = ee8/ee6;
    ee10 = ee9 + 1;
    ee11 = ee10 * ee6;
    ee12 = 1 + 1/ee5;
    ee13 = R_pow(ee3, 2);
    ee14 = ee8/ee11;
    ee15 = R_pow(ee5, 2);
    ee17 = ee10 * ee13 * ee6;
    ee18 = (ee12 * (1.5 - 1.5 * ee14) - 1.5/ee5) * ee2;
    ee19 = ee12 * ee5;
    ee20 = log1p(ee9);
    
    // g(0) = -(ee19/ee11);
    // g(1) = 1 - ee19 * ee7/ee11;
    // g(2) = (1.5 * (ee12 * ee7/ee11) - 1.5 * (ee20/ee15)) * ee2/ee13;
    
    h(0, 0) +=  - (ee12 * ee15/(R_pow(ee10, 2) * R_pow(ee6, 2)));
    h(1, 0) += (1 - ee14) * ee12 * ee5/ee11;
    h(2, 0) += -(ee18/ee17);
    h(1, 1) += -((ee14 - 1) * ee12 * ee5 * ee7/ee11);
    h(2, 1) += -(ee18 * ee7/ee17);
    h(2, 2) += -(((((2.25 * (ee7/(ee10 * ee3 * ee6)) - 3) * ee2/ee3 +
      1.5) * ee12 + 2.25 * (ee2/(ee13 * ee15))) * ee7/ee11 +
      (2.25 * (ee2 * ee7/ee17) - ((4.5/(ee3 * ee5) - 3) * ee2/ee3 +
      1.5) * ee20)/ee15) * ee2/ee13);
    
  }
    
  h(0, 1) = h(1, 0);
  h(0, 2) = h(2, 0);
  h(1, 2) = h(2, 1);
  
  return h;
  
}

  // [[Rcpp::export(.tpp0)]]
  double tpp0(arma::vec pars, arma::vec zv, double w, double u, double delta)
  {
    
    double out = tppu0(pars, u, w);
    out += tppz0(pars, zv);
    out += 0.5 * delta * (pars(2) - 0.7) * (pars(2) - 0.7);
    
    if(!arma::is_finite(out)) {
      out = 1e20;
    }
    
    return out;
    
  }
  
  
  
  // [[Rcpp::export(.tpp1)]]
  arma::vec tpp1(arma::vec pars, arma::vec zv, double w, double u, double delta)
  {
    
    arma::vec out = tppu1(pars, u, w);
    out += tppz1(pars, zv);
    out(2) += delta * (pars(2) - 0.7);
    
    return out;
    
  }
  
  // [[Rcpp::export(.tpp2)]]
  arma::mat tpp2(arma::vec pars, arma::vec zv, double w, double u, double delta)
  {
    
    arma::mat out = tppu2(pars, u, w);
    out += tppz2(pars, zv);
    out(2, 2) += delta;
    
    return out;
    
  }
  