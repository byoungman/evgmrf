// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export(.rlarged0)]]
double rlarged0(arma::vec pars, arma::mat yv)
{
  
  int n = yv.n_rows; // number of locations
  int r = yv.n_cols; // number of locations
    
  double nllh = 0.0;
    
  double mu = pars(0);
  double lpsi = pars(1);
  double txi = pars(2);
  double xi = 1.5 / (1.0 + exp(-txi)) - 1.0;
  arma::rowvec y(r), ee1(r);
    
  double ee2;
    
  for (int i=0; i < n; i++) {
      
    y = yv.row(i);
    ee1 = xi * (y - mu) / exp(lpsi);
    
    for (int l=0; l < r; l++) {
      
      ee2 = 1.0 / xi;
      
      if (ee1(l) <= -1.0) {
        
        nllh = 1e20;
        break;
        
      } else {
        
        nllh += lpsi + (ee2 + 1.0) * log1p(ee1(l));
        
      }

    }
    
    nllh += R_pow(1 + ee1[r - 1], -ee2);

  }
  
return(nllh);

}

// [[Rcpp::export(.rlarged1)]]
arma::vec rlarged1(arma::vec pars, arma::mat yv)
{

  int n = yv.n_rows; // number of locations
  int r = yv.n_cols; // number of locations
  
  double mu = pars(0);
  double lpsi = pars(1);
  double txi = pars(2);
  arma::rowvec y(r), ee1(r);
  
  double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21, ee24, ee25, ee26, ee29;
  double ee30, ee31, ee32, ee33;
  
  arma::vec g(3, arma::fill::zeros);

    for (int i=0; i < n; i++) {
    
    y = yv.row(i);
    
    for (int l=0; l < r; l++) {

      ee2 = exp(-txi);
      ee3 = 1 + ee2;
      ee5 = 1.5/ee3 - 1;
      ee6 = exp(lpsi);
      ee7 = y(l) - mu;
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
      g(2) += (1.5 * (ee12 * ee7/ee11) - 1.5 * (ee20/ee15)) *
        ee2/ee13;

    }
    
      ee2 = exp(-txi);
      ee3 = 1 + ee2;
      ee5 = 1.5/ee3 - 1;
      ee6 = exp(lpsi);
      ee7 = y(r - 1) - mu;
      ee9 = ee5 * ee7/ee6;
      ee10 = 1/ee5;
      ee11 = ee9 + 1;
      ee12 = 1 + ee10;
      ee13 = R_pow(ee11, ee12);
      ee14 = R_pow(ee3, 2);
      ee15 = log1p(ee9);
      ee16 = R_pow(ee11, (ee10 + 2));
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

      g(0) += 1/ee18;
      g(1) += ee26;
      g(2) += ee32/ee24;
      
  }
  
return g;

}

// [[Rcpp::export(.rlarged2)]]
arma::mat rlarged2(arma::vec pars, arma::mat yv)
{
  
  int n = yv.n_rows; // number of locations
  int r = yv.n_cols; // number of locations
  
  double mu = pars(0);
  double lpsi = pars(1);
  double txi = pars(2);
  arma::rowvec y(r), ee1(r);
  
  double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21, ee24, ee25, ee26, ee29;
  double ee30, ee31, ee32, ee33;
  
  arma::mat H(3, 3, arma::fill::zeros);
  
  for (int j=0; j < n; j++) {
    
      y = yv.row(j);
      
      for (int l=0; l < r; l++) {
        
        ee2 = exp(-txi);
        ee3 = 1 + ee2;
        ee5 = 1.5/ee3 - 1;
        ee6 = exp(lpsi);
        ee7 = y(l) - mu;
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
        
        H(0, 0) += -(ee12 * ee15/(R_pow(ee10, 2) * R_pow(ee6, 2)));
        H(1, 0) += (1 - ee14) * ee12 * ee5/ee11;
        H(2, 0) += -(ee18/ee17);
        H(1, 1) += -((ee14 - 1) * ee12 * ee5 * ee7/ee11);
        H(2, 1) += -(ee18 * ee7/ee17);
        H(2, 3) += -(((((2.25 * (ee7/(ee10 * ee3 * ee6)) - 3) * ee2/
          ee3 + 1.5) * ee12 + 2.25 * (ee2/(ee13 * ee15))) * ee7/ee11 +
            (2.25 * (ee2 * ee7/ee17) - ((4.5/(ee3 * ee5) - 3) * ee2/ee3 +
            1.5) * ee20)/ee15) * ee2/ee13);
    
      }
        
          ee2 = exp(-txi);
          ee3 = 1 + ee2;
          ee5 = 1.5/ee3 - 1;
          ee6 = exp(lpsi);
          ee7 = y(r - 1) - mu;
          ee9 = ee5 * ee7/ee6;
          ee10 = 1/ee5;
          ee11 = ee9 + 1;
          ee12 = 1 + ee10;
          ee13 = R_pow(ee11, ee12);
          ee14 = R_pow(ee3, 2);
          ee15 = log1p(ee9);
          ee16 = R_pow(ee11, (ee10 + 2));
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
          
          H(0, 0) += ee21/(ee16 * R_pow(ee6, 2));
          H(1, 0) += (ee29 - ee33)/ee6;
          H(2, 0) += ee25/ee31;
          H(1, 1) += -((ee33 - ee29) * ee7/ee6);
          H(2, 1) += ee25 * ee7/ee31;
          H(0, 0) += ((((2.25/ee30 - 3) * ee2/ee3 + 1.5)/ee13 - 1.5 *
            (ee25/ee14)) * ee7/ee6 + ((2.25 * (ee2 * ee7/(ee11 * ee14 *
            ee6)) - ((4.5/ee30 - 3) * ee2/ee3 + 1.5) * ee15)/ee20 + 1.5 *
            (ee32 * ee15/ee24))/ee5) * ee2/ee24;
          
      }  
    
H(0, 1) = H(1, 0);
H(0, 2) = H(2, 0);
H(1, 2) = H(2, 1);
  
  return H;
  
}
