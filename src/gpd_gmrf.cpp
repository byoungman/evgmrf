// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double xieps = 0.0001;

// [[Rcpp::export(.tgpdgmrfld0)]]
double tgpdgmrfld0(arma::mat pars, arma::field<arma::vec> yc)
{
    
int n = yc.n_rows; // number of locations
int m;

double y, lpsi, txi, xi;
double ee1, ee2;
double nllh = 0.0;

arma::vec yv;

for (int i=0; i < n; i++) {
  
  lpsi = pars(0, i);
  txi = pars(1, i);
  xi = 1.5 / (1.0 + exp(-txi)) - 1.0;
  
  m = yc(i).size(); // number of obs
  yv = yc(i);
  
  for (int k=0; k < m; k++) {

    y = yv(k);
    
    if(arma::is_finite(y)) {
      
      ee1 = xi * y / exp(lpsi);
      
      if (ee1 <= -1.0) {
        
        nllh = 1e20;
        break;
        
      } else {
      
        ee2 = 1.0 / xi;
        nllh += lpsi + (ee2 + 1.0) * log1p(ee1);
        
      }
      
    }
    
  } 
  
}


return(nllh);

}

// [[Rcpp::export(.tgpdgmrfld12)]]
arma::mat tgpdgmrfld12(arma::mat pars, arma::field<arma::vec> yc)
{
  
  int n = yc.n_rows; // number of locations
  int m;
  
  double y, pars1, pars2;
  
  arma::mat out = arma::mat(n, 5, arma::fill::zeros);
  
  double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee16, ee17, ee18, ee19;
  
  arma::vec yv;
  
  for (int j=0; j < n; j++) {
    
    pars1 = pars(0, j);
    pars2 = pars(1, j);
    
    m = yc(j).size(); // number of obs
    yv = yc(j);
    
    for (int k=0; k < m; k++) {
      
      y = yv(k);
      
      if(arma::is_finite(y)) {
      
      ee2 = exp(-pars2);
        ee3 = 1 + ee2;
        ee5 = 1.5/ee3 - 1;
        ee6 = exp(pars1);
        ee7 = y * ee5;
        ee8 = ee7/ee6;
        ee9 = 1 + ee8;
        ee10 = ee9 * ee6;
        ee11 = R_pow(ee3, 2);
        ee12 = 1 + 1/ee5;
        ee13 = R_pow(ee5, 2);
        ee14 = y * ee12;
        ee16 = ee11 * ee9 * ee6;
        ee17 = log1p(ee8);
        ee18 = ee14 * ee5;
        ee19 = ee7/ee10;
        
        out(j, 0) += 1 - ee18/ee10;
        out(j, 1) += (1.5 * (ee14/ee10) - 1.5 * (ee17/ee13)) * ee2/
          ee11;
        out(j, 2) += -(ee18 * (ee19 - 1)/ee10);
        out(j, 3) += -(y * (ee12 * (1.5 - 1.5 * ee19) - 1.5/ee5) * ee2/
          ee16);
        out(j, 4) += -(((2.25 * (y * ee2/ee16) - ((4.5/(ee3 * ee5) -
          3) * ee2/ee3 + 1.5) * ee17)/ee13 + y * (((2.25 * (y/(ee3 *
          ee9 * ee6)) - 3) * ee2/ee3 + 1.5) * ee12 + 2.25 * (ee2/(ee11 *
          ee13)))/ee10) * ee2/ee11);
      
      }
      
    }
    
  }
  
  return out;
  
}
