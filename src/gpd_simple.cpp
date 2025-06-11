#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double xieps = 0.0001;

// //' Generalized Pareto distribution (GPD) negative log-likelihood with constrained shape parameter
// //'
// //' @param pars a vector of the two GPD parameters
// //' @param yvec a vector
// //' @return tgpdd0 a scalar, the negative log-likelihood
// //' @return tgpdd1 a 3-vector
// //' @return tgpdd2 a 3x3 matrix
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export(.tgpd0)]]
double tgpd0(arma::vec pars, arma::vec yv)
{
    
int n = yv.size(); // number of locations

double nllh = 0.0;

double pars1 = pars(0);
double pars2 = pars(1);
double xi = 1.5 / (1.0 + exp(-pars2)) - 1.0;
double ixi = 1.0 / xi;
double y;

double ee1;

for (int i=0; i < n; i++) {
      
  y = yv(i);
  
  ee1 = xi * y / exp(pars1);
    
  if (ee1 <= -1.0) {
    nllh = 1e20;
    break;
  } 
    
  nllh += pars1 + (1.0 + ixi) * log1p(ee1);
    
}
  
return(nllh);

}

// //' @rdname tgpdd0
// [[Rcpp::export(.tgpd1)]]
arma::vec tgpd1(arma::vec pars, arma::vec yv)
{
  
int n = yv.size(); // number of locations
  
  arma::vec g(2, arma::fill::zeros);
  arma::vec h(3, arma::fill::zeros);
  arma::mat H(2, 2);
  
  double pars1 = pars(0);
  double pars2 = pars(1);
  double y;

  double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee16, ee17, ee18, ee19;
  
  for (int i=0; i < n; i++) {
    
    y = yv(i);
    
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
    
    g(0) += 1 - ee18/ee10;
    g(1) += (1.5 * (ee14/ee10) - 1.5 * (ee17/ee13)) * ee2/
      ee11;

  }
    
  return g;
  
}

// //' @rdname tgpdd0
// [[Rcpp::export(.tgpd2)]]
arma::mat tgpd2(arma::vec pars, arma::vec yv)
{
  
  int n = yv.size(); // number of locations
  
  arma::vec g(2, arma::fill::zeros);
  arma::vec h(3, arma::fill::zeros);
  arma::mat H(2, 2);
  
  double pars1 = pars(0);
  double pars2 = pars(1);
  double y;
  
  double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee16, ee17, ee18, ee19;
  
  for (int i=0; i < n; i++) {
    
    y = yv(i);
    
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
    
    h(0) += -(ee18 * (ee19 - 1)/ee10);
    h(1) += -(y * (ee12 * (1.5 - 1.5 * ee19) - 1.5/ee5) * ee2/
      ee16);
    h(2) += -(((2.25 * (y * ee2/ee16) - ((4.5/(ee3 * ee5) -
      3) * ee2/ee3 + 1.5) * ee17)/ee13 + y * (((2.25 * (y/(ee3 *
      ee9 * ee6)) - 3) * ee2/ee3 + 1.5) * ee12 + 2.25 * (ee2/(ee11 *
      ee13)))/ee10) * ee2/ee11);
    
  }
  
H(0, 0) = h(0);
H(1, 0) = h(1);
H(1, 1) = h(2);
H(0, 1) = H(1, 0);

return H;

}

