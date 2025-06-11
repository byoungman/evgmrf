#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double xieps = 0.0001;

// //' Generalized extreme value (GEV) distribution negative log-likelihood with constrained shape parameter
// //'
// //' @param pars a vector of the three GEV parameters
// //' @param yvec a vector
// //' @return tgevd0 a scalar, the negative log-likelihood
// //' @return tgevd1 a 3-vector
// //' @return tgevd2 a 3x3 matrix
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export(.tgev0)]]
double tgev0(arma::vec pars, arma::vec yv, double delta)
{
    
int n = yv.size(); // number of locations

double nllh = 0.0;

double mu = pars(0);
double lpsi = pars(1);
double txi = pars(2);
double xi = 1.5 / (1.0 + exp(-txi)) - 1.0;
double y;

double ee1, ee2;

for (int i=0; i < n; i++) {
      
  y = yv(i);
  
  if (fabs(xi) >= xieps) {

    ee1 = xi * (y - mu) / exp(lpsi);

    if (ee1 <= -1.0) {

      nllh = 1e20;
      break;

    } else {

      ee2 = 1.0 / xi;
      nllh += lpsi + (ee2 + 1.0) * log1p(ee1) + R_pow(1.0 + ee1, -ee2);

    }

  } else {

    ee1 = (y - mu) / exp(lpsi);
    nllh += lpsi + ee1 + exp(-ee1);

  }

}

nllh += 0.5 * delta * (pars(2) - 0.7) * (pars(2) - 0.7);

return(nllh);

}

// //' @rdname tgevd0
// [[Rcpp::export(.tgev1)]]
arma::vec tgev1(arma::vec pars, arma::vec yv, double delta)
{
  
int n = yv.size(); // number of locations
  
  arma::vec g(3, arma::fill::zeros);
  arma::vec h(6, arma::fill::zeros);
  arma::mat H(3, 3);
  
  double mu = pars(0);
  double lpsi = pars(1);
  double txi = pars(2);
  double y;
  
  double ee2, ee3, ee4, ee5, ee6, ee7, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
  double ee20, ee21, ee22, ee23, ee25, ee27, ee28, ee29;
  double ee30, ee31, ee33, ee34;
  
  for (int i=0; i < n; i++) {
    
    y = yv(i);
    
    ee2 = exp(-txi);
    ee3 = 1 + ee2;
    ee4 = 1.5/ee3;
    ee5 = ee4 - 1;
    ee6 = exp(lpsi);
    ee7 = y - mu;
    ee9 = ee5 * ee7/ee6;
    ee10 = ee9 + 1;
    ee11 = 1/ee5;
    ee12 = R_pow(ee10, ee11);
    ee13 = 1 + ee11;
    ee14 = ee10 * ee6;
    ee15 = R_pow(ee3, 2);
    ee16 = log1p(ee9);
    ee17 = 1/ee12;
    ee18 = R_pow(ee10, ee13);
    ee20 = ee3 * ee5;
    ee21 = 1 + ee17;
    ee22 = 1.5 * (ee16/(ee12 * ee5));
    ee23 = 1.5/ee12;
    ee25 = (((ee23 - 1.5 * ee5) * ee7/ee14 + 1.5) * ee13 - (1.5 +  ee22)/ee5)/ee10 * ee2;
    ee27 = ((4.5/ee20 - 3) * ee2/ee3 + 1.5) * ee16;
    ee28 = ee13 * ee5;
    ee29 = ee13 * ee7;
    ee30 = ee15 * ee6;
    ee31 = R_pow(ee5, 2);
    ee33 = 1.5 * (ee7/(ee18 * ee6));
    ee34 = ee4 - ee21;
    
    g(0) += -((ee28 - ee17)/ee14);
    g(1) += (ee17 - ee28) * ee7/ee14 + 1;
    g(2) += (((ee23 - 1.5) * ee16/ee5 - ee33)/ee5 + 1.5 * (ee29/ee14)) * ee2/ee15;
    
  }
  
  g(2) += delta * (pars(2) - 0.7);
    
  return g;
  
}

// //' @rdname tgevd0
// [[Rcpp::export(.tgev2)]]
arma::mat tgev2(arma::vec pars, arma::vec yv, double delta)
{
  
  int n = yv.size(); // number of locations
  
  arma::vec g(3, arma::fill::zeros);
  arma::vec h(6, arma::fill::zeros);
  arma::mat H(3, 3);
  
  double mu = pars(0);
  double lpsi = pars(1);
  double txi = pars(2);
  double y;
  
  double ee2, ee3, ee4, ee5, ee6, ee7, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
  double ee20, ee21, ee22, ee23, ee25, ee27, ee28, ee29;
  double ee30, ee31, ee33, ee34;
  
  for (int i=0; i < n; i++) {
    
    y = yv(i);
    
    ee2 = exp(-txi);
    ee3 = 1 + ee2;
    ee4 = 1.5/ee3;
    ee5 = ee4 - 1;
    ee6 = exp(lpsi);
    ee7 = y - mu;
    ee9 = ee5 * ee7/ee6;
    ee10 = ee9 + 1;
    ee11 = 1/ee5;
    ee12 = R_pow(ee10, ee11);
    ee13 = 1 + ee11;
    ee14 = ee10 * ee6;
    ee15 = R_pow(ee3, 2);
    ee16 = log1p(ee9);
    ee17 = 1/ee12;
    ee18 = R_pow(ee10, ee13);
    ee20 = ee3 * ee5;
    ee21 = 1 + ee17;
    ee22 = 1.5 * (ee16/(ee12 * ee5));
    ee23 = 1.5/ee12;
    ee25 = (((ee23 - 1.5 * ee5) * ee7/ee14 + 1.5) * ee13 - (1.5 +  ee22)/ee5)/ee10 * ee2;
    ee27 = ((4.5/ee20 - 3) * ee2/ee3 + 1.5) * ee16;
    ee28 = ee13 * ee5;
    ee29 = ee13 * ee7;
    ee30 = ee15 * ee6;
    ee31 = R_pow(ee5, 2);
    ee33 = 1.5 * (ee7/(ee18 * ee6));
    ee34 = ee4 - ee21;
    
    h(0) +=  - (ee13 * ee34 * ee5/(R_pow(ee10, 2) * R_pow(ee6, 2)));
    h(1) += (((ee21 - ee4) * ee7/ee14 + 1) * ee13 * ee5 - ee17)/ee10/ee6;
    h(2) += -(ee25/ee30);
    h(3) += -(((ee34 * ee7/ee14 - 1) * ee13 * ee5 + ee17)/ee10 * ee7/ee6);
    h(4) += -(ee25 * ee7/ee30);
    h(5) += (((((2.25/ee20 - 3) * ee2/ee3 + 1.5)/ee18 - 1.5 * ((1.5 * (ee16/(ee18 * ee31)) -
      1.5 * (ee29/(R_pow(ee10, ee11 +
      2) * ee6))) * ee2/ee15)) * ee7/ee6 + (ee27 + (1.5 * ((ee22 -
      ee33) * ee16/ee5) - 2.25 * (ee7/ee14)) * ee2/ee15 + (2.25 * (ee2 * ee7/(ee10 * ee15 * ee6)) -
      ee27)/ee12)/ee5)/ee5 -
      (((2.25 * (ee7/(ee10 * ee3 * ee6)) - 3) * ee2/ee3 + 1.5) * ee13 +
      2.25 * (ee2/(ee15 * ee31))) * ee7/ee14) * ee2/ee15;    

}

H(0, 0) = h(0);
H(1, 0) = h(1);
H(2, 0) = h(2);
H(1, 1) = h(3);
H(2, 1) = h(4);
H(2, 2) = h(5);
H(0, 1) = H(1, 0);
H(0, 2) = H(2, 0);
H(1, 2) = H(2, 1);

H(2, 2) += delta;

return H;

}

