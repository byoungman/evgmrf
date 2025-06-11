// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export(.tpprocgmrfld0)]]
double tpprocgmrfld0(arma::mat pars, arma::field<double> yc, double mult)
{
    
int n = yc.n_rows; // number of locations

double y, llambda, lambda;
double nllh = 0.0;

double yv;

for (int i=0; i < n; i++) {
  
  llambda = pars(0, i);
  lambda = exp(llambda);
  yv = yc(i);
  
  nllh += mult * lambda - llambda * yv;
        
}

return(nllh);

}

// [[Rcpp::export(.tpprocgmrfld12)]]
arma::mat tpprocgmrfld12(arma::mat pars, arma::field<double> yc, double mult)
{
  
  int n = yc.n_rows; // number of locations

  double y, llambda, lambda;
  arma::mat out = arma::mat(n, 2, arma::fill::zeros);
  
  double yv;
  
  for (int i=0; i < n; i++) {
    
    llambda = pars(0, i);
    lambda = exp(llambda);
    yv = yc(i);
    
    out(i, 0) += mult * lambda - yv;
    out(i, 1) += mult * lambda;

  }
  
  return out;
  
}

// // [[Rcpp::export(.tpprocgmrfld0t)]]
// double tpprocgmrfld0t(arma::mat pars, arma::field<arma::vec> yc, double mult)
// {
//   
//   int n = yc.n_rows; // number of locations
//   
//   double y, llambda, lambda;
//   double nllh = 0.0;
//   int m;
//   
//   arma::vec yv;
//   
//   for (int i=0; i < n; i++) {
//     
//     m = yc(i).size(); // number of obs
//     yv = yc(i);
//     
//     for (int k=0; k < m; k++) {
//       
//       llambda = pars(i, k);
//       lambda = exp(llambda);
//       
//       nllh += mult * lambda - llambda * yv(k);
//       
//     }
//     
//   }
//   
//   return(nllh);
//   
// }
// 
// // [[Rcpp::export(.tpprocgmrfld12t)]]
// arma::mat tpprocgmrfld12t(arma::mat pars, arma::field<arma::vec> yc, double mult)
// {
//   
//   int n = yc.n_rows; // number of locations
//   
//   double y, llambda, lambda;
//   
//   arma::vec yv;
//   int m = yc(0).size(); // number of obs
//   
//   arma::mat out = arma::mat(n * m, 2, arma::fill::zeros);
// 
//   for (int i=0; i < n; i++) {
//     
//     llambda = pars(0, i);
//     lambda = exp(llambda);
//     yv = yc(i);
//     
//     for (int k=0; k < m; k++) {
//       
//       llambda = pars(i, k);
//       lambda = exp(llambda);
//       
//       out(k * n + i, 0) += mult * lambda - yv(k);
//       out(k * n + i, 1) += mult * lambda;
//       
//     }
//     
//   }
//   
//   return out;
//   
// }
