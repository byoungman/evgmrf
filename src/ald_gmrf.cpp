#define ARMA_USE_SUPERLU 1
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export(.aldgmrfld0)]]
double aldgmrfld0(arma::mat pars, arma::field<arma::vec> yc, double C, double tau)
{
  
  int n = yc.n_rows; // number of locations
  int m;
  
  double nllh = 0.0;
  
  double y, mu, lpsi;
  double res;
  arma::vec yv;

  for (int i=0; i < n; i++) {
    
    mu = pars(0, i);
    lpsi = pars(1, i);

    m = yc(i).size(); // number of obs
    yv = yc(i);
    
    for (int k=0; k < m; k++) {
      
      y = yv(k);
      
      res = (y - mu) / exp(lpsi);
      nllh += lpsi;
      
      if (res <= -C) {
        
        nllh += (tau - 1.0) * (2.0 * res + C);
        
      } else {
        
        if (res < 0.0) {
          
          nllh += (1.0 - tau) * res * res / C;
          
        } else {
          
          if (res <= C) {
            
            nllh += tau * res * res / C;
            
          } else {
            
            nllh += tau * (2.0 * res - C);
            
          }
        }
      }
      
    }
    
  }
  
  return(nllh);
  
}

// [[Rcpp::export(.aldgmrfld1)]]
arma::vec aldgmrfld1(arma::mat pars, arma::field<arma::vec> yc, double C, double tau)
{
  
  int n = yc.n_rows; // number of locations
  int m;
  
  double y, mu, lpsi;
  double res;
  double ee1, ee2, ee3, ee5, ee6, ee7;

  arma::uvec id(2);
  id(0) = 0;
  id(1) = n;
  arma::vec g(n * 2, arma::fill::zeros);
  // arma::vec h(3, arma::fill::zeros);
  // arma::sp_mat H(n * 2, n * 2);
  arma::vec yv;
  
  for (int i=0; i < n; i++) {
    
    mu = pars(0, i);
    lpsi = pars(1, i);
    
    m = yc(i).size(); // number of obs
    yv = yc(i);
    
    // h.zeros();
    
    for (int k=0; k < m; k++) {
      
      y = yv(k);
      res = (y - mu) / exp(lpsi);
      
      if (res <= -C) {
        
        ee1 = exp(lpsi);
        ee2 = tau - 1;
        ee6 = 2 * (ee2 * (y - mu)/ee1);
        ee7 = 2 * (ee2/ee1);
        
        g(id[0]) += -ee7;
        g(id[1]) += 1 - ee6;
        
        // h(0) += 0;
        // h(1) += ee7;
        // h(2) += ee6;
        
      } else {
        
        if (res < 0.0) {
          
          ee1 = 1 - tau;
          ee2 = C * exp(2 * lpsi);
          ee3 = y - mu;
          ee5 = ee1 * ee3/ee2;
          ee7 = ee1 * ee3 * ee3/ee2;
          
          g(id[0]) += -(2 * ee5);
          g(id[1]) += 1 - 2 * ee7;
          
          // h(0) += 2 * (ee1/ee2);
          // h(1) += 4 * ee5;
          // h(2) += 4 * ee7;
          
        } else {
          
          if (res <= C) {
            
            ee1 = C * exp(2 * lpsi);
            ee2 = y - mu;
            ee5 = tau * ee2/ee1;
            ee7 = tau * ee2 * ee2/ee1;
            
            g(id[0]) += -(2 * ee5);
            g(id[1]) += 1 - 2 * ee7;
            
            // h(0) += 2 * (tau/ee1);
            // h(1) += 4 * ee5;
            // h(2) += 4 * ee7;
            
          } else {
            
            ee1 = exp(lpsi);
            ee2 = 2 * (tau * (y - mu)/ee1);
            ee3 = 2 * (tau/ee1);
            
            g(id[0]) += -ee3;
            g(id[1]) += 1 - ee2;
            
            // h(0) += 0;
            // h(1) += ee3;
            // h(2) += ee2;
            
          }
        }
      }
      
    }
    
    // H(id(0), id(0)) = h(0);
    // H(id(1), id(0)) = h(1);
    // H(id(0), id(1)) = h(1);
    // H(id(1), id(1)) = h(2);

    id += 1;
    
  }
  
  return g;
  
}

// [[Rcpp::export(.aldgmrfld2)]]
arma::sp_mat aldgmrfld2(arma::mat pars, arma::field<arma::vec> yc, double C, double tau)
{
  
  int n = yc.n_rows; // number of locations
  int m;
  
  double y, mu, lpsi;
  double res;
  double ee1, ee2, ee3, ee5, ee6, ee7;
  arma::vec yv;
  
  arma::uvec id(2);
  id(0) = 0;
  id(1) = n;
  // arma::vec g(n * 2, arma::fill::zeros);
  arma::vec h(3, arma::fill::zeros);
  arma::sp_mat H(n * 2, n * 2);

  for (int i=0; i < n; i++) {
    
    mu = pars(0, i);
    lpsi = pars(1, i);
    
    m = yc(i).size(); // number of obs
    yv = yc(i);
    
    h.zeros();
    
    for (int k=0; k < m; k++) {
      
      y = yv(k);
      res = (y - mu) / exp(lpsi);
      
      if (res <= -C) {
        
        ee1 = exp(lpsi);
        ee2 = tau - 1;
        ee6 = 2 * (ee2 * (y - mu)/ee1);
        ee7 = 2 * (ee2/ee1);
        
        // g(id[0]) += -ee7;
        // g[id[1]) += 1 - ee6;
        
        h(0) += 0;
        h(1) += ee7;
        h(2) += ee6;
        
      } else {
        
        if (res < 0.0) {
          
          ee1 = 1 - tau;
          ee2 = C * exp(2 * lpsi);
          ee3 = y - mu;
          ee5 = ee1 * ee3/ee2;
          ee7 = ee1 * ee3 * ee3/ee2;
          
          // g(id[0]) += -(2 * ee5);
          // g(id[1]) += 1 - 2 * ee7;
          
          h(0) += 2 * (ee1/ee2);
          h(1) += 4 * ee5;
          h(2) += 4 * ee7;
          
        } else {
          
          if (res <= C) {
            
            ee1 = C * exp(2 * lpsi);
            ee2 = y - mu;
            ee5 = tau * ee2/ee1;
            ee7 = tau * ee2 * ee2/ee1;
            
            // g(id[0]) += -(2 * ee5);
            // g(id[1]) += 1 - 2 * ee7;
            
            h(0) += 2 * (tau/ee1);
            h(1) += 4 * ee5;
            h(2) += 4 * ee7;
            
          } else {
            
            ee1 = exp(lpsi);
            ee2 = 2 * (tau * (y - mu)/ee1);
            ee3 = 2 * (tau/ee1);
            
            // g(id[0]) += -ee3;
            // g(id[1]) += 1 - ee2;
            
            h(0) += 0;
            h(1) += ee3;
            h(2) += ee2;
            
          }
        }
      }
      
    }
    
    H(id(0), id(0)) = h(0);
    H(id(1), id(0)) = h(1);
    H(id(0), id(1)) = h(1);
    H(id(1), id(1)) = h(2);
    
    id += 1;
    
  }
  
  return H;
  
}
