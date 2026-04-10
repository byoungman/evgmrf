// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export(.aldgmrfld0)]]
double aldgmrfld0(arma::mat pars, arma::field<arma::vec> yc, arma::field<arma::vec> wc, double C, double tau)
{
  
  int n = yc.n_rows; // number of locations
  int m;
  
  double nllh = 0.0;
  
  double y, w, mu, lpsi;
  double res;
  arma::vec yv, wv;

  for (int j=0; j < n; j++) {
    
    mu = pars(0, j);
    lpsi = pars(1, j);

    m = yc(j).size(); // number of obs
    yv = yc(j);
    wv = wc(j);
    
    for (int k=0; k < m; k++) {
      
      y = yv(k);
      w = wv(k);
      
      res = (y - mu) / exp(lpsi);
      nllh += w * lpsi;
      
      if (res <= -C) {
        
        nllh += w * (tau - 1.0) * (2.0 * res + C);
        
      } else {
        
        if (res < 0.0) {
          
          nllh += w * (1.0 - tau) * res * res / C;
          
        } else {
          
          if (res <= C) {
            
            nllh += w * tau * res * res / C;
            
          } else {
            
            nllh += w * tau * (2.0 * res - C);
            
          }
        }
      }
      
    }
    
  }
  
  return(nllh);
  
}

// [[Rcpp::export(.aldgmrfld12)]]
arma::mat aldgmrfld12(arma::mat pars, arma::field<arma::vec> yc, arma::field<arma::vec> wc, double C, double tau)
{
  
  int n = yc.n_rows; // number of locations
  int m;
  
  double y, w, mu, lpsi;
  double res;
  double ee1, ee2, ee3, ee5, ee6, ee7;
  arma::vec yv, wv;
  
  arma::mat out = arma::mat(n, 5, arma::fill::zeros);

  for (int j=0; j < n; j++) {
    
    mu = pars(0, j);
    lpsi = pars(1, j);
    
    m = yc(j).size(); // number of obs
    yv = yc(j);
    wv = wc(j);
    
    for (int k=0; k < m; k++) {
      
      y = yv(k);
      w = wv(k);
      res = (y - mu) / exp(lpsi);
      
      if (res <= -C) {
        
        ee1 = exp(lpsi);
        ee2 = tau - 1;
        ee6 = 2 * (ee2 * (y - mu)/ee1);
        ee7 = 2 * (ee2/ee1);
        
        out(j, 0) += w * (-ee7);
        out(j, 1) += w * (1 - ee6);
        out(j, 2) += w * (0);
        out(j, 3) += w * (ee7);
        out(j, 4) += w * (ee6);
        
      } else {
        
        if (res < 0.0) {
          
          ee1 = 1 - tau;
          ee2 = C * exp(2 * lpsi);
          ee3 = y - mu;
          ee5 = ee1 * ee3/ee2;
          ee7 = ee1 * ee3 * ee3/ee2;
          
          out(j, 0) += w * (-(2 * ee5));
          out(j, 1) += w * (1 - 2 * ee7);
          out(j, 2) += w * (2 * (ee1/ee2));
          out(j, 3) += w * (4 * ee5);
          out(j, 4) += w * (4 * ee7);
          
        } else {
          
          if (res <= C) {
            
            ee1 = C * exp(2 * lpsi);
            ee2 = y - mu;
            ee5 = tau * ee2/ee1;
            ee7 = tau * ee2 * ee2/ee1;
            
            out(j, 0) += w * (-(2 * ee5));
            out(j, 1) += w * (1 - 2 * ee7);
            out(j, 2) += w * (2 * (tau/ee1));
            out(j, 3) += w * (4 * ee5);
            out(j, 4) += w * (4 * ee7);
            
          } else {
            
            ee1 = exp(lpsi);
            ee2 = 2 * (tau * (y - mu)/ee1);
            ee3 = 2 * (tau/ee1);
            
            out(j, 0) += w * (-ee3);
            out(j, 1) += w * (1 - ee2);
            out(j, 2) += w * (0);
            out(j, 3) += w * (ee3);
            out(j, 4) += w * (ee2);
            
          }
        }
      }
      
    }
    
  }
  
  return out;
  
}

// [[Rcpp::export(.aldgmrfldJ)]]
arma::mat aldgmrfldJ(arma::mat pars, arma::field<arma::vec> yc, arma::field<arma::vec> wc, double C, double tau)
{
  
  int n = yc.n_rows; // number of locations
  int m;
  
  double y, w, mu, lpsi;
  double res;
  double ee1, ee2, ee3, ee5, ee6, ee7;
  arma::vec yv, wv;
  
  arma::uvec id(2);
  id(0) = 0;
  id(1) = n;
  
  arma::mat h = arma::mat(n, 3, arma::fill::zeros);
  arma::cube g = arma::cube(n, m, 2, arma::fill::zeros);

  for (int j=0; j < n; j++) {
    
    mu = pars(0, j);
    lpsi = pars(1, j);
    
    m = yc(j).size(); // number of obs
    yv = yc(j);
    wv = wc(j);
    
    h.zeros();
    
    for (int k=0; k < m; k++) {
      
      y = yv(k);
      w = wv(k);
      res = (y - mu) / exp(lpsi);
      
      if (res <= -C) {
        
        ee1 = exp(lpsi);
        ee2 = tau - 1;
        ee6 = 2 * (ee2 * (y - mu)/ee1);
        ee7 = 2 * (ee2/ee1);
        
        g(j, k, 0) = w * (-ee7);
        g(j, k, 1) = w * (1 - ee6);
        h(j, 0) += w * (0);
        h(j, 1) += w * (ee7);
        h(j, 2) += w * (ee6);
        
      } else {
        
        if (res < 0.0) {
          
          ee1 = 1 - tau;
          ee2 = C * exp(2 * lpsi);
          ee3 = y - mu;
          ee5 = ee1 * ee3/ee2;
          ee7 = ee1 * ee3 * ee3/ee2;
          
          g(j, k, 0) = w * (-(2 * ee5));
          g(j, k, 1) = w * (1 - 2 * ee7);
          h(j, 2) += w * (2 * (ee1/ee2));
          h(j, 3) += w * (4 * ee5);
          h(j, 4) += w * (4 * ee7);
          
        } else {
          
          if (res <= C) {
            
            ee1 = C * exp(2 * lpsi);
            ee2 = y - mu;
            ee5 = tau * ee2/ee1;
            ee7 = tau * ee2 * ee2/ee1;
            
            g(j, k, 0) = w * (-(2 * ee5));
            g(j, k, 1) = w * (1 - 2 * ee7);
            h(j, 2) += w * (2 * (tau/ee1));
            h(j, 3) += w * (4 * ee5);
            h(j, 4) += w * (4 * ee7);
            
          } else {
            
            ee1 = exp(lpsi);
            ee2 = 2 * (tau * (y - mu)/ee1);
            ee3 = 2 * (tau/ee1);
            
            g(j, k, 0) = w * (-ee3);
            g(j, k, 1) = w * (1 - ee2);
            h(j, 2) += w * (0);
            h(j, 3) += w * (ee3);
            h(j, 4) += w * (ee2);
            
          }
        }
      }
      
    }
    
  }
  
  arma::mat G = arma::mat(2 * n, 2 * n, arma::fill::zeros);
  arma::mat Gk;
  arma::vec Gk2;
    
  for (int k=0; k < m; k++) {  
    Gk = g.col(k);
    Gk2 = arma::vectorise(Gk);
    G += Gk2 * Gk2.t();
  }
    
  return G;

}


