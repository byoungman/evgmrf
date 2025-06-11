// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double xieps = 0.0001;

// //' Generalized extreme value (GEV) distribution negative log-likelihood
// //'
// //' @param pars a nx x ny x 3 cube of GEV parameters
// //' @param y a nx x ny x m cube of block maxima
// //' @return tppugmrfd0 a scalar, the negative log-likelihood
// //' @return tppugmrfd12 a matrix, first then second derivatives w.r.t. GEV parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export(.tppugmrfld0)]]
double tppugmrfld0(arma::mat pars, arma::vec uv, double w)
{
    
int n = uv.size(); // number of locations

double nllh = 0.0;
double u, mu, lpsi, txi, xi, ee1;

for (int j=0; j < n; j++) {
  
  mu = pars(0, j);
  lpsi = pars(1, j);
  txi = pars(2, j);
  xi = 1.5 / (1.0 + exp(-txi)) - 1.0;
  
  u = uv(j);

  ee1 = xi * (u - mu) / exp(lpsi);
  nllh += w * R_pow(1.0 + ee1, -1.0/xi);
  
}

return(nllh);
  
}

// // [[Rcpp::export]]
// double tppugmrfld0(arma::mat pars, arma::vec uv, double w)
// {
//   
//   int n = uv.size(); // number of locations
//   
//   double nllh = 0.0;
//   
//   double u, mu, lpsi, txi, xi;
//   double ee1;
//   
//   for (int i=0; i < n; i++) {
//     
//     mu = pars(0, i);
//     lpsi = pars(1, i);
//     txi = pars(2, i);
//     xi = 1.5 / (1.0 + exp(-txi)) - 1.0;
//     
//     u = uv(i);
//     
//     ee1 = xi * (u - mu) / exp(lpsi);
//     nllh += w * R_pow(1.0 + ee1, -1/xi);
//     
//   }
//   
//   return(nllh);
//   
// }

// // //' @rdname tppugmrfld0
// // [[Rcpp::export]]
// arma::vec tppugmrfld1(arma::mat pars, arma::vec uv, double w)
// {
// 
// int n = uv.size(); // number of locations
//   
// double u, mu, lpsi, txi;
// 
// double ee2, ee3, ee5, ee6, ee7, ee9;
// double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
// double ee20, ee21, ee24, ee25, ee26, ee29;
// double ee30, ee31, ee32, ee33;
// 
// arma::uvec id(3);
// id(0) = 0;
// id(1) = n;
// id(2) = 2 * n;
// arma::vec g(n * 3, arma::fill::zeros);
// 
// for (int i=0; i < n; i++) {
// 
//   mu = pars(0, i);
//   lpsi = pars(1, i);
//   txi = pars(2, i);
//   
//   u = uv(i);
//   
//   ee2 = exp(-txi);
//   ee3 = 1 + ee2;
//   ee5 = 1.5/ee3 - 1;
//   ee6 = exp(lpsi);
//   ee7 = u - mu;
//   ee9 = ee5 * ee7/ee6;
//   ee10 = 1/ee5;
//   ee11 = ee9 + 1;
//   ee12 = 1 + ee10;
//   ee13 = R_pow(ee11, ee12);
//   ee14 = R_pow(ee3, 2);
//   ee15 = log1p(ee9);
//   ee16 = R_pow(ee11, ee10 + 2);
//   ee17 = ee16 * ee6;
//   ee18 = ee13 * ee6;
//   ee20 = R_pow(ee11, ee10);
//   ee21 = ee12 * ee5;
//   ee24 = ee14 * ee5;
//   ee25 = (1.5 * (ee15/(ee13 * R_pow(ee5, 2))) - 1.5 * (ee12 * ee7/ee17)) * ee2;
//   ee26 = ee7/ee18;
//   ee29 = ee21 * ee7/ee17;
//   ee30 = ee3 * ee5;
//   ee31 = ee14 * ee6;
//   ee32 = (1.5 * (ee15/(ee20 * ee5)) - 1.5 * ee26) * ee2;
//   ee33 = 1/ee13;
//   
//   g(id[0]) = 1/ee18;
//   g(id[1]) = ee26;
//   g(id[2]) = ee32/ee24; 
//     
//   id += 1;
// 
// }
// 
// return w * g;
// 
// }
// 
// // //' @rdname tppugmrfld0
// // [[Rcpp::export]]
// arma::sp_mat tppugmrfld2(arma::mat pars, arma::vec uv, double w)
// {
// 
// int n = uv.size(); // number of locations
//   
// arma::sp_mat H(n * 3, n * 3);
// 
// double u, mu, lpsi, txi;
// 
// double ee2, ee3, ee5, ee6, ee7, ee9;
// double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
// double ee20, ee21, ee24, ee25, ee26, ee29;
// double ee30, ee31, ee32, ee33;
// 
// arma::uvec id(3);
// id(0) = 0;
// id(1) = n;
// id(2) = 2 * n;
// arma::vec g(n * 3, arma::fill::zeros);
// 
// arma::vec h(6, arma::fill::zeros);
// 
// for (int i=0; i < n; i++) {
// 
//   mu = pars(0, i);
//   lpsi = pars(1, i);
//   txi = pars(2, i);
// 
//   h.zeros();
// 
//   u = uv(i);
//   
//   ee2 = exp(-txi);
//   ee3 = 1 + ee2;
//   ee5 = 1.5/ee3 - 1;
//   ee6 = exp(lpsi);
//   ee7 = u - mu;
//   ee9 = ee5 * ee7/ee6;
//   ee10 = 1/ee5;
//   ee11 = ee9 + 1;
//   ee12 = 1 + ee10;
//   ee13 = R_pow(ee11, ee12);
//   ee14 = R_pow(ee3, 2);
//   ee15 = log1p(ee9);
//   ee16 = R_pow(ee11, ee10 + 2);
//   ee17 = ee16 * ee6;
//   ee18 = ee13 * ee6;
//   ee20 = R_pow(ee11, ee10);
//   ee21 = ee12 * ee5;
//   ee24 = ee14 * ee5;
//   ee25 = (1.5 * (ee15/(ee13 * R_pow(ee5, 2))) - 1.5 * (ee12 * ee7/ee17)) * ee2;
//   ee26 = ee7/ee18;
//   ee29 = ee21 * ee7/ee17;
//   ee30 = ee3 * ee5;
//   ee31 = ee14 * ee6;
//   ee32 = (1.5 * (ee15/(ee20 * ee5)) - 1.5 * ee26) * ee2;
//   ee33 = 1/ee13;
//   
//   h(0) = ee21/(ee16 * R_pow(ee6, 2));
//   h(1) = (ee29 - ee33)/ee6;
//   h(2) = ee25/ee31;
//   h(3) = -((ee33 - ee29) * ee7/ee6);
//   h(4) = ee25 * ee7/ee31;
//   h(5) = ((((2.25/ee30 - 3) * ee2/ee3 + 1.5)/ee13 - 1.5 * (ee25/ee14)) * ee7/ee6 +
//     ((2.25 * (ee2 * ee7/(ee11 * ee14 * ee6)) -
//     ((4.5/ee30 - 3) * ee2/ee3 + 1.5) * ee15)/ee20 + 1.5 * (ee32 * ee15/ee24))/ee5) * ee2/ee24;
//   
//   H(id(0), id(0)) = h(0);
//   H(id(1), id(0)) = h(1);
//   H(id(2), id(0)) = h(2);
//   H(id(0), id(1)) = h(1);
//   H(id(1), id(1)) = h(3);
//   H(id(2), id(1)) = h(4);
//   H(id(0), id(2)) = h(2);
//   H(id(1), id(2)) = h(4);
//   H(id(2), id(2)) = h(5);
// 
//   id += 1;
// 
// }
// 
// return w * H;
// 
// }

// //' @rdname tppugmrfld0
// [[Rcpp::export(.tppugmrfld12)]]
arma::mat tppugmrfld12(arma::mat pars, arma::vec uv, double w)
{
  
  int n = uv.size(); // number of locations
  
  double u, mu, lpsi, txi;
  
  double ee2, ee3, ee5, ee6, ee7, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
  double ee20, ee21, ee24, ee25, ee26, ee29;
  double ee30, ee31, ee32, ee33;
  
  arma::mat out = arma::mat(n, 9, arma::fill::zeros);
  
  for (int j=0; j < n; j++) {
    
    mu = pars(0, j);
    lpsi = pars(1, j);
    txi = pars(2, j);
    
    u = uv(j);
    
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
    
    out(j, 0) += 1/ee18;
    out(j, 1) += ee26;
    out(j, 2) += ee32/ee24;
    out(j, 3) += ee21/(ee16 * R_pow(ee6, 2));
    out(j, 4) += (ee29 - ee33)/ee6;
    out(j, 5) += ee25/ee31;
    out(j, 6) += -((ee33 - ee29) * ee7/ee6);
    out(j, 7) += ee25 * ee7/ee31;
    out(j, 8) += ((((2.25/ee30 - 3) * ee2/ee3 + 1.5)/ee13 - 1.5 * (ee25/ee14)) * ee7/ee6 +
      ((2.25 * (ee2 * ee7/(ee11 * ee14 * ee6)) -
      ((4.5/ee30 - 3) * ee2/ee3 + 1.5) * ee15)/ee20 + 1.5 * (ee32 * ee15/ee24))/ee5) * ee2/ee24;    
  }
  
  return w * out;
  
}
// //' Generalized extreme value (GEV) distribution negative log-likelihood
// //'
// //' @param pars a nx x ny x 3 cube of GEV parameters
// //' @param y a nx x ny x m cube of block maxima
// //' @return tppzgmrfd0 a scalar, the negative log-likelihood
// //' @return tppzgmrfd12 a matrix, first then second derivatives w.r.t. GEV parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export(.tppzgmrfld0)]]
double tppzgmrfld0(arma::mat pars, arma::field<arma::vec> yc)
{
  
  int n = yc.n_rows; // number of locations
  int m;
  
  double nllh = 0.0;
  
  double y, mu, lpsi, txi, xi;
  double ee1, ee2;
  arma::vec yv;
  double w = 1.0;
  
  for (int i=0; i < n; i++) {
    
    mu = pars(0, i);
    lpsi = pars(1, i);
    txi = pars(2, i);
    xi = 1.5 / (1.0 + exp(-txi)) - 1.0;
    
    m = yc(i).size(); // number of obs
    yv = yc(i);
    
    for (int k=0; k < m; k++) {
      
      y = yv(k);
      
      if(arma::is_finite(y)) {
      
        ee1 = xi * (y - mu) / exp(lpsi);
      
        if (ee1 <= -1.0) {
          nllh = 1e20;
          break;
        } else {
        
          nllh += w * (lpsi + (1 / xi + 1) * log1p(ee1));
        }
     
      }   
       
    }
    
  }
  
  return(nllh);
  
}

// // //' @rdname tppzgmrfld0
// // [[Rcpp::export]]
// arma::vec tppzgmrfld1(arma::mat pars, arma::field<arma::vec> yc)
// {
//   
//   int n = yc.n_rows; // number of locations
//   int m;
//   
//   double y, mu, lpsi, txi;
//   
//   double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
//   double ee10, ee11, ee12, ee13, ee14, ee15, ee17, ee18, ee19;
//   double ee20;
//   
//   arma::uvec id(3);
//   id(0) = 0;
//   id(1) = n;
//   id(2) = 2 * n;
//   arma::vec g(n * 3, arma::fill::zeros);
//   arma::vec yv;
//   
//   for (int i=0; i < n; i++) {
//     
//     mu = pars(0, i);
//     lpsi = pars(1, i);
//     txi = pars(2, i);
//     
//     m = yc(i).size(); // number of obs
//     yv = yc(i);
//     
//     for (int k=0; k < m; k++) {
//       
//       y = yv(k);
//       
//       if(arma::is_finite(y)) {
//       
//       ee2 = exp(-txi);
//       ee3 = 1 + ee2;
//       ee5 = 1.5/ee3 - 1;
//       ee6 = exp(lpsi);
//       ee7 = y - mu;
//       ee8 = ee5 * ee7;
//       ee9 = ee8/ee6;
//       ee10 = ee9 + 1;
//       ee11 = ee10 * ee6;
//       ee12 = 1 + 1/ee5;
//       ee13 = R_pow(ee3, 2);
//       ee14 = ee8/ee11;
//       ee15 = R_pow(ee5, 2);
//       ee17 = ee10 * ee13 * ee6;
//       ee18 = (ee12 * (1.5 - 1.5 * ee14) - 1.5/ee5) * ee2;
//       ee19 = ee12 * ee5;
//       ee20 = log1p(ee9);
//       
//       g(id[0]) += -(ee19/ee11);
//       g(id[1]) += 1 - ee19 * ee7/ee11;
//       g(id[2]) += (1.5 * (ee12 * ee7/ee11) - 1.5 * (ee20/ee15)) * ee2/ee13;
//       
//       }
//   
//     }
//     
//     id += 1;
//     
//   }
//   
//   return g;
//   
// }
// 
// // //' @rdname tppzgmrfld0
// // [[Rcpp::export]]
// arma::sp_mat tppzgmrfld2(arma::mat pars, arma::field<arma::vec> yc)
// {
//   
//   int n = yc.n_rows; // number of locations
//   int m;
//   
//   arma::sp_mat H(n * 3, n * 3);
//   
//   double y, mu, lpsi, txi;
//   
//   double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
//   double ee10, ee11, ee12, ee13, ee14, ee15, ee17, ee18, ee19;
//   double ee20;
//   
//   arma::uvec id(3);
//   id(0) = 0;
//   id(1) = n;
//   id(2) = 2 * n;
//   
//   arma::vec g(n * 3);
//   
//   arma::vec h(6, arma::fill::zeros);
//   arma::vec yv;
//   
//   for (int i=0; i < n; i++) {
//     
//     mu = pars(0, i);
//     lpsi = pars(1, i);
//     txi = pars(2, i);
//     
//     h.zeros();
//     
//     m = yc(i).size(); // number of obs
//     yv = yc(i);
//     
//     for (int k=0; k < m; k++) {
//       
//       y = yv(k);
//       
//       if(arma::is_finite(y)) {
//       
//       ee2 = exp(-txi);
//       ee3 = 1 + ee2;
//       ee5 = 1.5/ee3 - 1;
//       ee6 = exp(lpsi);
//       ee7 = y - mu;
//       ee8 = ee5 * ee7;
//       ee9 = ee8/ee6;
//       ee10 = ee9 + 1;
//       ee11 = ee10 * ee6;
//       ee12 = 1 + 1/ee5;
//       ee13 = R_pow(ee3, 2);
//       ee14 = ee8/ee11;
//       ee15 = R_pow(ee5, 2);
//       ee17 = ee10 * ee13 * ee6;
//       ee18 = (ee12 * (1.5 - 1.5 * ee14) - 1.5/ee5) * ee2;
//       ee19 = ee12 * ee5;
//       ee20 = log1p(ee9);
//       
//       h(0) +=  - (ee12 * ee15/(R_pow(ee10, 2) * R_pow(ee6, 2)));
//       h(1) += (1 - ee14) * ee12 * ee5/ee11;
//       h(2) += -(ee18/ee17);
//       h(3) += -((ee14 - 1) * ee12 * ee5 * ee7/ee11);
//       h(4) += -(ee18 * ee7/ee17);
//       h(5) += -(((((2.25 * (ee7/(ee10 * ee3 * ee6)) - 3) * ee2/ee3 +
//         1.5) * ee12 + 2.25 * (ee2/(ee13 * ee15))) * ee7/ee11 +
//         (2.25 * (ee2 * ee7/ee17) - ((4.5/(ee3 * ee5) - 3) * ee2/ee3 +
//         1.5) * ee20)/ee15) * ee2/ee13);
//       
//       }
//       
//     }
//     
//     H(id(0), id(0)) = h(0);
//     H(id(1), id(0)) = h(1);
//     H(id(2), id(0)) = h(2);
//     H(id(0), id(1)) = h(1);
//     H(id(1), id(1)) = h(3);
//     H(id(2), id(1)) = h(4);
//     H(id(0), id(2)) = h(2);
//     H(id(1), id(2)) = h(4);
//     H(id(2), id(2)) = h(5);
//     
//     id += 1;
//     
//   }
//   
//   return H;
//   
// }

// //' @rdname tppzgmrfld0
// [[Rcpp::export(.tppzgmrfld12)]]
arma::mat tppzgmrfld12(arma::mat pars, arma::field<arma::vec> yc)
{
  
  int n = yc.n_rows; // number of locations
  int m;
  
  double y, mu, lpsi, txi, xi;
  
  arma::mat out = arma::mat(n, 9, arma::fill::zeros);
  
  double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee17, ee18, ee19;
  double ee20;
  
  arma::vec yv;
  
  for (int j=0; j < n; j++) {
    
    mu = pars(0, j);
    lpsi = pars(1, j);
    txi = pars(2, j);
    
    m = yc(j).size(); // number of obs
    yv = yc(j);
    
    for (int k=0; k < m; k++) {
      
      y = yv(k);
      
      if(arma::is_finite(y)) {
        
        ee2 = exp(-txi);
        ee3 = 1 + ee2;
        ee5 = 1.5/ee3 - 1;
        ee6 = exp(lpsi);
        ee7 = y - mu;
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
        
        out(j, 0) += -(ee19/ee11);
        out(j, 1) += 1 - ee19 * ee7/ee11;
        out(j, 2) += (1.5 * (ee12 * ee7/ee11) - 1.5 * (ee20/ee15)) * ee2/ee13;

        out(j, 3) +=  - (ee12 * ee15/(R_pow(ee10, 2) * R_pow(ee6, 2)));
        out(j, 4) +=  (1 - ee14) * ee12 * ee5/ee11;
        out(j, 5) +=  -(ee18/ee17);
        out(j, 6) +=  -((ee14 - 1) * ee12 * ee5 * ee7/ee11);
        out(j, 7) +=  -(ee18 * ee7/ee17);
        out(j, 8) +=  -(((((2.25 * (ee7/(ee10 * ee3 * ee6)) - 3) * ee2/ee3 +
          1.5) * ee12 + 2.25 * (ee2/(ee13 * ee15))) * ee7/ee11 +
          (2.25 * (ee2 * ee7/ee17) - ((4.5/(ee3 * ee5) - 3) * ee2/ee3 +
          1.5) * ee20)/ee15) * ee2/ee13);
        
      }
    
    }
    
  }
  
  return out;
  
}

