// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double xieps = 0.0001;

// [[Rcpp::export(.rlargegmrfld0)]]
double rlargegmrfld0(arma::mat pars, arma::field<arma::mat> yc, arma::field<arma::mat> wc)
{
  
int n = yc.n_rows; // number of locations
int m, r;
    
double nllh = 0.0;

double mu, lpsi, txi, xi;
double ee2;
arma::mat yv, wv;
arma::rowvec y, w, ee1;

for (int j=0; j < n; j++) {
      
  mu = pars(0, j);
  lpsi = pars(1, j);
  txi = pars(2, j);
  xi = 1.5 / (1.0 + exp(-txi)) - 1.0;
  
  m = yc(j).n_rows; // number of obs
  r = yc(j).n_cols; // number of order statistics
  yv = yc(j);
  wv = wc(j);
  
  for (int k=0; k < m; k++) {
    
    y = yv.row(k);
    w = wv.row(k);
    
    if (arma::is_finite(y)) {
      
    ee1 = xi * (y - mu) / exp(lpsi);
    
    for (int l=0; l < r; l++) {
      
      ee2 = 1.0 / xi;
      
      if (ee1(l) <= -1.0) {
        
        nllh = 1e20;
        break;
        
      } else {
        
        nllh += w(l) * (lpsi + (ee2 + 1.0) * log1p(ee1(l)));
        
      }

    }
    
    nllh += w[r - 1] * R_pow(1 + ee1[r - 1], -ee2);
    
    }

  }
  
}

return(nllh);

}

// // [[Rcpp::export]]
// arma::vec rlargegmrfld1(arma::mat pars, arma::field<arma::mat> yc)
// {
// 
// int n = yc.n_rows; // number of locations
// int m, r;
// 
// arma::uvec id = {0, 1, 2};
// arma::vec g(n * 3, arma::fill::zeros);
// 
// double mu, lpsi, txi, xi;
// 
// arma::mat yv;
// arma::rowvec y, ee1;
// 
// double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
// double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
// double ee20, ee21, ee24, ee25, ee26, ee29;
// double ee30, ee31, ee32, ee33;
// 
// for (int i=0; i < n; i++) {
//   
//   mu = pars(0, i);
//   lpsi = pars(1, i);
//   txi = pars(2, i);
//   
//   m = yc(i).n_rows; // number of obs
//   r = yc(i).n_cols; // number of order statistics
//   yv = yc(i);
//   
//   for (int k=0; k < m; k++) {
//     
//     y = yv.row(k);
//     
//     if(arma::is_finite(y)) {
//     
//     for (int l=0; l < r; l++) {
// 
//       ee2 = exp(-txi);
//       ee3 = 1 + ee2;
//       ee5 = 1.5/ee3 - 1;
//       ee6 = exp(lpsi);
//       ee7 = y(l) - mu;
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
//       g(i) += -(ee19/ee11);
//       g(i + n) += 1 - ee19 * ee7/ee11;
//       g(i + 2 * n) += (1.5 * (ee12 * ee7/ee11) - 1.5 * (ee20/ee15)) *
//         ee2/ee13;
// 
//     }
//     
//       ee2 = exp(-txi);
//       ee3 = 1 + ee2;
//       ee5 = 1.5/ee3 - 1;
//       ee6 = exp(lpsi);
//       ee7 = y[r - 1] - mu;
//       ee9 = ee5 * ee7/ee6;
//       ee10 = 1/ee5;
//       ee11 = ee9 + 1;
//       ee12 = 1 + ee10;
//       ee13 = R_pow(ee11, ee12);
//       ee14 = R_pow(ee3, 2);
//       ee15 = log1p(ee9);
//       ee16 = R_pow(ee11, (ee10 + 2));
//       ee17 = ee16 * ee6;
//       ee18 = ee13 * ee6;
//       ee20 = R_pow(ee11, ee10);
//       ee21 = ee12 * ee5;
//       ee24 = ee14 * ee5;
//       ee25 = (1.5 * (ee15/(ee13 * R_pow(ee5, 2))) - 1.5 * (ee12 * ee7/ee17)) * ee2;
//       ee26 = ee7/ee18;
//       ee29 = ee21 * ee7/ee17;
//       ee30 = ee3 * ee5;
//       ee31 = ee14 * ee6;
//       ee32 = (1.5 * (ee15/(ee20 * ee5)) - 1.5 * ee26) * ee2;
//       ee33 = 1/ee13;
// 
//       g(i) += 1/ee18;
//       g(i + n) += ee26;
//       g(i + 2 * n) += ee32/ee24;
//       
//   }
//     
//   }
//   
//   id += 1;
// 
// }
// 
// return g;
// 
// }

// [[Rcpp::export(.rlargegmrfld12)]]
arma::mat rlargegmrfld12(arma::mat pars, arma::field<arma::mat> yc, arma::field<arma::mat> wc)
{
  
  int n = yc.n_rows; // number of locations
  int m, r;
  
  double mu, lpsi, txi, xi;
  
  arma::mat yv, wv;
  arma::rowvec y, w, ee1;
  
  double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21, ee24, ee25, ee26, ee29, ee30, ee31, ee32, ee33;
  
  arma::mat out = arma::mat(n, 9, arma::fill::zeros);
  
  for (int j=0; j < n; j++) {
    
    mu = pars(0, j);
    lpsi = pars(1, j);
    txi = pars(2, j);
    xi = 1.5 / (1.0 + exp(-txi)) - 1.0;
    
    m = yc(j).n_rows; // number of obs
    r = yc(j).n_cols; // number of order statistics
    yv = yc(j);
    wv = wc(j);
    
    for (int k=0; k < m; k++) {
      
      y = yv.row(k);
      w = wv.row(k);
      
      if(arma::is_finite(y)) {
      
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
        
        out(j, 0) += w(l) * (-(ee19/ee11));
        out(j, 1) += w(l) * (1 - ee19 * ee7/ee11);
        out(j, 2) += w(l) * ((1.5 * (ee12 * ee7/ee11) - 1.5 * (ee20/ee15)) *
          ee2/ee13);
        out(j, 3) += w(l) * (-(ee12 * ee15/(R_pow(ee10, 2) * R_pow(ee6, 2))));
        out(j, 4) += w(l) * ((1 - ee14) * ee12 * ee5/ee11);
        out(j, 5) += w(l) * (-(ee18/ee17));
        out(j, 6) += w(l) * (-((ee14 - 1) * ee12 * ee5 * ee7/ee11));
        out(j, 7) += w(l) * (-(ee18 * ee7/ee17));
        out(j, 8) += w(l) * (-(((((2.25 * (ee7/(ee10 * ee3 * ee6)) - 3) * ee2/
          ee3 + 1.5) * ee12 + 2.25 * (ee2/(ee13 * ee15))) * ee7/ee11 +
            (2.25 * (ee2 * ee7/ee17) - ((4.5/(ee3 * ee5) - 3) * ee2/ee3 +
            1.5) * ee20)/ee15) * ee2/ee13));
    
      }
        
          ee2 = exp(-txi);
          ee3 = 1 + ee2;
          ee5 = 1.5/ee3 - 1;
          ee6 = exp(lpsi);
          ee7 = y[r - 1] - mu;
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
          
          out(j, 0) += w[r - 1] * (1/ee18);
          out(j, 1) += w[r - 1] * (ee26);
          out(j, 2) += w[r - 1] * (ee32/ee24);
          out(j, 3) += w[r - 1] * (ee21/(ee16 * R_pow(ee6, 2)));
          out(j, 4) += w[r - 1] * ((ee29 - ee33)/ee6);
          out(j, 5) += w[r - 1] * (ee25/ee31);
          out(j, 6) += w[r - 1] * (-((ee33 - ee29) * ee7/ee6));
          out(j, 7) += w[r - 1] * (ee25 * ee7/ee31);
          out(j, 8) += w[r - 1] * (((((2.25/ee30 - 3) * ee2/ee3 + 1.5)/ee13 - 1.5 *
            (ee25/ee14)) * ee7/ee6 + ((2.25 * (ee2 * ee7/(ee11 * ee14 *
            ee6)) - ((4.5/ee30 - 3) * ee2/ee3 + 1.5) * ee15)/ee20 + 1.5 *
            (ee32 * ee15/ee24))/ee5) * ee2/ee24);
          
      }
      
    }
    
  }
  
  return out;
  
}

// // [[Rcpp::export]]
// arma::sp_mat rlargegmrfld2(arma::mat pars, arma::field<arma::mat> yc)
// {
// 
// int n = yc.n_rows; // number of locations
// int m, r;
// 
// arma::sp_mat H(n * 3, n * 3);
// arma::uvec id = {0, 1, 2};
// arma::vec h(6, arma::fill::zeros);
// 
// double mu, lpsi, txi, xi;
// 
// arma::mat yv;
// arma::rowvec y, ee1;
// 
// double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
// double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
// double ee20, ee21, ee24, ee25, ee26, ee29;
// double ee30, ee31, ee32, ee33;
// 
// for (int i=0; i < n; i++) {
//   
//   mu = pars(0, i);
//   lpsi = pars(1, i);
//   txi = pars(2, i);
//   
//   h.zeros();
//   
//   m = yc(i).n_rows; // number of obs
//   r = yc(i).n_cols; // number of order statistics
//   yv = yc(i);
//   
//   for (int k=0; k < m; k++) {
//     
//     y = yv.row(k);
//     
//     if(arma::is_finite(y)) {
// 
//     for (int l=0; l < r; l++) {
// 
//       ee2 = exp(-txi);
//       ee3 = 1 + ee2;
//       ee5 = 1.5/ee3 - 1;
//       ee6 = exp(lpsi);
//       ee7 = y(l) - mu;
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
//       h(0) += -(ee12 * ee15/(R_pow(ee10, 2) * R_pow(ee6, 2)));
//       h(1) += (1 - ee14) * ee12 * ee5/ee11;
//       h(2) += -(ee18/ee17);
//       h(3) += -((ee14 - 1) * ee12 * ee5 * ee7/ee11);
//       h(4) += -(ee18 * ee7/ee17);
//       h(5) += -(((((2.25 * (ee7/(ee10 * ee3 * ee6)) - 3) * ee2/
//         ee3 + 1.5) * ee12 + 2.25 * (ee2/(ee13 * ee15))) * ee7/ee11 +
//           (2.25 * (ee2 * ee7/ee17) - ((4.5/(ee3 * ee5) - 3) * ee2/
//             ee3 + 1.5) * ee20)/ee15) * ee2/ee13);
//     
//     }
// 
//     ee2 = exp(-txi);
//     ee3 = 1 + ee2;
//     ee5 = 1.5/ee3 - 1;
//     ee6 = exp(lpsi);
//     ee7 = y[r - 1] - mu;
//     ee9 = ee5 * ee7/ee6;
//     ee10 = 1/ee5;
//     ee11 = ee9 + 1;
//     ee12 = 1 + ee10;
//     ee13 = R_pow(ee11, ee12);
//     ee14 = R_pow(ee3, 2);
//     ee15 = log1p(ee9);
//     ee16 = R_pow(ee11, (ee10 + 2));
//     ee17 = ee16 * ee6;
//     ee18 = ee13 * ee6;
//     ee20 = R_pow(ee11, ee10);
//     ee21 = ee12 * ee5;
//     ee24 = ee14 * ee5;
//     ee25 = (1.5 * (ee15/(ee13 * R_pow(ee5, 2))) - 1.5 * (ee12 * ee7/ee17)) * ee2;
//     ee26 = ee7/ee18;
//     ee29 = ee21 * ee7/ee17;
//     ee30 = ee3 * ee5;
//     ee31 = ee14 * ee6;
//     ee32 = (1.5 * (ee15/(ee20 * ee5)) - 1.5 * ee26) * ee2;
//     ee33 = 1/ee13;
//     
//     h(0) += ee21/(ee16 * R_pow(ee6, 2));
//     h(1) += (ee29 - ee33)/ee6;
//     h(2) += ee25/ee31;
//     h(3) += -((ee33 - ee29) * ee7/ee6);
//     h(4) += ee25 * ee7/ee31;
//     h(5) += ((((2.25/ee30 - 3) * ee2/ee3 + 1.5)/ee13 - 1.5 *
//       (ee25/ee14)) * ee7/ee6 + ((2.25 * (ee2 * ee7/(ee11 * ee14 *
//       ee6)) - ((4.5/ee30 - 3) * ee2/ee3 + 1.5) * ee15)/ee20 + 1.5 *
//       (ee32 * ee15/ee24))/ee5) * ee2/ee24;
//     
//     }
//     
//   }
//         
//   H(i, i) = h(0);
//   H(i + n, i) = h(1);
//   H(i, i + n) = h(1);
//   H(i + 2 * n, i) = h(2);
//   H(i, i + 2 * n) = h(2);
//   H(i + n, i + n) = h(3);
//   H(i + n, i + 2 * n) = h(4);
//   H(i + 2 * n, i + n) = h(4);
//   H(i + 2 * n, i + 2 * n) = h(5);
//   
//   id += 3;
// 
// }
// 
// return H;
// 
// }
// 
// // [[Rcpp::export]]
// arma::mat rlargegmrfld2mat(arma::mat pars, arma::field<arma::mat> yc)
// {
//   
//   int n = yc.n_rows; // number of locations
//   int m, r;
//   
//   arma::mat H(n, 6, arma::fill::zeros);
// 
//   double mu, lpsi, txi, xi;
//   
//   arma::mat yv;
//   arma::rowvec y, ee1;
//   
//   double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
//   double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
//   double ee20, ee21, ee24, ee25, ee26, ee29;
//   double ee30, ee31, ee32, ee33;
//   
//   for (int i=0; i < n; i++) {
//     
//     mu = pars(0, i);
//     lpsi = pars(1, i);
//     txi = pars(2, i);
//     
//     m = yc(i).n_rows; // number of obs
//     r = yc(i).n_cols; // number of order statistics
//     yv = yc(i);
//     
//     for (int k=0; k < m; k++) {
//       
//       y = yv.row(k);
//       
//       if(arma::is_finite(y)) {
//       
//       for (int l=0; l < r; l++) {
// 
//         ee2 = exp(-txi);
//         ee3 = 1 + ee2;
//         ee5 = 1.5/ee3 - 1;
//         ee6 = exp(lpsi);
//         ee7 = y(l) - mu;
//         ee8 = ee5 * ee7;
//         ee9 = ee8/ee6;
//         ee10 = ee9 + 1;
//         ee11 = ee10 * ee6;
//         ee12 = 1 + 1/ee5;
//         ee13 = R_pow(ee3, 2);
//         ee14 = ee8/ee11;
//         ee15 = R_pow(ee5, 2);
//         ee17 = ee10 * ee13 * ee6;
//         ee18 = (ee12 * (1.5 - 1.5 * ee14) - 1.5/ee5) * ee2;
//         ee19 = ee12 * ee5;
//         ee20 = log1p(ee9);
//         
//         H(i, 0) += -(ee12 * ee15/(R_pow(ee10, 2) * R_pow(ee6, 2)));
//         H(i, 1) += (1 - ee14) * ee12 * ee5/ee11;
//         H(i, 2) += -(ee18/ee17);
//         H(i, 3) += -((ee14 - 1) * ee12 * ee5 * ee7/ee11);
//         H(i, 4) += -(ee18 * ee7/ee17);
//         H(i, 5) += -(((((2.25 * (ee7/(ee10 * ee3 * ee6)) - 3) * ee2/
//           ee3 + 1.5) * ee12 + 2.25 * (ee2/(ee13 * ee15))) * ee7/ee11 +
//             (2.25 * (ee2 * ee7/ee17) - ((4.5/(ee3 * ee5) - 3) * ee2/
//               ee3 + 1.5) * ee20)/ee15) * ee2/ee13);
// 
//       }
//       
//           ee2 = exp(-txi);
//           ee3 = 1 + ee2;
//           ee5 = 1.5/ee3 - 1;
//           ee6 = exp(lpsi);
//           ee7 = y[r - 1] - mu;
//           ee9 = ee5 * ee7/ee6;
//           ee10 = 1/ee5;
//           ee11 = ee9 + 1;
//           ee12 = 1 + ee10;
//           ee13 = R_pow(ee11, ee12);
//           ee14 = R_pow(ee3, 2);
//           ee15 = log1p(ee9);
//           ee16 = R_pow(ee11, (ee10 + 2));
//           ee17 = ee16 * ee6;
//           ee18 = ee13 * ee6;
//           ee20 = R_pow(ee11, ee10);
//           ee21 = ee12 * ee5;
//           ee24 = ee14 * ee5;
//           ee25 = (1.5 * (ee15/(ee13 * R_pow(ee5, 2))) - 1.5 * (ee12 * ee7/ee17)) * ee2;
//           ee26 = ee7/ee18;
//           ee29 = ee21 * ee7/ee17;
//           ee30 = ee3 * ee5;
//           ee31 = ee14 * ee6;
//           ee32 = (1.5 * (ee15/(ee20 * ee5)) - 1.5 * ee26) * ee2;
//           ee33 = 1/ee13;
// 
//           H(i, 0) += ee21/(ee16 * R_pow(ee6, 2));
//           H(i, 1) += (ee29 - ee33)/ee6;
//           H(i, 2) += ee25/ee31;
//           H(i, 3) += -((ee33 - ee29) * ee7/ee6);
//           H(i, 4) += ee25 * ee7/ee31;
//           H(i, 5) += ((((2.25/ee30 - 3) * ee2/ee3 + 1.5)/ee13 - 1.5 *
//             (ee25/ee14)) * ee7/ee6 + ((2.25 * (ee2 * ee7/(ee11 * ee14 *
//             ee6)) - ((4.5/ee30 - 3) * ee2/ee3 + 1.5) * ee15)/ee20 + 1.5 *
//             (ee32 * ee15/ee24))/ee5) * ee2/ee24;
//           
//     }
//       
//     }
//     
//   }
//     
//   return H;
//   
// }

