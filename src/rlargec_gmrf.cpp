// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double xieps = 0.0001;

// [[Rcpp::export(.rlargecgmrfld0)]]
double rlargecgmrfld0(arma::mat pars, arma::field<arma::mat> yc, int drop)
{
  
  int n = yc.n_rows; // number of locations
  int m;
  int r = yc(0).n_cols - drop;
  // Rcpp::Rcout << "The value of r is: " << r << std::endl;
  
  double nllh = 0.0;
  
  double mu, lpsi, txi, xi;
  arma::mat yv;
  arma::rowvec y(r + drop), z(r + drop), w(r + drop);
  
  for (int i=0; i < n; i++) {
    
    mu = pars(0, i);
    lpsi = pars(1, i);
    txi = pars(2, i);
    xi = 1.5 / (1.0 + exp(-txi)) - 1.0;
    
    m = yc(i).n_rows; // number of obs
    yv = yc(i);
    
    for (int k=0; k < m; k++) {
      
      y = yv.row(k);
      z = 1 + xi * (y - mu) / exp(lpsi);
      
      for (int l=0; l < r + drop; l++) {
        w(l) = R_pow(z(l), -1.0/xi);
        if (l < drop) {
          if (z(l) < 0.0) {
            w(l) = 0.0;
          }
        } else {
          if (z(l) < 0.0) {
            nllh = 1e20;
            break;
          }
        }
      }
      
      for (int l=0; l < r; l++) {
        
        if (l == r - 1) {
          
          nllh -= log(exp(-w(l)) - exp(-w(l + drop)));
          
        } else {
          
          if (l < drop && z(l) <= 0.0) {
            
            nllh += (1 / xi) * log(z(l + drop));
            
          } else {
            
            nllh -= log(w(l + drop) - w(l));
            
          }
          
        }
        
      } 
      
    }
    
  }
  
  return(nllh);
  
}

// [[Rcpp::export(.rlargecgmrfld12)]]
arma::mat rlargecgmrfld12(arma::mat pars, arma::field<arma::mat> yc, int drop)
{
  
  int n = yc.n_rows; // number of locations
  int m;
  int r = yc(0).n_cols - drop;
  
  arma::mat yv;
  
  double mu, lpsi, txi, xi;
  arma::rowvec y(r + drop), z(r + drop);
  
  double ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9, ee10;
  double ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
  double ee34, ee35, ee36, ee37, ee38, ee39, ee30, ee31, ee32, ee33;
  double ee43, ee45, ee46, ee47, ee48, ee40, ee41, ee42;
  double ee50, ee56, ee58, ee51, ee53, ee59;
  double ee63, ee64, ee66, ee67, ee68, ee60, ee61, ee65;
  double ee70, ee72, ee73, ee74, ee75, ee76, ee78, ee79;
  double ee80, ee83, ee86;
  
  arma::mat out = arma::mat(n, 9, arma::fill::zeros);
  
  for (int j=0; j < n; j++) {
    
    mu = pars(0, j);
    lpsi = pars(1, j);
    txi = pars(2, j);
    
    m = yc(j).n_rows; // number of obs
    yv = yc(j);
    
    for (int k=0; k < m; k++) {
      
      y = yv.row(k);
      z = 1 + xi * (y - mu) / exp(lpsi);
      
      for (int l=0; l < r; l++) {
        
        if (l == r - 1) {
          
          ee2 = exp(-txi);
          ee3 = 1 + ee2;
          ee5 = 1.5/ee3 - 1;
          ee6 = exp(lpsi);
          ee7 = 1/ee5;
          ee8 = y(l) - mu;
          ee9 = y(l + drop) - mu;
          ee11 = ee5 * ee8/ee6;
          ee13 = ee5 * ee9/ee6;
          ee14 = ee11 + 1;
          ee15 = ee13 + 1;
          ee16 = -ee7;
          ee17 = 1 + ee7;
          ee20 = exp(-R_pow(ee14, ee16));
          ee21 = exp(-R_pow(ee15, ee16));
          ee22 = R_pow(ee14, ee17);
          ee23 = R_pow(ee15, ee17);
          ee24 = ee20 - ee21;
          ee25 = R_pow(ee3, 2);
          ee26 = log1p(ee11);
          ee27 = log1p(ee13);
          ee28 = R_pow(ee14, ee7);
          ee29 = R_pow(ee15, ee7);
          ee34 = 1.5 * (ee8/(ee22 * ee6));
          ee35 = 1.5 * (ee9/(ee23 * ee6));
          ee36 = ee7 + 2;
          ee37 = ee17 * ee5;
          ee38 = ee24 * ee6;
          ee39 = 2 * ee17;
          ee43 = 1.5 * (ee26/(ee28 * ee5)) - ee34;
          ee45 = 1.5 * (ee27/(ee29 * ee5)) - ee35;
          ee46 = R_pow(ee14, ee36);
          ee47 = R_pow(ee15, ee36);
          ee48 = ee3 * ee5;
          ee50 = ee43 * ee20 - ee45 * ee21;
          ee56 = ee21 * ee9/ee23 - ee20 * ee8/ee22;
          ee58 = ee21/ee23 - ee20/ee22;
          ee63 = ee17 * ee8;
          ee64 = ee17 * ee9;
          ee66 = ee25 * ee5 * ee24;
          ee67 = 1.5/ee28;
          ee68 = 1.5/ee29;
          ee70 = 1/R_pow(ee14, ee39) - ee37/ee46;
          ee72 = 1/R_pow(ee15, ee39) - ee37/ee47;
          ee73 = (((ee67 - 1.5) * ee26/ee5 - ee34)/ee5 + 1.5 * (ee63/(ee14 *  ee6))) * ee20;
          ee74 = (((ee68 - 1.5) * ee27/ee5 - ee35)/ee5 + 1.5 * (ee64/(ee15 *  ee6))) * ee21;
          ee75 = (ee70 * ee8/ee6 + 1/ee22) * ee20;
          ee76 = (ee72 * ee9/ee6 + 1/ee23) * ee21;
          ee78 = ee25 * ee24 * ee6;
          ee79 = ee5 * ee24;
          ee80 = R_pow(ee5, 2);
          ee83 = (2.25/ee48 - 3) * ee2/ee3 + 1.5;
          ee86 = (4.5/ee48 - 3) * ee2/ee3 + 1.5;
          
          out(j, 0) += -(ee58/ee38);
          out(j, 1) += -(ee56/ee38);
          out(j, 2) += ee50 * ee2/ee66;
          out(j, 3) += -((ee70 * ee20 - (ee72 * ee21 + R_pow(ee58, 2)/
            ee24))/(ee24 * R_pow(ee6, 2)));
          out(j, 4) += -((ee75 - (ee76 + ee56 * ee58/ee38))/ee38);
          out(j, 5) += -((ee73/ee22 + ee50 * ee58/ee79 - ee74/ee23) * ee2/
            ee78);
          out(j, 6) += -((ee75 * ee8 - (ee76 * ee9 + R_pow(ee56, 2)/ee38))/
            ee38);
          out(j, 7) += -((ee73 * ee8/ee22 + ee50 * ee56/ee79 - ee74 * ee9/
            ee23) * ee2/ee78);
          out(j, 8) += (((((1.5 - ee67) * ee26/ee5 + ee34) * ee43 * ee2/
            ee25 + (2.25 * (ee2 * ee8/(ee14 * ee25 * ee6)) - ee86 * ee26)/
              ee28)/ee5 + (ee83/ee22 - 1.5 * ((1.5 * (ee26/(ee22 * ee80)) -
                1.5 * (ee63/(ee46 * ee6))) * ee2/ee25)) * ee8/ee6) * ee20 +
                R_pow(ee50, 2) * ee2/ee66 - ((((1.5 - ee68) * ee27/ee5 +
                ee35) * ee45 * ee2/ee25 + (2.25 * (ee2 * ee9/(ee15 * ee25 *
                ee6)) - ee86 * ee27)/ee29)/ee5 + (ee83/ee23 - 1.5 * ((1.5 *
                (ee27/(ee23 * ee80)) - 1.5 * (ee64/(ee47 * ee6))) * ee2/ee25)) *
                ee9/ee6) * ee21) * ee2/ee66;
          
          
        } else {
          
          if (l < drop && z(l) <= 0.0) {
            
            ee2 = exp(-txi);
            ee3 = 1 + ee2;
            ee4 = exp(lpsi);
            ee6 = 1.5/ee3 - 1;
            ee7 = y(l + drop) - mu;
            ee8 = ee6 * ee7;
            ee9 = ee8/ee4;
            ee10 = ee9 + 1;
            ee11 = ee10 * ee4;
            ee12 = R_pow(ee3, 2);
            ee13 = ee8/ee11;
            ee14 = ee10 * ee12;
            ee15 = ee7/ee11;
            ee17 = ee14 * ee6 * ee4;
            ee18 = ee12 * ee6;
            ee20 = 1.5 * ee13 * ee2;
            ee21 = log1p(ee9);
            
            out(j, 0) += -(1/ee11);
            out(j, 1) += -ee15;
            out(j, 2) += (1.5 * ee15 - 1.5 * (ee21/ee6)) * ee2/ee18;
            out(j, 3) += -(ee6/(R_pow(ee10, 2) * R_pow(ee4, 2)));
            out(j, 4) += (1 - ee13)/ee11;
            out(j, 5) += ee20/ee17;
            out(j, 6) += -((ee13 - 1) * ee7/ee11);
            out(j, 7) += ee20 * ee7/ee17;
            out(j, 8) += -(((((2.25 * ee15 + 2.25/ee6)/ee3 - 3) * ee2/ee3 +
              1.5) * ee7/ee11 + (2.25 * (ee2 * ee7/(ee14 * ee4)) - ((4.5/
                (ee3 * ee6) - 3) * ee2/ee3 + 1.5) * ee21)/ee6) * ee2/ee18);
            
          } else {
            
            ee2 = exp(-txi);
            ee3 = 1 + ee2;
            ee5 = 1.5/ee3 - 1;
            ee6 = exp(lpsi);
            ee7 = 1/ee5;
            ee8 = y(l) - mu;
            ee9 = y(l + drop) - mu;
            ee11 = ee5 * ee8/ee6;
            ee13 = ee5 * ee9/ee6;
            ee14 = ee11 + 1;
            ee15 = ee13 + 1;
            ee16 = 1 + ee7;
            ee17 = R_pow(ee14, ee7);
            ee18 = R_pow(ee15, ee7);
            ee19 = R_pow(ee14, ee16);
            ee20 = R_pow(ee15, ee16);
            ee23 = 1/ee18 - 1/ee17;
            ee24 = ee7 + 2;
            ee25 = log1p(ee11);
            ee26 = log1p(ee13);
            ee27 = R_pow(ee3, 2);
            ee28 = ee8/ee19;
            ee29 = ee9/ee20;
            ee30 = R_pow(ee14, ee24);
            ee31 = R_pow(ee15, ee24);
            ee32 = 1/ee19;
            ee33 = 1/ee20;
            ee34 = ee23 * ee6;
            ee35 = ee3 * ee5;
            ee37 = (1.5 * ee28 - 1.5 * ee29)/ee6 + (1.5 * (ee26/ee18) -  1.5 * (ee25/ee17))/ee5;
            ee38 = R_pow(ee5, 2);
            ee39 = ee28 - ee29;
            ee40 = ee32 - ee33;
            ee41 = ee30 * ee6;
            ee42 = ee31 * ee6;
            ee43 = ee16 * ee5;
            ee51 = ee27 * ee5 * ee23;
            ee53 = ee27 * ee23 * ee6;
            ee56 = (2.25/ee35 - 3) * ee2/ee3 + 1.5;
            ee59 = (4.5/ee35 - 3) * ee2/ee3 + 1.5;
            ee60 = ee8/ee30;
            ee61 = ee9/ee31;
            ee65 = 1.5 * (ee25/(ee19 * ee38)) - 1.5 * (ee16 * ee8/ee41);
            ee67 = 1.5 * (ee26/(ee20 * ee38)) - 1.5 * (ee16 * ee9/ee42);
            
            out(j, 0) += ee40/ee34;
            out(j, 1) += ee39/ee34;
            out(j, 2) += -(ee37 * ee2/ee51);
            out(j, 3) += -((ee43 * (1/ee31 - 1/ee30) - R_pow(ee40, 2)/ee23)/
              (ee23 * R_pow(ee6, 2)));
            out(j, 4) += -((((ee61 - ee60) * ee16 * ee5 - ee39 * ee40/ee23)/
              ee6 + ee32 - ee33)/ee34);
            out(j, 5) += (((1.5 * (ee25/ee19) - 1.5 * (ee26/ee20))/ee5 -
              ee37 * ee40/ee23)/ee5 + ee16 * (1.5 * ee61 - 1.5 * ee60)/ee6) *
              ee2/ee53;
            out(j, 6) += (R_pow(ee39, 2)/ee34 + (ee33 - ee43 * ee9/ee42) *
              ee9 - (ee32 - ee43 * ee8/ee41) * ee8)/ee34;
            out(j, 7) += (ee65 * ee8 - (ee37 * ee39/(ee5 * ee23) + ee67 *
              ee9)) * ee2/ee53;
            out(j, 8) += -((((ee56/ee20 - 1.5 * (ee67 * ee2/ee27)) * ee9 -
              (ee56/ee19 - 1.5 * (ee65 * ee2/ee27)) * ee8)/ee6 + (((1.5 *
              ((1.5 * (ee26/(ee18 * ee5)) - 1.5 * (ee9/(ee20 * ee6))) *
              ee26) - 1.5 * ((1.5 * (ee25/(ee17 * ee5)) - 1.5 * (ee8/(ee19 *
              ee6))) * ee25))/ee5 - R_pow(ee37, 2)/ee23) * ee2/ee27 + (2.25 *
              (ee2 * ee9/(ee15 * ee27 * ee6)) - ee59 * ee26)/ee18 -
              (2.25 * (ee2 * ee8/(ee14 * ee27 * ee6)) - ee59 * ee25)/ee17)/
                ee5) * ee2/ee51);
            
          }  
          
        }
      
      }
      
    }
    
  }
  
  return out;
  
}
