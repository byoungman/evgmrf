// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

const double xieps = 0.0001;

// [[Rcpp::export(.tgevgmrfld0)]]
double tgevgmrfld0(arma::mat pars, arma::field<arma::vec> yc)
{
    
int n = yc.n_rows; // number of locations
int m;

double y, mu, lpsi, txi, xi;
double ee1, ee2;
double nllh = 0.0;

arma::vec yv;

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
      
        ee2 = 1.0 / xi;
        nllh += lpsi + (ee2 + 1.0) * log1p(ee1) + R_pow(1.0 + ee1, -ee2);
        
      }
      
    }
    
  } 
  
}


return(nllh);

}

// // [[Rcpp::export(.tgevgmrfld1)]]
// arma::vec tgevgmrfld1(const arma::mat& pars, const arma::field<arma::vec>& yc)
// {
// 
// int n = yc.n_rows; // number of locations
// int m;
// 
// double y, mu, lpsi, txi;
// 
// double ee2, ee3, ee4, ee5, ee6, ee7, ee9;
// double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
// double ee20, ee21, ee22, ee23, ee25, ee27, ee28, ee29;
// double ee30, ee31, ee33, ee34;
// 
// arma::uvec id(3);
// id(0) = 0;
// id(1) = n;
// id(2) = 2 * n;
// arma::vec g(n * 3, arma::fill::zeros);
// arma::vec yv;
// 
// for (int i=0; i < n; i++) {
// 
//   mu = pars(0, i);
//   lpsi = pars(1, i);
//   txi = pars(2, i);
//   
//   m = yc(i).size(); // number of obs
//   yv = yc(i);
// 
//   for (int k=0; k < m; k++) {
// 
//     y = yv(k);
//     
//     if(arma::is_finite(y)) {
// 
//     ee2 = exp(-txi);
//     ee3 = 1 + ee2;
//     ee4 = 1.5/ee3;
//     ee5 = ee4 - 1;
//     ee6 = exp(lpsi);
//     ee7 = y - mu;
//     ee9 = ee5 * ee7/ee6;
//     ee10 = ee9 + 1;
//     ee11 = 1/ee5;
//     ee12 = R_pow(ee10, ee11);
//     ee13 = 1 + ee11;
//     ee14 = ee10 * ee6;
//     ee15 = R_pow(ee3, 2);
//     ee16 = log1p(ee9);
//     ee17 = 1/ee12;
//     ee18 = R_pow(ee10, ee13);
//     ee20 = ee3 * ee5;
//     ee21 = 1 + ee17;
//     ee22 = 1.5 * (ee16/(ee12 * ee5));
//     ee23 = 1.5/ee12;
//     ee25 = (((ee23 - 1.5 * ee5) * ee7/ee14 + 1.5) * ee13 - (1.5 +  ee22)/ee5)/ee10 * ee2;
//     ee27 = ((4.5/ee20 - 3) * ee2/ee3 + 1.5) * ee16;
//     ee28 = ee13 * ee5;
//     ee29 = ee13 * ee7;
//     ee30 = ee15 * ee6;
//     ee31 = R_pow(ee5, 2);
//     ee33 = 1.5 * (ee7/(ee18 * ee6));
//     ee34 = ee4 - ee21;
//     
//     g(i) += -((ee28 - ee17)/ee14);
//     g(i + n) += (ee17 - ee28) * ee7/ee14 + 1;
//     g(i + 2 * n) += (((ee23 - 1.5) * ee16/ee5 - ee33)/ee5 + 1.5 * (ee29/ee14)) * ee2/ee15;
// 
//     // g(id[0]) += -((ee28 - ee17)/ee14);
//     // g(id[1]) += (ee17 - ee28) * ee7/ee14 + 1;
//     // g(id[2]) += (((ee23 - 1.5) * ee16/ee5 - ee33)/ee5 + 1.5 * (ee29/ee14)) * ee2/ee15;
//     
//     }
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
// 
// // [[Rcpp::export(.tgevgmrfld2)]]
// arma::sp_mat tgevgmrfld2(const arma::mat& pars, const arma::field<arma::vec>& yc)
// {
// 
// int n = yc.n_rows; // number of locations
// int m;
//   
// arma::sp_mat H(n * 3, n * 3);
// 
// double y, mu, lpsi, txi;
// 
// double ee2, ee3, ee4, ee5, ee6, ee7, ee9;
// double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
// double ee20, ee21, ee22, ee23, ee25, ee27, ee28, ee29;
// double ee30, ee31, ee33, ee34;
// 
// arma::uvec id(3);
// id(0) = 0;
// id(1) = n;
// id(2) = 2 * n;
// 
// arma::vec h(6, arma::fill::zeros);
// arma::vec yv;
// 
// for (int i=0; i < n; i++) {
// 
//   mu = pars(0, i);
//   lpsi = pars(1, i);
//   txi = pars(2, i);
// 
//   h.zeros();
// 
//   m = yc(i).size(); // number of obs
//   yv = yc(i);
//   
//   for (int k=0; k < m; k++) {
//     
//     y = yv(k);
//     
//     if(arma::is_finite(y)) {
//     
//     ee2 = exp(-txi);
//     ee3 = 1 + ee2;
//     ee4 = 1.5/ee3;
//     ee5 = ee4 - 1;
//     ee6 = exp(lpsi);
//     ee7 = y - mu;
//     ee9 = ee5 * ee7/ee6;
//     ee10 = ee9 + 1;
//     ee11 = 1/ee5;
//     ee12 = R_pow(ee10, ee11);
//     ee13 = 1 + ee11;
//     ee14 = ee10 * ee6;
//     ee15 = R_pow(ee3, 2);
//     ee16 = log1p(ee9);
//     ee17 = 1/ee12;
//     ee18 = R_pow(ee10, ee13);
//     ee20 = ee3 * ee5;
//     ee21 = 1 + ee17;
//     ee22 = 1.5 * (ee16/(ee12 * ee5));
//     ee23 = 1.5/ee12;
//     ee25 = (((ee23 - 1.5 * ee5) * ee7/ee14 + 1.5) * ee13 - (1.5 +  ee22)/ee5)/ee10 * ee2;
//     ee27 = ((4.5/ee20 - 3) * ee2/ee3 + 1.5) * ee16;
//     ee28 = ee13 * ee5;
//     ee29 = ee13 * ee7;
//     ee30 = ee15 * ee6;
//     ee31 = R_pow(ee5, 2);
//     ee33 = 1.5 * (ee7/(ee18 * ee6));
//     ee34 = ee4 - ee21;
//     
//     h(0) +=  - (ee13 * ee34 * ee5/(R_pow(ee10, 2) * R_pow(ee6, 2)));
//     h(1) += (((ee21 - ee4) * ee7/ee14 + 1) * ee13 * ee5 - ee17)/ee10/ee6;
//     h(2) += -(ee25/ee30);
//     h(3) += -(((ee34 * ee7/ee14 - 1) * ee13 * ee5 + ee17)/ee10 * ee7/ee6);
//     h(4) += -(ee25 * ee7/ee30);
//     h(5) += (((((2.25/ee20 - 3) * ee2/ee3 + 1.5)/ee18 - 1.5 * ((1.5 * (ee16/(ee18 * ee31)) -
//       1.5 * (ee29/(R_pow(ee10, ee11 +
//       2) * ee6))) * ee2/ee15)) * ee7/ee6 + (ee27 + (1.5 * ((ee22 -
//       ee33) * ee16/ee5) - 2.25 * (ee7/ee14)) * ee2/ee15 + (2.25 * (ee2 * ee7/(ee10 * ee15 * ee6)) -
//       ee27)/ee12)/ee5)/ee5 -
//       (((2.25 * (ee7/(ee10 * ee3 * ee6)) - 3) * ee2/ee3 + 1.5) * ee13 +
//       2.25 * (ee2/(ee15 * ee31))) * ee7/ee14) * ee2/ee15;
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
//   // H(id(0), id(0)) = h(0);
//   // H(id(1), id(0)) = h(1);
//   // H(id(2), id(0)) = h(2);
//   // H(id(0), id(1)) = h(1);
//   // H(id(1), id(1)) = h(3);
//   // H(id(2), id(1)) = h(4);
//   // H(id(0), id(2)) = h(2);
//   // H(id(1), id(2)) = h(4);
//   // H(id(2), id(2)) = h(5);
//   // 
//   // id += 1;
// 
// }
// 
// return H;
// 
// }
// 
// // [[Rcpp::export(.tgevgmrfld2mat)]]
// arma::mat tgevgmrfld2mat(const arma::mat& pars, const arma::field<arma::vec>& yc)
// {
//   
//   int n = yc.n_rows; // number of locations
//   int m;
//   
//   arma::mat H(n, 6, arma::fill::zeros);
//   
//   double y, mu, lpsi, txi;
//   
//   double ee2, ee3, ee4, ee5, ee6, ee7, ee9;
//   double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
//   double ee20, ee21, ee22, ee23, ee25, ee27, ee28, ee29;
//   double ee30, ee31, ee33, ee34;
//   
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
//       ee2 = exp(-txi);
//       ee3 = 1 + ee2;
//       ee4 = 1.5/ee3;
//       ee5 = ee4 - 1;
//       ee6 = exp(lpsi);
//       ee7 = y - mu;
//       ee9 = ee5 * ee7/ee6;
//       ee10 = ee9 + 1;
//       ee11 = 1/ee5;
//       ee12 = R_pow(ee10, ee11);
//       ee13 = 1 + ee11;
//       ee14 = ee10 * ee6;
//       ee15 = R_pow(ee3, 2);
//       ee16 = log1p(ee9);
//       ee17 = 1/ee12;
//       ee18 = R_pow(ee10, ee13);
//       ee20 = ee3 * ee5;
//       ee21 = 1 + ee17;
//       ee22 = 1.5 * (ee16/(ee12 * ee5));
//       ee23 = 1.5/ee12;
//       ee25 = (((ee23 - 1.5 * ee5) * ee7/ee14 + 1.5) * ee13 - (1.5 +  ee22)/ee5)/ee10 * ee2;
//       ee27 = ((4.5/ee20 - 3) * ee2/ee3 + 1.5) * ee16;
//       ee28 = ee13 * ee5;
//       ee29 = ee13 * ee7;
//       ee30 = ee15 * ee6;
//       ee31 = R_pow(ee5, 2);
//       ee33 = 1.5 * (ee7/(ee18 * ee6));
//       ee34 = ee4 - ee21;
//       
//       H(i, 0) +=  - (ee13 * ee34 * ee5/(R_pow(ee10, 2) * R_pow(ee6, 2)));
//       H(i, 1) += (((ee21 - ee4) * ee7/ee14 + 1) * ee13 * ee5 - ee17)/ee10/ee6;
//       H(i, 2) += -(ee25/ee30);
//       H(i, 3) += -(((ee34 * ee7/ee14 - 1) * ee13 * ee5 + ee17)/ee10 * ee7/ee6);
//       H(i, 4) += -(ee25 * ee7/ee30);
//       H(i, 5) += (((((2.25/ee20 - 3) * ee2/ee3 + 1.5)/ee18 - 1.5 * ((1.5 * (ee16/(ee18 * ee31)) -
//         1.5 * (ee29/(R_pow(ee10, ee11 +
//         2) * ee6))) * ee2/ee15)) * ee7/ee6 + (ee27 + (1.5 * ((ee22 -
//         ee33) * ee16/ee5) - 2.25 * (ee7/ee14)) * ee2/ee15 + (2.25 * (ee2 * ee7/(ee10 * ee15 * ee6)) -
//         ee27)/ee12)/ee5)/ee5 -
//         (((2.25 * (ee7/(ee10 * ee3 * ee6)) - 3) * ee2/ee3 + 1.5) * ee13 +
//         2.25 * (ee2/(ee15 * ee31))) * ee7/ee14) * ee2/ee15;    
//       
//     }
//     
//     // H(i, i) = h(0);
//     // H(i + n, i) = h(1);
//     // H(i, i + n) = h(1);
//     // H(i + 2 * n, i) = h(2);
//     // H(i, i + 2 * n) = h(2);
//     // H(i + n, i + n) = h(3);
//     // H(i + n, i + 2 * n) = h(4);
//     // H(i + 2 * n, i + n) = h(4);
//     // H(i + 2 * n, i + 2 * n) = h(5);
//     // 
//     // H(id(0), id(0)) = h(0);
//     // H(id(1), id(0)) = h(1);
//     // H(id(2), id(0)) = h(2);
//     // H(id(0), id(1)) = h(1);
//     // H(id(1), id(1)) = h(3);
//     // H(id(2), id(1)) = h(4);
//     // H(id(0), id(2)) = h(2);
//     // H(id(1), id(2)) = h(4);
//     // H(id(2), id(2)) = h(5);
//     // 
//     // id += 1;
//     
//   }
//   
//   return H;
//   
// }

// [[Rcpp::export(.tgevgmrfld12)]]
arma::mat tgevgmrfld12(arma::mat pars, arma::field<arma::vec> yc)
{
  
  int n = yc.n_rows; // number of locations
  int m;
  
  double y, mu, lpsi, txi, xi;
  
  arma::mat out = arma::mat(n, 9, arma::fill::zeros);
  
  double ee2, ee3, ee4, ee5, ee6, ee7, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
  double ee20, ee21, ee22, ee23, ee25, ee27, ee28, ee29;
  double ee30, ee31, ee33, ee34;
  
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
      
      out(j, 0) += -((ee28 - ee17)/ee14);
      out(j, 1) += (ee17 - ee28) * ee7/ee14 + 1;
      out(j, 2) += (((ee23 - 1.5) * ee16/ee5 - ee33)/ee5 + 1.5 * (ee29/
        ee14)) * ee2/ee15;
      out(j, 3) += -(ee13 * ee34 * ee5/(R_pow(ee10, 2) * R_pow(ee6, 2)));
      out(j, 4) += (((ee21 - ee4) * ee7/ee14 + 1) * ee13 * ee5 - ee17)/
        ee10/ee6;
      out(j, 5) += -(ee25/ee30);
      out(j, 6) += -(((ee34 * ee7/ee14 - 1) * ee13 * ee5 + ee17)/ee10 *
        ee7/ee6);
      out(j, 7) += -(ee25 * ee7/ee30);
      out(j, 8) += (((((2.25/ee20 - 3) * ee2/ee3 + 1.5)/ee18 - 1.5 *
        ((1.5 * (ee16/(ee18 * ee31)) - 1.5 * (ee29/(R_pow(ee10, (ee11 +
        2)) * ee6))) * ee2/ee15)) * ee7/ee6 + (ee27 + (1.5 * ((ee22 -
        ee33) * ee16/ee5) - 2.25 * (ee7/ee14)) * ee2/ee15 +
        (2.25 * (ee2 * ee7/(ee10 * ee15 * ee6)) - ee27)/ee12)/ee5)/
          ee5 - (((2.25 * (ee7/(ee10 * ee3 * ee6)) - 3) * ee2/ee3 + 1.5) *
            ee13 + 2.25 * (ee2/(ee15 * ee31))) * ee7/ee14) * ee2/ee15;
      
      }
      
    }
    
  }
  
  return out;
  
}

// [[Rcpp::export(.tgevgmrfld0_omp)]]
double tgevgmrfld0_omp(arma::mat pars, arma::field<arma::vec> yc, int threads = 0)
{
  
  int n = yc.n_rows; // number of locations
  
  // Set up OpenMP threads
#ifdef _OPENMP
  if (threads == 0) {
    threads = omp_get_max_threads();
  }
  omp_set_num_threads(threads);
#endif
  
  // Vector to store local results from each thread
  std::vector<double> partial_nllh(threads, 0.0);
  double nllh = 0.0;
  bool invalid = false;
  
#pragma omp parallel shared(invalid)
{
  int thread_id = 0;
#ifdef _OPENMP
  thread_id = omp_get_thread_num();
#endif
  
  double local_nllh = 0.0;
  int m;
  double y, mu, lpsi, txi, xi;
  double ee1, ee2;
  arma::vec yv;
  
  // Use static scheduling to ensure deterministic order
#pragma omp for schedule(static) ordered
  for (int i=0; i < n; i++) {
    
#pragma omp ordered
{
  // Only process if no invalid condition has been found yet
  if (!invalid) {
    
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
          invalid = true;
          break; // Can break here since we're in ordered section
        } else {
          
          ee2 = 1.0 / xi;
          local_nllh += lpsi + (ee2 + 1.0) * log1p(ee1) + R_pow(1.0 + ee1, -ee2);
          
        }
        
      }
      
    }
    
  }
}
  }
  
  partial_nllh[thread_id] = local_nllh;
}

// Check if we encountered an invalid situation
if (invalid) {
  return 1e20;
}

// Sum up results from all threads
for (int t = 0; t < threads; t++) {
  nllh += partial_nllh[t];
}

return nllh;
}

// [[Rcpp::export(.tgevgmrfld12_omp)]]
arma::mat tgevgmrfld12_omp(arma::mat pars, arma::field<arma::vec> yc, int threads = 0)
{
  
  int n = yc.n_rows; // number of locations
  
  // Set up OpenMP threads
#ifdef _OPENMP
  if (threads == 0) {
    threads = omp_get_max_threads();
  }
  omp_set_num_threads(threads);
#endif
  
  arma::mat out = arma::mat(n, 9, arma::fill::zeros);
  
  // Create thread-local matrices that will be combined later
  std::vector<arma::mat> thread_outputs;
  
#ifdef _OPENMP
  thread_outputs.resize(threads, arma::mat(n, 9, arma::fill::zeros));
#else
  thread_outputs.resize(1, arma::mat(n, 9, arma::fill::zeros));
#endif
  
  // Add early termination flag for consistency
  bool invalid = false;
  
#pragma omp parallel shared(invalid)
{
  int thread_id = 0;
#ifdef _OPENMP
  thread_id = omp_get_thread_num();
#endif
  
  arma::mat& local_out = thread_outputs[thread_id];
  int m;
  double y, mu, lpsi, txi, xi;
  double ee2, ee3, ee4, ee5, ee6, ee7, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
  double ee20, ee21, ee22, ee23, ee25, ee27, ee28, ee29;
  double ee30, ee31, ee33, ee34;
  arma::vec yv;
  
  // Use static scheduling to ensure deterministic order
#pragma omp for schedule(static) ordered
  for (int j=0; j < n; j++) {
    
#pragma omp ordered
{
  // Only process if no invalid condition has been found yet
  if (!invalid) {
    
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
        ee4 = 1.5/ee3;
        ee5 = ee4 - 1;
        ee6 = exp(lpsi);
        ee7 = y - mu;
        ee9 = ee5 * ee7/ee6;
        ee10 = ee9 + 1;
        
        // Check for invalid condition (same as in tgevgmrfld0_omp)
        if (ee5 * (y - mu) / ee6 <= -1.0) {
          invalid = true;
          break; // Can break here since we're in ordered section
        } else {
          
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
          
          local_out(j, 0) += -((ee28 - ee17)/ee14);
          local_out(j, 1) += (ee17 - ee28) * ee7/ee14 + 1;
          local_out(j, 2) += (((ee23 - 1.5) * ee16/ee5 - ee33)/ee5 + 1.5 * (ee29/
            ee14)) * ee2/ee15;
          local_out(j, 3) += -(ee13 * ee34 * ee5/(R_pow(ee10, 2) * R_pow(ee6, 2)));
          local_out(j, 4) += (((ee21 - ee4) * ee7/ee14 + 1) * ee13 * ee5 - ee17)/
            ee10/ee6;
          local_out(j, 5) += -(ee25/ee30);
          local_out(j, 6) += -(((ee34 * ee7/ee14 - 1) * ee13 * ee5 + ee17)/ee10 *
            ee7/ee6);
          local_out(j, 7) += -(ee25 * ee7/ee30);
          local_out(j, 8) += (((((2.25/ee20 - 3) * ee2/ee3 + 1.5)/ee18 - 1.5 *
            ((1.5 * (ee16/(ee18 * ee31)) - 1.5 * (ee29/(R_pow(ee10, (ee11 +
            2)) * ee6))) * ee2/ee15)) * ee7/ee6 + (ee27 + (1.5 * ((ee22 -
            ee33) * ee16/ee5) - 2.25 * (ee7/ee14)) * ee2/ee15 +
            (2.25 * (ee2 * ee7/(ee10 * ee15 * ee6)) - ee27)/ee12)/ee5)/
              ee5 - (((2.25 * (ee7/(ee10 * ee3 * ee6)) - 3) * ee2/ee3 + 1.5) *
                ee13 + 2.25 * (ee2/(ee15 * ee31))) * ee7/ee14) * ee2/ee15;
          
        }
        
      }
      
    }
    
  }
}
  }
}

// If invalid condition was encountered, return zero matrix (or handle appropriately)
if (invalid) {
  // You might want to return a special value or throw an error here
  // For now, returning zero matrix to match the early termination behavior
  return arma::mat(n, 9, arma::fill::zeros);
}

// Combine results from all threads
for (size_t t = 0; t < thread_outputs.size(); t++) {
  out += thread_outputs[t];
}

return out;

}