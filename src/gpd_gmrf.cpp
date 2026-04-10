// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

const double xieps = 0.0001;

// [[Rcpp::export(.tgpdgmrfld0)]]
double tgpdgmrfld0(arma::mat pars, arma::field<arma::vec> yc, arma::field<arma::vec> wc)
{
    
int n = yc.n_rows; // number of locations
int m;

double y, w, lpsi, txi, xi;
double ee1, ee2;
double nllh = 0.0;

arma::vec yv, wv;

for (int j=0; j < n; j++) {
  
  lpsi = pars(0, j);
  txi = pars(1, j);
  xi = 1.5 / (1.0 + exp(-txi)) - 1.0;
  
  m = yc(j).size(); // number of obs
  yv = yc(j);
  wv = wc(j);
  
  for (int k=0; k < m; k++) {

    y = yv(k);
    w = wv(k);
    
    if(std::isfinite(y)) {
      
      ee1 = xi * y / exp(lpsi);
      
      if (ee1 <= -1.0) {
        
        nllh = 1e20;
        break;
        
      } else {
      
        ee2 = 1.0 / xi;
        nllh += w * (lpsi + (ee2 + 1.0) * log1p(ee1));
        
      }
      
    }
    
  } 
  
}


return(nllh);

}

// [[Rcpp::export(.tgpdgmrfld12)]]
arma::mat tgpdgmrfld12(arma::mat pars, arma::field<arma::vec> yc, arma::field<arma::vec> wc)
{
  
  int n = yc.n_rows; // number of locations
  int m;
  
  double y, w, pars1, pars2;
  
  arma::mat out = arma::mat(n, 5, arma::fill::zeros);
  
  double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee16, ee17, ee18, ee19;
  
  arma::vec yv, wv;
  
  for (int j=0; j < n; j++) {
    
    pars1 = pars(0, j);
    pars2 = pars(1, j);
    
    m = yc(j).size(); // number of obs
    yv = yc(j);
    wv = wc(j);
    
    for (int k=0; k < m; k++) {
      
      y = yv(k);
      w = wv(k);
      
      if(std::isfinite(y)) {
      
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
        
        out(j, 0) += w * (1 - ee18/ee10);
        out(j, 1) += w * ((1.5 * (ee14/ee10) - 1.5 * (ee17/ee13)) * ee2/
          ee11);
        out(j, 2) += w * (-(ee18 * (ee19 - 1)/ee10));
        out(j, 3) += w * (-(y * (ee12 * (1.5 - 1.5 * ee19) - 1.5/ee5) * ee2/
          ee16));
        out(j, 4) += w * (-(((2.25 * (y * ee2/ee16) - ((4.5/(ee3 * ee5) -
          3) * ee2/ee3 + 1.5) * ee17)/ee13 + y * (((2.25 * (y/(ee3 *
          ee9 * ee6)) - 3) * ee2/ee3 + 1.5) * ee12 + 2.25 * (ee2/(ee11 *
          ee13)))/ee10) * ee2/ee11));
      
      }
      
    }
    
  }
  
  return out;
  
}

// [[Rcpp::export(.tgpdgmrfldJ)]]
arma::mat tgpdgmrfldJ(arma::mat pars, arma::field<arma::vec> yc, arma::field<arma::vec> wc)
{
  
  int n = yc.n_rows; // number of locations
  int m;
  
  double y, w, pars1, pars2;
  
  arma::mat h = arma::mat(n, 3, arma::fill::zeros);
  arma::cube g = arma::cube(n, m, 2, arma::fill::zeros);

  
  double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee16, ee17, ee18, ee19;
  
  arma::vec yv, wv;
  
  for (int j=0; j < n; j++) {
    
    pars1 = pars(0, j);
    pars2 = pars(1, j);
    
    m = yc(j).size(); // number of obs
    yv = yc(j);
    wv = wc(j);
    
    for (int k=0; k < m; k++) {
      
      y = yv(k);
      w = wv(k);
      
      if(std::isfinite(y)) {
      
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
        
        g(j, k, 0) += w * (1 - ee18/ee10);
        g(j, k, 1) += w * ((1.5 * (ee14/ee10) - 1.5 * (ee17/ee13)) * ee2/
          ee11);
        h(j, 0) += w * (-(ee18 * (ee19 - 1)/ee10));
        h(j, 1) += w * (-(y * (ee12 * (1.5 - 1.5 * ee19) - 1.5/ee5) * ee2/
          ee16));
        h(j, 2) += w * (-(((2.25 * (y * ee2/ee16) - ((4.5/(ee3 * ee5) -
          3) * ee2/ee3 + 1.5) * ee17)/ee13 + y * (((2.25 * (y/(ee3 *
          ee9 * ee6)) - 3) * ee2/ee3 + 1.5) * ee12 + 2.25 * (ee2/(ee11 *
          ee13)))/ee10) * ee2/ee11));
      
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

// [[Rcpp::export(.tgpdgmrfld0_omp)]]
double tgpdgmrfld0_omp(arma::mat pars, arma::field<arma::vec> yc, arma::field<arma::vec> wc, int threads = 0)
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
  double y, w, lpsi, txi, xi;
  double ee1, ee2;

  arma::vec yv, wv;
  
#pragma omp for schedule(static) ordered
  for (int j=0; j < n; j++) {
    
#pragma omp ordered
{
    
  if (!invalid) {
      
    lpsi = pars(0, j);
    txi = pars(1, j);
    xi = 1.5 / (1.0 + exp(-txi)) - 1.0;
    
    m = yc(j).size(); // number of obs
    yv = yc(j);
    wv = wc(j);
    
    for (int k=0; k < m; k++) {
      
      y = yv(k);
      w = wv(k);
      
      if(std::isfinite(y)) {
        
        ee1 = xi * y / exp(lpsi);
        
        if (ee1 <= -1.0) {
          
          invalid = true;
          break;
          
        } else {
          
          ee2 = 1.0 / xi;
          local_nllh += w * (lpsi + (ee2 + 1.0) * log1p(ee1));
          
        }
        
      }
      
    }
    
}
  
}}
  
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

// [[Rcpp::export(.tgpdgmrfld0_omp2)]]
double tgpdgmrfld0_omp2(arma::mat pars, arma::field<arma::vec> yc, arma::field<arma::vec> wc, int threads)
{
  
  int n = yc.n_rows; // number of locations
  int m;
  
  double y, w, lpsi, txi, xi;
  double ee1, ee2;
  arma::vec nllh = arma::vec(n, arma::fill::zeros);
  
  arma::vec yv, wv;

  // Set up OpenMP threads
#ifdef _OPENMP
  if (threads == 0) {
    threads = omp_get_max_threads();
  }
  omp_set_num_threads(threads);
#endif

#pragma omp parallel for private(m, y, w, lpsi, txi, xi, ee1, ee2, yv, wv) shared(nllh, pars, yc, wc, n)
  for (int j=0; j < n; j++) {
    
    lpsi = pars(0, j);
    txi = pars(1, j);
    xi = 1.5 / (1.0 + exp(-txi)) - 1.0;
    
    m = yc(j).size(); // number of obs
    yv = yc(j);
    wv = wc(j);
    
    for (int k=0; k < m; k++) {
      
      y = yv(k);
      w = wv(k);
      
      if(std::isfinite(y)) {
        
        ee1 = xi * y / exp(lpsi);
        
        if (ee1 <= -1.0) {
          
          nllh = 1e20;
          break;
          
        } else {
          
          ee2 = 1.0 / xi;
          nllh(j) += w * (lpsi + (ee2 + 1.0) * log1p(ee1));
          
        }
        
      }
      
    } 
    
  }
  
  double out = arma::accu(nllh);
  
  return out;
  
}

// [[Rcpp::export(.tgpdgmrfld12_omp2)]]
arma::mat tgpdgmrfld12_omp2(arma::mat pars, arma::field<arma::vec> yc, arma::field<arma::vec> wc, int threads)
{
  
  int n = yc.n_rows; // number of locations
  int m;
  
  double y, w, pars1, pars2;
  
  arma::mat out = arma::mat(n, 5, arma::fill::zeros);
  
  double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee16, ee17, ee18, ee19;
  
  arma::vec yv, wv;
  
  // Set up OpenMP threads
#ifdef _OPENMP
  if (threads == 0) {
    threads = omp_get_max_threads();
  }
  omp_set_num_threads(threads);
#endif
  
#pragma omp parallel for private(m, y, w, pars1, pars2, ee2, ee3, yv, wv) shared(pars, yc, wc, n)
  for (int j=0; j < n; j++) {
    
    pars1 = pars(0, j);
    pars2 = pars(1, j);
    
    m = yc(j).size(); // number of obs
    yv = yc(j);
    wv = wc(j);
    
    for (int k=0; k < m; k++) {
      
      y = yv(k);
      w = wv(k);
      
      if(std::isfinite(y)) {
        
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
        
        out(j, 0) += w * (1 - ee18/ee10);
        out(j, 1) += w * ((1.5 * (ee14/ee10) - 1.5 * (ee17/ee13)) * ee2/
          ee11);
        out(j, 2) += w * (-(ee18 * (ee19 - 1)/ee10));
        out(j, 3) += w * (-(y * (ee12 * (1.5 - 1.5 * ee19) - 1.5/ee5) * ee2/
          ee16));
        out(j, 4) += w * (-(((2.25 * (y * ee2/ee16) - ((4.5/(ee3 * ee5) -
          3) * ee2/ee3 + 1.5) * ee17)/ee13 + y * (((2.25 * (y/(ee3 *
          ee9 * ee6)) - 3) * ee2/ee3 + 1.5) * ee12 + 2.25 * (ee2/(ee11 *
          ee13)))/ee10) * ee2/ee11));
        
      }
      
    }
    
  }
  
  return out;
  
}


// [[Rcpp::export(.tgpdgmrfld0_omp3)]]
double tgpdgmrfld0_omp3(arma::mat pars, arma::field<arma::vec> yc, arma::field<arma::vec> wc, int threads)
{
  
  int n = yc.n_rows; // number of locations
  // int m;
  // 
  // double y, w, lpsi, txi, xi;
  // double ee1, ee2;
  arma::vec nllh = arma::vec(n, arma::fill::zeros);
  
  arma::vec yv, wv;
  
  // Set up OpenMP threads
#ifdef _OPENMP
  if (threads == 0) {
    threads = omp_get_max_threads();
  }
  omp_set_num_threads(threads);
#endif
  
#pragma omp parallel shared(nllh, pars, yc, wc, n)
{
  int m;
  double y, w, lpsi, txi, xi, ee1, ee2;
  arma::vec yv, wv;
  
#pragma omp for schedule(dynamic)
  for (int j=0; j < n; j++) {
    
    lpsi = pars(0, j);
    txi = pars(1, j);
    xi = 1.5 / (1.0 + exp(-txi)) - 1.0;
    
    m = yc(j).size(); // number of obs
    yv = yc(j);
    wv = wc(j);
    
    for (int k=0; k < m; k++) {
      
      y = yv(k);
      w = wv(k);
      
      if(std::isfinite(y)) {
        
        ee1 = xi * y / exp(lpsi);
        
        if (ee1 <= -1.0) {
          
          nllh = 1e20;
          break;
          
        } else {
          
          ee2 = 1.0 / xi;
          nllh(j) += w * (lpsi + (ee2 + 1.0) * log1p(ee1));
          
        }
        
      }
      
    } 
    
  }
}
  
  double out = arma::accu(nllh);
  
  return out;
  
}

// [[Rcpp::export(.tgpdgmrfld12_omp)]]
arma::mat tgpdgmrfld12_omp(arma::mat pars, arma::field<arma::vec> yc, arma::field<arma::vec> wc, int threads = 0)
{
  
  int n = yc.n_rows; // number of locations

  // Set up OpenMP threads
#ifdef _OPENMP
  if (threads == 0) {
    threads = omp_get_max_threads();
  }
  omp_set_num_threads(threads);
#endif
  
  arma::mat out = arma::mat(n, 5, arma::fill::zeros);
  
  // Create thread-local matrices that will be combined later
  std::vector<arma::mat> thread_outputs;
  
#ifdef _OPENMP
  thread_outputs.resize(threads, arma::mat(n, 5, arma::fill::zeros));
#else
  thread_outputs.resize(1, arma::mat(n, 5, arma::fill::zeros));
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
  double y, w, pars1, pars2;
  
  arma::mat out = arma::mat(n, 5, arma::fill::zeros);
  
  double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee16, ee17, ee18, ee19;
  
  arma::vec yv, wv;
  
#pragma omp for schedule(static) ordered
  for (int j=0; j < n; j++) {
    
#pragma omp ordered
{
  // Only process if no invalid condition has been found yet
  if (!invalid) {
    
    pars1 = pars(0, j);
    pars2 = pars(1, j);
    
    m = yc(j).size(); // number of obs
    yv = yc(j);
    wv = wc(j);
    
    for (int k=0; k < m; k++) {
      
      y = yv(k);
      w = wv(k);
      
      if(std::isfinite(y)) {
        
        ee2 = exp(-pars2);
        ee3 = 1 + ee2;
        ee5 = 1.5/ee3 - 1;
        ee6 = exp(pars1);
        ee7 = y * ee5;
        ee8 = ee7/ee6;
        
        if (ee8 <= -1.0) {
          invalid = true;
          break; // Can break here since we're in ordered section
        } else {
        
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
        
        local_out(j, 0) += w * (1 - ee18/ee10);
        local_out(j, 1) += w * ((1.5 * (ee14/ee10) - 1.5 * (ee17/ee13)) * ee2/
          ee11);
        local_out(j, 2) += w * (-(ee18 * (ee19 - 1)/ee10));
        local_out(j, 3) += w * (-(y * (ee12 * (1.5 - 1.5 * ee19) - 1.5/ee5) * ee2/
          ee16));
        local_out(j, 4) += w * (-(((2.25 * (y * ee2/ee16) - ((4.5/(ee3 * ee5) -
          3) * ee2/ee3 + 1.5) * ee17)/ee13 + y * (((2.25 * (y/(ee3 *
          ee9 * ee6)) - 3) * ee2/ee3 + 1.5) * ee12 + 2.25 * (ee2/(ee11 *
          ee13)))/ee10) * ee2/ee11));
        
      }
      
    }
    
  }
    
  }}}}
    
    // If invalid condition was encountered, return zero matrix (or handle appropriately)
    if (invalid) {
      // You might want to return a special value or throw an error here
      // For now, returning zero matrix to match the early termination behavior
      return arma::mat(n, 5, arma::fill::zeros);
    }
    
    // Combine results from all threads
    for (size_t t = 0; t < thread_outputs.size(); t++) {
      out += thread_outputs[t];
    }
    
  
  return out;
  
}

// // [[Rcpp::export(.tgpdgmrfld0_omp)]]
// double tgpdgmrfld0_omp(arma::mat pars, arma::field<arma::vec> yc, arma::field<arma::vec> wc, int threads = 0)
// {
//   int n = yc.n_rows; // number of locations
//   double nllh = 0.0;
//   
// #ifdef _OPENMP
//   if (threads == 0) {
//     threads = omp_get_max_threads();
//   }
//   omp_set_num_threads(threads);
// #endif
//   
// #pragma omp parallel for reduction(+:nllh)
//   for (int j = 0; j < n; j++) {
//     double lpsi = pars(0, j);
//     double txi  = pars(1, j);
//     double xi   = 1.5 / (1.0 + exp(-txi)) - 1.0;
//     
//     arma::vec yv = yc(j);
//     arma::vec wv = wc(j);
//     int m = yv.size();
//     
//     double ee1, ee2;
//     double local_nllh = 0.0;
//     
//     for (int i = 0; i < m; i++) {
//       double y = yv(i);
//       if (std::isfinite(y)) {
//       double w = wv(i);
//       
//       ee1 = 1.0 + xi * exp(lpsi) * y;
//       if (ee1 <= 0.0) continue;
//       ee2 = pow(ee1, -1.0 / xi);
//       local_nllh += w * (lpsi + (1.0 + 1.0 / xi) * log(ee1) + ee2);
//     }
//     }
//     
//     nllh += local_nllh;
//   }
//   
//   return nllh;
// }
// 
// // [[Rcpp::export(.tgpdgmrfld12_omp)]]
// arma::mat tgpdgmrfld12_omp(arma::mat pars, arma::field<arma::vec> yc, arma::field<arma::vec> wc, int threads = 0)
// {
//   int n = yc.n_rows; // number of locations
//   arma::mat out = arma::zeros(2, n);
//   
// #ifdef _OPENMP
//   if (threads == 0) {
//     threads = omp_get_max_threads();
//   }
//   omp_set_num_threads(threads);
// #endif
//   
// #pragma omp parallel for
//   for (int j = 0; j < n; j++) {
//     double lpsi = pars(0, j);
//     double txi  = pars(1, j);
//     double xi   = 1.5 / (1.0 + exp(-txi)) - 1.0;
//     double dxi  = 1.5 * exp(-txi) / pow(1.0 + exp(-txi), 2.0);
//     
//     arma::vec yv = yc(j);
//     arma::vec wv = wc(j);
//     int m = yv.size();
//     
//     double ee1, ee2;
//     double l1 = 0.0, l2 = 0.0;
//     
//     for (int i = 0; i < m; i++) {
//       double y = yv(i);
//       double w = wv(i);
//       
//       ee1 = 1.0 + xi * exp(lpsi) * y;
//       if (ee1 <= 0.0) continue;
//       ee2 = pow(ee1, -1.0 / xi);
//       
//       double dlpsi = w * (1.0 - (1.0 + 1.0 / xi) * (xi * exp(lpsi) * y / ee1) - exp(lpsi) * y * ee2);
//       double dtxi = w * dxi * (
//         ((1.0 + 1.0 / xi) * (exp(lpsi) * y / ee1) - (1.0 / xi / xi) * log(ee1)) + (exp(lpsi) * y / ee1) * ee2 * log(ee1)
//       );
//       
//       l1 += dlpsi;
//       l2 += dtxi;
//     }
//     
//     out(0, j) = l1;
//     out(1, j) = l2;
//   }
//   
//   return out;
// }
// 
