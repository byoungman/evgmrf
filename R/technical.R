#' Control parameters for evgmrf model optimization
#'
#' Define convergence thresholds, finite difference structures, and regularisation 
#' adjustments for the internal numerical solvers used by \code{\link{evgmrf}}.
#'
#' @param eps A numeric value specifying the convergence tolerance for the 
#'   inner optimization loop. Defaults to \code{5e-3}.
#' @param it0 An integer setting the maximum number of initial iterations 
#'   allowed for the primary optimization routine. Defaults to \code{20}.
#' @param step_size A numeric value defining the step size scaling factor 
#'   applied during optimization line searches. Defaults to \code{0.2}.
#' @param reml_eps A numeric value specifying the convergence tolerance for 
#'   the Restricted Maximum Likelihood (REML) outer loop. Defaults to \code{5e-3}.
#' @param reml_direction A character string or vector detailing the numerical 
#'   differentiation scheme used to evaluate REML parameter gradients. Supported 
#'   options are \code{"forward"}, \code{"backward"}, \code{"central"}, or 
#'   \code{"ad-hoc"}. If \code{"ad-hoc"} is selected, a backward difference is chosen 
#'   for parameters greater than 1, and a forward difference is applied otherwise. 
#'   Defaults to \code{"ad-hoc"}.
#' @param reml_steptol A numeric value defining the minimum step length tolerance 
#'   for REML updates. Defaults to \code{5e-3}.
#' @param reml_stepmax A numeric value specifying the maximum permissible step 
#'   length taken during a single REML update iteration. Defaults to \code{1}.
#' @param reml_itlim An integer setting the hard iteration limit for the REML 
#'   optimization process. Defaults to \code{100}.
#' @param inner_optim A character string selecting the linear algebraic solver 
#'   for inner loop parameters; defaults to \code{"chol"} (Cholesky). If explicitly 
#'   set to \code{"Cholesky"}, \code{super} is automatically forced to \code{TRUE}.
#' @param alpha.tol A numeric convergence tolerance threshold for the model's 
#'   spatial \code{alpha} parameters. Defaults to \code{1e-6}.
#' @param grad_mult A numeric penalty coefficient applied during outer REML optimization. 
#'   If any unpenalised inner optimization gradient coordinate exceeds an absolute 
#'   threshold of 1, this value is added directly as a discrete penalty to the objective 
#'   function score. Defaults to \code{0}.
#' @param par_mult A numeric scaling factor controlling an L2 ridge regularisation penalty 
#'   that shrinks outer REML parameters toward their target baseline vectors. 
#'   Defaults to \code{0.1}.
#' @param update Logical; indicates whether spatial weight structures should be 
#'   iteratively updated during convergence passes. Defaults to \code{FALSE}.
#' @param openmp Logical; toggles whether matrix solving routines exploit multi-core 
#'   CPU parallelism via OpenMP extensions. Defaults to \code{FALSE}.
#' @param threads An integer specifying the target system CPU core thread allocation limit 
#'   used if \code{openmp = TRUE}. A value of \code{0} defaults to system automatic selection.
#' @param perturb.tol A numeric baseline ridge penalty added directly to the diagonal of 
#'   near-singular Hessian matrices to guarantee positive-definiteness during standard 
#'   Cholesky factorization loops. Defaults to \code{1e-2}.
#' @param perturb.mult A scaling multiplier used to exponentially step up the diagonal 
#'   \code{perturb.tol} ridge penalty during iterative Cholesky factorization failures. 
#'   Defaults to \code{5}.
#' @param perturb.method A character string defining the algebraic framework used to apply 
#'   perturbations; defaults to \code{"chol"}. If set to \code{"Cholesky"}, \code{super} 
#'   is forced to \code{TRUE}.
#' @param perturb.tol.eigen A numeric baseline tolerance threshold checking the smallest 
#'   eigenvalue during dedicated Eigen-based matrix perturbation loops. Defaults to \code{1e-3}.
#' @param perturb.mult.eigen A scaling multiplier used to exponentially step up the 
#'   eigenvalue stabilization penalty threshold upon iterative Cholesky failure. 
#'   Defaults to \code{10}.
#' @param cv_eps A numeric convergence tolerance value for cross-validation evaluation loops. 
#'   Defaults to \code{5e-2}.
#' @param cv_gradtol A numeric gradient tolerance threshold used to evaluate cross-validation stability. 
#'   Defaults to \code{0.05}.
#' @param cv_steptol A numeric step-length boundary tolerance used during cross-validation loops. 
#'   Defaults to \code{1e-4}.
#' @param super Logical; activates CHOLMOD supernodal sparse matrix factorization settings. 
#'   Automatically forced to \code{TRUE} if \code{inner_optim} or \code{perturb.method} 
#'   equal \code{"Cholesky"}. Defaults to \code{FALSE}.
#' @param sandwich Logical; if \code{TRUE}, evaluates Huber-White sandwich covariance estimators 
#'   to protect standard error calculations against structural model misspecification. Defaults to \code{FALSE}.
#' @param refine.inits Logical; if \code{TRUE}, initial values for parameters are refined 
#'   using a diagonal quasi-Newton method. Defaults to \code{FALSE}.
#' 
#' @details 
#' Defining explicit arguments simplifies code verification and surfaces automatic 
#' variable completions during interactive package sessions.
#'
#' @return A structured, named \code{list} encompassing all validated logical indicators, 
#'   character frameworks, and numeric convergence boundaries.
#' 
#' @seealso \code{\link{evgmrf}}
#'
#' @examples
#' # Generate default tuning structures
#' config <- evgmrf.control()
#' 
#' # Configure custom parallel threads and adjust the parameter shrinkage weight
#' custom_config <- evgmrf.control(openmp = TRUE, threads = 4, par_mult = 0.5)
#' 
#' @export
evgmrf.control <- function(eps = 5e-3, it0 = 20, step_size = 0.2, reml_eps = 5e-3, 
                           reml_direction = 'ad-hoc', reml_steptol = 5e-3, 
                           reml_stepmax = 1, reml_itlim = 1e2, inner_optim = 'chol',
                           alpha.tol = 1e-6, grad_mult = 0, par_mult = .1, 
                           update = FALSE, openmp = FALSE, threads = 0, 
                           perturb.tol = 1e-2, perturb.mult = 5, perturb.method = 'chol', 
                           perturb.tol.eigen = 1e-3, perturb.mult.eigen = 10,
                           cv_eps = 5e-2, cv_gradtol = .05, cv_steptol = 1e-4, 
                           super = FALSE, sandwich = FALSE, refine.inits = FALSE) {
  
  # 1. Gather all explicit configurations into a named list
  out <- list(
    eps = eps, it0 = it0, step_size = step_size, reml_eps = reml_eps, 
    reml_direction = reml_direction, reml_steptol = reml_steptol, 
    reml_stepmax = reml_stepmax, reml_itlim = reml_itlim, inner_optim = inner_optim,
    alpha.tol = alpha.tol, grad_mult = grad_mult, par_mult = par_mult, 
    update = update, openmp = openmp, threads = threads, 
    perturb.tol = perturb.tol, perturb.mult = perturb.mult, perturb.method = perturb.method, 
    perturb.tol.eigen = perturb.tol.eigen, perturb.mult.eigen = perturb.mult.eigen,
    cv_eps = cv_eps, cv_gradtol = cv_gradtol, cv_steptol = cv_steptol, 
    super = super, sandwich = sandwich, refine.inits = refine.inits
  )
  
  # 2. Apply conditional overrides based on assigned options
  if (out$inner_optim == 'Cholesky' || out$perturb.method == 'Cholesky') {
    out$super <- TRUE
  }
  
  out
}
