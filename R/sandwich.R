.J_Q <- function(pars, likdata, likfns, Q, mult = 1) {
  pl <- split(pars, likdata$psplit)
  pm <- t(sapply(seq_along(pl), function(i) as.vector(likdata$Xl[[i]] %*% pl[[i]])))
  likfns$J(pm, likdata) + mult * Q
}

.G_Q <- function(pars, likdata, likfns, Q, mult = 1) {
  pl <- split(pars, likdata$psplit)
  pm <- t(sapply(seq_along(pl), function(i) as.vector(likdata$Xl[[i]] %*% pl[[i]])))
  G <- apply(likfns$G(pm, likdata), 2, c) 
  G <- G + mult * as.vector(Q %*% as.vector(t(pm))) / ncol(G)
  G
}

.G_Q0 <- function(pars, likdata, likfns) {
  pl <- split(pars, likdata$psplit)
  pm <- t(sapply(seq_along(pl), function(i) as.vector(likdata$Xl[[i]] %*% pl[[i]])))
  apply(likfns$G(pm, likdata), 2, c) 
}

sandsim <- function(object, nsim = 1, seed = NULL, type = 'link', prob = NULL,
                            simplify2array = TRUE, decompose = FALSE, ...) {
  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1) # initialize the RNG if necessary
  if(is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  type0 <- type
  if (!is.null(prob))
    type <- 'quantile'
  if (type == 'quantile')
    type0 <- 'response'
  iH <- solve(object$Hessian)
  V <- iH %*% object$J %*% iH
  R <- suppressWarnings(chol(V, pivot = TRUE))
  piv <- order(attr(R, "pivot"))  ## reverse pivoting index
  r <- attr(R, "rank")  ## numerical rank
  R2 <- R[1:r, piv]
  z <- matrix(rnorm(nsim * r), ncol = nsim)
  mat <- object$beta + crossprod(R2, z)
  A <- solve(object$Hessian, object$Gc)
  mat <- object$beta + A %*% matrix(rnorm(ncol(A) * nrow(A)), ncol(A))
  lst <- list()
  for (i in 1:object$np) {
    lst[[i]] <- mat[attr(object$beta, 'split') == i, , drop = FALSE]
    if (decompose &  object$model[i] %in% paste('bym', 2:4, sep = '')) {
      Xlc <- object$likdata$Xlc[[i]]
      Xlc <- Xlc[sapply(Xlc, ncol) > 0]
      ind <- 1:length(Xlc)
      spl <- rep(ind, sapply(Xlc, ncol))
      pl <- lapply(ind, function(j) lst[[i]][spl == j, , drop = FALSE])
      lst[[i]] <- lapply(ind, function(j) Xlc[[j]] %*% pl[[j]])
    } else {
      lst[[i]] <- object$X[[i]] %*% lst[[i]]
    }
  }
  if (decompose) {
    lst <- unlist(lst, recursive = FALSE)
    return(lst)
  }
  if (type %in% c('response', 'quantile')) {
    for (i in 1:object$np) {
      lst[[i]] <- object$unlink[[i]](lst[[i]])
    }
  }
  if (simplify2array) {
    for (i in 1:object$np)
      lst[[i]] <- drop(array(lst[[i]], c(object$nx, object$ny, nsim)))
  }
  names(lst) <- unlist(object$names[type0])
  if (type == 'quantile') {
    nprob <- length(prob)
    if (nprob == 1) {
      lst$p <- prob[1]
      lst <- do.call(object$quantile, lst)
    } else {
      temp <- list()
      for (i in 1:nprob) {
        lst$p <- prob[i]
        temp[[i]] <- do.call(object$quantile, lst)
      }
      lst <- temp
    }
    names(lst) <- paste('q', prob, sep = '_')
  }
  lst
}

sandsim2 <- function(object, nsim = 1, seed = NULL, type = 'link', prob = NULL,
                    simplify2array = TRUE, decompose = FALSE, ...) {
  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1) # initialize the RNG if necessary
  if(is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  type0 <- type
  if (!is.null(prob))
    type <- 'quantile'
  if (type == 'quantile')
    type0 <- 'response'
  iH <- solve(object$Hessian)
  V <- iH# %*% object$J %*% iH
  R <- suppressWarnings(chol(V, pivot = TRUE))
  piv <- order(attr(R, "pivot"))  ## reverse pivoting index
  r <- attr(R, "rank")  ## numerical rank
  R2 <- R[1:r, piv]
  z <- matrix(rnorm(nsim * r), ncol = nsim)
  mat <- object$beta + crossprod(R2, z)
  lst <- list()
  for (i in 1:object$np) {
    lst[[i]] <- mat[attr(object$beta, 'split') == i, , drop = FALSE]
    if (decompose &  object$model[i] %in% paste('bym', 2:4, sep = '')) {
      Xlc <- object$likdata$Xlc[[i]]
      Xlc <- Xlc[sapply(Xlc, ncol) > 0]
      ind <- 1:length(Xlc)
      spl <- rep(ind, sapply(Xlc, ncol))
      pl <- lapply(ind, function(j) lst[[i]][spl == j, , drop = FALSE])
      lst[[i]] <- lapply(ind, function(j) Xlc[[j]] %*% pl[[j]])
    } else {
      lst[[i]] <- object$X[[i]] %*% lst[[i]]
    }
  }
  if (decompose) {
    lst <- unlist(lst, recursive = FALSE)
    return(lst)
  }
  if (type %in% c('response', 'quantile')) {
    for (i in 1:object$np) {
      lst[[i]] <- object$unlink[[i]](lst[[i]])
    }
  }
  if (simplify2array) {
    for (i in 1:object$np)
      lst[[i]] <- drop(array(lst[[i]], c(object$nx, object$ny, nsim)))
  }
  names(lst) <- unlist(object$names[type0])
  if (type == 'quantile') {
    nprob <- length(prob)
    if (nprob == 1) {
      lst$p <- prob[1]
      lst <- do.call(object$quantile, lst)
    } else {
      temp <- list()
      for (i in 1:nprob) {
        lst$p <- prob[i]
        temp[[i]] <- do.call(object$quantile, lst)
      }
      lst <- temp
    }
    names(lst) <- paste('q', prob, sep = '_')
  }
  lst
}
