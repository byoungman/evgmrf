.model2n <- function(model) {
  if (is.na(model))
    return(0)
  if (model == 'icar') {
    return(1)
  } else {
    return(2)
  }
}

.make_D1 <- function(n) {
  i <- 1:(n - 1)
  j <- c(1:(n - 1), 2:n)
  x <- c(rep(-1, n - 1), rep(1, n - 1))
  Matrix::sparseMatrix(i = rep(i, 2), j = j, x = x, dims = c(n - 1, n))
}

.make_D2 <- function(n) {
  i <- c(1:(n-2), 1:(n-2), 1:(n-2))
  j <- c(1:(n-2), 2:(n-1), 3:n)
  x <- c(rep(1, n-2), rep(-2, n-2), rep(1, n-2))
  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n-2, n))
}

.make_D3 <- function(n) {
  i <- rep(1:(n - 3), each = 4)
  j <- c(1:(n - 3), 2:(n - 2), 3:(n - 1), 4:n)
  x <- rep(c(1, -3, 3, -1), n - 3)
  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n - 3, n))
}

.make_Q <- function(nx, ny, order) {
  if (order == 1) {
    Dx <- .make_D1(nx)
    Dy <- .make_D1(ny)
  }
  if (order == 2) {
    Dx <- .make_D2(nx)
    Dy <- .make_D2(ny)
  }
  if (order == 3) {
    Dx <- .make_D3(nx)
    Dy <- .make_D3(ny)
  }
  Qx <- Matrix::kronecker(Matrix::Diagonal(ny), Matrix::crossprod(Dx))
  Qy <- Matrix::kronecker(Matrix::crossprod(Dy), Matrix::Diagonal(nx))
  as(Qx + Qy, "generalMatrix")
}

.makeQ_data <- function(nx, ny, model, order, n_null, W = NULL) {
  if (is.null(W)) {
    n <- nx * ny
  } else {
    n <- nrow(W)
  }
  # d <- Matrix::colSums(W)
  # ldetD <- sum(log(d))
  # D <- Matrix::Diagonal(n, d)
  # B <- Matrix::Diagonal(n, 1 / Matrix::colSums(W)) %*% W
  id <- seq_along(n_null)
  mods <- ords <- rep(NA, length(id))
  mods[id] <- model
  ords[id] <- order
  rho_spl <- rep(id, sapply(mods, .model2n))
  target <- unlist(list(icar = 1, car = c(1, -4), bym = c(1, -4))[mods])
  Ql <- list()
  for (i in 1:3) {
    if (i %in% order)
      Ql[[i]] <- .make_Q(nx, ny, 1)
  }
  # list(D = D, B = B, W = W, n = n, ldetD = ldetD, n_null = n_null, 
  #      ord = ords, mod = mods, spl = rho_spl, np = length(id), 
  #      target = target, Ql = Ql, nx = nx, ny = ny)
  list(n = n, n_null = n_null, 
       ord = ords, mod = mods, spl = rho_spl, np = length(id), 
       target = target, Ql = Ql, nx = nx, ny = ny)
}

.makeQ_any <- function(pars, Qd, model, order, n_null, R, alpha.tol = 1e-6) {
  if (!is.na(model)) {
    lambda <- exp(pars[1])
    if (model == 'car') {
      alpha.tol <- alpha.tol + (1 - alpha.tol) * pnorm(pars[2])
    } 
    alpha <- 1 - alpha.tol  
    Q <- alpha * Qd$Ql[[order]]
    if (model != 'bym')
      Q <- Q + Matrix::Diagonal(n = Qd$n, x = 1 - alpha)
    if (model == 'icar') {
      logdet <- (Qd$n - order) * pars[1]
    } else {
      if (model == 'car') {
        logdet <- .ldchol(Q)#Matrix::determinant(Q)$modulus
        logdet <- logdet + Qd$n * pars[1]
      }
    }
    Q <- lambda * Q
    if (model == 'bym') {
      Q <- Q + Matrix::Diagonal(n = Qd$n, x = alpha.tol + exp(pars[2]))
      logdet <- .ldchol(Q)
    }
  } else {
    Q <- Matrix::Diagonal(0)
  }
  if (n_null > 0) {
    Q <- Matrix::.bdiag(list(Q, Matrix::Diagonal(n_null, numeric(n_null))))
    Q <- as(Q, 'CsparseMatrix')
  }
  attr(Q, 'logdet') <- logdet
  Q
}

.mQ <- function(rho, Qd, alpha.tol = 1e-6) {
  mods <- Qd$mod
  ords <- Qd$ord
  rho_spl <- .split2(rho, Qd$spl)
  Ql <- list()
  for (i in 1:Qd$np) {
    Ql[[i]] <- .makeQ_any(rho_spl[[i]], Qd, mods[i], ords[i], Qd$n_null[i], Qd$R[[i]], alpha.tol)
  }
  Q <- Matrix::.bdiag(Ql)
  attr(Q, 'logdet') <- sum(sapply(Ql, attr, 'logdet'))
  Q
}

# 
# .orderQ <- function(Q, ord) {
#   Q0 <- Q
#   if (ord > 1) {
#     for (i in 2:ord)
#       Q <- Q %*% Q0
#   }
#   Q
# }
# 
# .orderW <- function(W, ord) {
#   out <- list()
#   out[[1]] <- W
#   if (ord > 1) {
#     for (i in 2:ord) {
#       out[[i]] <- out[[i - 1]] %*% W
#     }
#   }
#   out
# }
# 
# .Wlist2mat <- function(lst, alpha) {
#   ord <- length(lst)
#   if (ord == 1) {
#     return(alpha * lst[[1]])
#   } else {
#     if (ord == 2) {
#       return(2 * alpha * lst[[1]] - (alpha^2) * lst[[2]])
#     } else {
#       if (ord == 3) {
#         return(3 * alpha * lst[[1]] - 3 * (alpha^2) * lst[[2]] + (alpha^3) * lst[[3]])
#       } else {
#         stop("Can't have neighbourhood order greater than 3")
#       }
#     }
#   }
# }
# 
