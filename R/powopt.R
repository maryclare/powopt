#' Penalized likelihood function, used to differentiate between local modes
powObj <- function(beta, X, y, sigma.sq.z, lambda, q) {
  log.lik <- -(1/(2*sigma.sq.z))*(crossprod(crossprod(t(X), beta)) - 2*crossprod(beta, crossprod(X, y))) - lambda*sum(abs(beta)^q)
  return(-1*log.lik)
}
#' Thresholding function for penalized regression with power penalty
#'
#' @name powThresh
#'
#' @description Gives the value of beta that minimizes (beta - z)^2/2 + lambda|beta|^q for fixed z, lambda and q. For q <= 1, uses thresholding function given in Marjanovic and Solo (2014). For q > 1, computed using fixed point iteration.
#'
#' @usage \code{powThresh(z =  1, lambda = 1, q = 1)}
#'
#' @param \code{z} vector of values to find the minimizing beta for (\code{powThresh} is applied coordinatewise when \code{length(z)} > 1)
#' @param \code{lambda} scalar multiplier applied to penalty
#' @param \code{q} power of penalty function
#'
#' @source For q <= 1, uses thresholding function given in Marjanovic and Solo (2014).
#' Marjanovic, Goran, and Victor Solo. "$ l_ {q} $ Sparsity Penalized Linear Regression With Cyclic Descent." IEEE Transactions on Signal Processing 62.6 (2014): 1464-1475.
#'
#' @export
powThresh <- function(z,
                     lambda,
                     q) {

  js <- numeric(length(z))

  if (q <= 1) {

    beta.lambda <- (2*lambda*(1 - q))^(1/(2 - q))
    h.lambda <- beta.lambda + lambda*q*beta.lambda^(q - 1)

    for (j in 1:length(js)) {
      if (abs(z[j]) < h.lambda) {
        js[j] <- 0
      } else if (abs(z[j]) == h.lambda) {
        js[j] <- max(0, beta.lambda)
      } else {

        beta.bar.0 <- Inf
        beta.bar.1 <- abs(z[j])
        while(abs(beta.bar.0 - beta.bar.1) > 10^(-7)) {
          beta.bar.0 <- beta.bar.1
          beta.bar.1 <- abs(z[j]) - lambda*q*beta.bar.0^(q - 1)
        }
        beta.bar <- beta.bar.1

        js[j] <- sign(z[j])*beta.bar
      }
    }
  } else {
    for (j in 1:length(js)) {
      beta.bar.0 <- 10^(-14)
      beta.bar.1 <- beta.bar.0 - f(beta.bar.0, q = q, lambda = lambda, z = z[j])/f.p(beta.bar.0, q = q, lambda = lambda)

      while (abs(beta.bar.1 - beta.bar.0) > 10^(-14)) {
        beta.bar.0 <- beta.bar.1
        beta.bar.1 <- beta.bar.0 - (beta.bar.0 + q*lambda*beta.bar.0^(q - 1) - abs(z[j]))/(1 + (q - 1)*q*lambda*beta.bar.0^(q - 2))
      }
      beta.bar <- beta.bar.1
      js[j] <- sign(z[j])*beta.bar
    }
  }
  return(js)
}

powCD <- function(X, y, sigma.sq.z, lambda, q, max.iter = 10000,
                   print.iter = FALSE, tol = 10^(-7), eps = 10^(-7), ridge.start = TRUE, rand.restart = 0) {

  lambda <- sigma.sq.z*lambda

  p <- ncol(X)

  bb.tmp <- rep(NA, p)
  obj.tmp <- Inf

  iter <- ifelse(ridge.start | rand.restart == 0, 1, rand.restart)

  for (l in 1:iter) {
    if (ridge.start) {
      bb <- crossprod(solve(crossprod(X) + eps*diag(p)), crossprod(X, y))
    } else {
      bb <- rnorm(p)
    }

    rr <- y - crossprod(t(X), bb)
    k <- 1
    opt.cond <- FALSE
    while (k <= max.iter & !opt.cond) {

      if (print.iter) {
        cat("k = ", k, "\n")
      }

      for (i in 1:p) {
        X.ii <- crossprod(X[, i])
        z.k.i <- (crossprod(X[, i], y - crossprod(t(X[, -i]), bb[-i])))/X.ii
        lambda.ii <- lambda/X.ii
        b.kp.i <- powThresh(z.k.i, lambda = lambda.ii, q = q)
        rr <- rr - (b.kp.i - bb[i])*X[, i]
        bb[i] <- b.kp.i
      }
      k <- k + 1

      opt <- numeric(p)

      for (i in 1:p) {
        X.ii <- crossprod(X[, i])
        z.k.i <- (crossprod(X[, i], y - crossprod(t(X[, -i]), bb[-i])))/X.ii
        lambda.ii <- lambda/X.ii
        if (q <= 1) {
          beta.lambda.ii <- (2*lambda.ii*(1 - q))^(1/(2 - q))
          h.lambda.ii <- beta.lambda.ii + lambda.ii*q*beta.lambda.ii^(q - 1)
        }
        if (bb[i] == 0) {
          opt[i] <- abs(bb[i] - z.k.i) <= h.lambda.ii
        } else {
          if (q <= 1) {
            opt[i] <- (abs(bb[i]) >= beta.lambda.ii)*((bb[i] - z.k.i + sign(bb[i])*lambda.ii*q*abs(bb[i])^(q - 1)) <=  tol)
          } else {
            opt[i] <- ((bb[i] - z.k.i + sign(bb[i])*lambda.ii*q*abs(bb[i])^(q - 1)) <=  tol)
          }
        }
      }

      opt.cond <- (sum(opt) == p)

    }
    if (!opt.cond) {bb <- rep(Inf, p)}
    obj.bb <- powObj(beta = bb, X = X, y = y, sigma.sq.z = 1, lambda = lambda, q = q)
    if (obj.bb <= obj.tmp) {
      bb.tmp <- bb
      obj.tmp <- obj.bb
    }
  }
  if (is.infinite(bb.tmp[1])) {bb.tmp <- rep(NA, p)}
  return(bb.tmp)
}
