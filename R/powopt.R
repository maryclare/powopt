#' Penalized likelihood function, used to differentiate between local modes
powObj <- function(beta, X, y, sigma.sq, lambda, q, Q = NULL, l = NULL) {

  if (is.null(Q)) {Q <- crossprod(X)}
  if (is.null(l)) {l <- crossprod(X, y)}
  log.lik <- -(1/(2*sigma.sq))*(crossprod(t(crossprod(beta, Q)), beta) - 2*crossprod(beta, l)) - lambda*sum(abs(beta)^q)

  return(-1*log.lik)
}
#' Thresholding function for penalized regression with power penalty
#'
#' @name powThresh
#'
#' @description Gives the value of \eqn{\beta} that minimizes: \cr
#' \deqn{(\beta - z)^2/2 + \lambda|\beta|^q} \cr
#' for fixed \eqn{z}, \eqn{\lambda \ge 0} and \eqn{q > 0}. \cr \cr
#' For \eqn{q \le 1}, uses thresholding function given in Marjanovic and Solo (2014).
#' For \eqn{q > 1}, computed using fixed point iteration.
#'
#' @usage \code{powThresh(z =  1, lambda = 1, q = 1)}
#'
#' @param \code{z} vector of values to find the minimizing beta for (\code{powThresh} is applied coordinatewise when \code{length(z)} > 1)
#' @param \code{lambda} scalar multiplier applied to penalty
#' @param \code{q} scalar power of penalty function
#'
#' @source For \eqn{q\le 1}, uses thresholding function given in Marjanovic and Solo (2014).
#' Marjanovic, Goran, and Victor Solo. "\eqn{l_q} Sparsity Penalized Linear Regression With Cyclic Descent." IEEE Transactions on Signal Processing 62.6 (2014): 1464-1475.
#'
#' @export
powThresh <- function(z,
                     lambda,
                     q) {
  if (q <= 0) {
    cat("Values of q less than or equal to zero not supported!\n")
    break
  }
  if (lambda < 0) {
    cat("Values of lambda less zero not supported!\n")
    break
  }
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
#' Coordinate descent for penalized regression with power penalty
#'
#' @name powCD
#'
#' @description Gives the value of of the length \eqn{p} vector \eqn{\beta} that minimizes: \cr
#' \deqn{(y - X\beta)^2/(2\sigma^2) + \lambda ||\beta||^q_q} \cr
#' for fixed \eqn{y}, \eqn{X}, \eqn{\sigma^2 > 0}, \eqn{\lambda \ge 0} and \eqn{q > 0}. \cr \cr
#'
#' This corresponds to finding the posterior mode for \eqn{beta} given \eqn{X}, \eqn{y}, \eqn{\sigma^2}, \eqn{\lambda} and \eqn{q} under the model:
#' \deqn{y = X\beta + \sigma Z}
#' \deqn{p(\beta_i) = q\lambda^(1/q) exp(-\lambda|\beta_i|^q)/(2 \Gamma(1/q)),}
#' where elements of Z are independent, standard normal random variates.
#'
#' This distribution for elements of \eqn{\beta} is a generalized normal distribution
#' for \eqn{\beta} with scale \eqn{\alpha = \lambda^(-1/q)} and shape \eqn{\beta = q} (Box and Tiao, 1973).
#'
#' For \eqn{q \le 1}, uses coordinate descent algorithm given in Marjanovic and Solo (2014), modified to accomodate X that do not have standardized columns.
#'
#' @usage \code{powCD(X = X, y = y, sigma.sq = 1, lambda = 1,
#'       q = 1, max.iter = 10000,
#'       print.iter = FALSE,
#'       tol = 10^(-7), ridge.eps = 10^(-7),
#'       rand.restart = 0)}
#'
#' @param \code{X} design matrix
#' @param \code{y} response vector
#' @param \code{sigma.sq} scalar value \eqn{\sigma^2}
#' @param \code{lambda} scalar value \eqn{\lambda}
#' @param \code{q} scalar value \eqn{q}
#' @param \code{max.iter} maximum number of iterations for coordinate descent
#' @param \code{print.iter} logical value indicating whether iteration count for coordinate descent should be printed
#' @param \code{tol} scalar tolerance value for assessing whether or not gradient of objective function is sufficiently close to zero
#' @param \code{ridge.eps} ridge regression tuning parameter for obtaining starting value of \eqn{\beta} in coordinate descent
#' @param \code{rand.restart} number of times coordinate descent should be repeated from random starting value for \eqn{\beta} after an initial application of coordinate descent starting from ridge solution, needed when \eqn{X} is not orthogonal because the coordinate descent algorithm is not guaranteed to converge to the global optimum for all non-orthogonal \eqn{X}
#'
#' @return Returns a vector of optimal values for \eqn{\beta}. If the coordinate descent algorithm does not meet the optimality conditions given in Marjanovic and Solo (2014), a vector of \code{NA}'s is returned.
#'
#' Note that a non-\code{NA} solution for \eqn{\beta} guarantees the global minimum has been attained when \eqn{X} is full rank and orthogonal. Otherwise, the solution will only correspond to a global minimum when the conditions on \eqn{X} given in Marjanovic and Solo (2014) are satisfied.
#'
#' @source For \eqn{q\le 1}, uses coordinate descent algorithm given by Marjanovic and Solo (2014), modified to accomodate X that do not have standardized columns. \cr
#'
#' Box, G. E. P., and G. C. Tiao. "Bayesian inference in statistical inference." Adison-Wesley, Reading, Mass (1973). \cr
#'
#' Marjanovic, Goran, and Victor Solo. "\eqn{l_q} Sparsity Penalized Linear Regression With Cyclic Descent." IEEE Transactions on Signal Processing 62.6 (2014): 1464-1475. \cr
#'
#' @export
powCD <- function(X, y, sigma.sq, lambda, q, max.iter = 10000,
                   print.iter = FALSE, tol = 10^(-7), ridge.eps = 10^(-7), rand.restart = 0) {


  if (q <= 0) {
    cat("Values of q less than or equal to zero not supported!\n")
    break
  }
  if (lambda < 0) {
    cat("Values of lambda less zero not supported!\n")
    break
  }
  if (sigma.sq <= 0) {
    cat("Values of sigma.sq less zero or equal to zero not supported!\n")
    break
  }

  l <- crossprod(X, y)
  Q <- crossprod(X)

  orthx <- mean(abs(Q[lower.tri(Q, diag = FALSE)])) <= 10^(-14)
  fullx <- min(eigen(Q)$values) > 0
  if (!orthx & rand.restart == 0) {
    cat("The design matrix is not orthogonal. It is possible that the coordinate descent algorithm will not converge to the global minimum regardless of the starting value. Setting rand.restart > 0 and examining the solution is strongly recommended.\n")
  }
  if (orthx & fullx & rand.restart > 0) {
    cat("The design matrix is full rank and orthogonal and coordinate descent will always converge to the global minimum, so there is no need to repeat coordinate descent from random starting values. The value of rand.restart will be reset to zero.\n")
    rand.restart <- 0
  }

  lambda <- sigma.sq*lambda

  p <- ncol(X)

  bb.tmp <- rep(NA, p)
  obj.tmp <- Inf

  iter <- ifelse(rand.restart == 0, 1, rand.restart + 1)

  for (m in 1:iter) {
    if (m == 1) {
      bb <- crossprod(solve(Q + ridge.eps*diag(p)), l)
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
        X.ii <- Q[i, i]
        z.k.i <- (crossprod(X[, i], y - crossprod(t(X[, -i]), bb[-i])))/X.ii
        lambda.ii <- lambda/X.ii
        b.kp.i <- powThresh(z.k.i, lambda = lambda.ii, q = q)
        rr <- rr - (b.kp.i - bb[i])*X[, i]
        bb[i] <- b.kp.i
      }
      k <- k + 1

      opt <- numeric(p)

      for (i in 1:p) {
        X.ii <- Q[i, i]
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
    obj.bb <- powObj(beta = bb, X = X, y = y, sigma.sq = 1, lambda = lambda, q = q, Q = Q, l = l)
    if (obj.bb <= obj.tmp) {
      bb.tmp <- bb
      obj.tmp <- obj.bb
    }
  }
  if (is.infinite(bb.tmp[1])) {bb.tmp <- rep(NA, p)}
  return(bb.tmp)
}
