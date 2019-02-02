#' Coordinate descent for penalized regression with power penalty
#'
#' @name powCD
#'
#' @description Gives the value of of the length \eqn{p} vector \eqn{\beta} that minimizes: \cr
#' \deqn{(y - X\beta)^2/(2\sigma^2) + \lambda ||\beta||^q_q} \cr
#' for fixed \eqn{y}, \eqn{X}, \eqn{\sigma^2 > 0}, \eqn{\lambda > 0 0} and \eqn{q > 0}. \cr \cr
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
#'       tol = 10^(-7), ridge.eps = 0,
#'       rand.restart = 0)}
#'
#' @param \code{X} design matrix
#' @param \code{y} response vector
#' @param \code{sigma.sq} scalar value \eqn{\sigma^2}
#' @param \code{lambda} scalar value \eqn{\lambda}
#' @param \code{q} scalar value \eqn{q}
#' @param \code{max.iter} maximum number of iterations for coordinate descent
#' @param \code{print.iter} logical value indicating whether iteration count for coordinate descent should be printed
#' @param \code{tol} scalar tolerance value for assessing convergence of objective function
#' @param \code{ridge.eps} ridge regression tuning parameter for obtaining starting value of \eqn{\beta} in coordinate descent, defult is zero
#' @param \code{rand.restart} number of times coordinate descent should be repeated from random starting value for \eqn{\beta} after an initial application of coordinate descent starting from ridge solution, needed when \eqn{X} is not orthogonal because the coordinate descent algorithm is not guaranteed to converge to the global optimum for all non-orthogonal \eqn{X} when \eqn{q \le 1}
#' @param \code{start} starting value, set to null by default
#' @param \code{return.obj.iter} TRUE/FALSE value indicating whether or not objectives at convergence and # of coordinate descent iterations should be returned
#' @param \code{order} Vector indicating order of variables in coordinate descent, defaults to 1:p
#'
#' @return Returns a vector of optimal values for \eqn{\beta}. If the coordinate descent algorithm does not meet the optimality conditions given in Marjanovic and Solo (2014), a vector of \code{NA}'s is returned.
#'
#' Note that a non-\code{NA} solution for \eqn{\beta} for any value of \eqn{q} guarantees the global minimum has been attained when \eqn{X} is full rank and orthogonal. Otherwise, when \eqn{q \le 1} the solution will only correspond to a global minimum when the conditions on \eqn{X} given in Marjanovic and Solo (2014) are satisfied.
#'
#' @source For \eqn{q\le 1}, uses coordinate descent algorithm given by Marjanovic and Solo (2014), modified to accomodate X that do not have standardized columns. \cr
#'
#' Box, G. E. P., and G. C. Tiao. "Bayesian inference in statistical inference." Adison-Wesley, Reading, Mass (1973). \cr
#'
#' Marjanovic, Goran, and Victor Solo. "\eqn{l_q} Sparsity Penalized Linear Regression With Cyclic Descent." IEEE Transactions on Signal Processing 62.6 (2014): 1464-1475. \cr
#'
#' @export
powCD <- function(X, y, sigma.sq, lambda, q, max.iter = 10000,
                  print.iter = FALSE, tol = 10^(-7), ridge.eps = 0, rand.restart = 0,
                  start = NULL, return.obj.iter = FALSE, order = 1:p) {


  if (q <= 0) {
    if (print.iter) {cat("Values of q less than or equal to zero not supported!\n")}
    break
  }
  if (lambda < 0) {
    if (print.iter) {cat("Values of lambda less than or equal to zero not supported!\n")}
    break
  }
  if (sigma.sq <= 0) {
    if (print.iter) {cat("Values of sigma.sq less zero or equal to zero not supported!\n")}
    break
  }

  n <- length(y)
  p <- ncol(X)
  l <- crossprod(X, y)
  Q <- crossprod(X)
  yty <- crossprod(y)

  fullx <- min(eigen(Q)$values) > 0
  orthx <- max(abs(Q[lower.tri(Q)])) < 10^(-14)
  if (!fullx & rand.restart == 0 & q == 1) {
    if (print.iter) {cat("The design matrix is not full rank and q = 1. It is possible that the coordinate descent algorithm will not converge to the global minimum regardless of the starting value. Setting rand.restart > 0 and examining the solution is strongly recommended.\n")}
  }
  if (rand.restart == 0 & q < 1) {
    if (print.iter) {cat("The problem is nonconvex. It is possible that the coordinate descent algorithm will not converge to the global minimum regardless of the starting value. Setting rand.restart > 0 and examining the solution is strongly recommended.\n")}
  }
  if (orthx & fullx & (rand.restart > 0 | ridge.eps > 0)) {
    if (print.iter) {cat("The design matrix is full rank and orthogonal and coordinate descent starting from the OLS solution will always converge to the global minimum, so there is no need to repeat coordinate descent from random starting values. The value of rand.restart should be reset to zero.\n")}
  }

  lambda <- sigma.sq*lambda

  bb.tmp <- rep(NA, p)
  obj.tmp <- Inf

  iter <- ifelse(rand.restart == 0, 1, rand.restart + 1)
  if (return.obj.iter) {
    objs <- numeric(iter)
    iters <- numeric(iter)
  }
  obj <- matrix(NA, nrow = max.iter, ncol = iter)
  for (m in 1:iter) {
    if (m == 1) {
      if (is.null(start)) {
        bb <- crossprod(solve(Q + ridge.eps*diag(p)), l)
      } else {
        bb <- start
      }
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

      for (i in order) {
        X.ii <- Q[i, i]
        z.k.i <- (crossprod(X[, i], y - crossprod(t(X[, -i]), bb[-i])))/X.ii
        lambda.ii <- lambda/X.ii
        b.kp.i <- powThresh(z.k.i, lambda = lambda.ii, q = q)
        if (is.na(b.kp.i)) {return(rep(NA, p))}
        rr <- rr - (b.kp.i - bb[i])*X[, i]
        bb[i] <- b.kp.i
      }

      obj[k, m] <- powObj(beta = bb, sigmasq = 1, lambda = lambda, q = q, Q = Q, l = l, yty = yty)/n
      if (k > 1) {
        obj.diff <- obj[k, m] - obj[k - 1, m]
        if (print.iter) {
          cat("Objective=", round(obj[k, m], 7), "\n")
        }
        if (abs(obj.diff) < tol) {
          opt.cond <- TRUE
        }
      }


      k <- k + 1

    }
    if (return.obj.iter) {
      if (k < max.iter) {
        objs[m] <- obj[k-1, m]
        iters[m] <- k-1
      } else {
        objs[m] <- NA
        iters[m] <- NA
      }
    }
    if (!opt.cond) {
      bb <- rep(Inf, p)
      obj.bb <- Inf
    } else if (iter > 1) {
      obj.bb <- obj[k - 1, m]
      if (obj.bb <= obj.tmp) {
        bb.tmp <- bb
        obj.tmp <- obj.bb
      }
    } else if (iter == 1) {
      bb.tmp <- bb
    }
  }
  if (is.infinite(bb.tmp[1])) {bb.tmp <- rep(NA, p)}

  if (!return.obj.iter) {
    return(bb.tmp)
  } else {
    return(list("opt.b" = bb.tmp, "iter" = iters, "obj" = objs))
  }
}

#' Pathwise coordinate descent for penalized regression with power penalty
#'
#' @name from.zero
#'
#' @description Gives the value of of the length \eqn{p} vector \eqn{\beta} that minimizes: \cr
#' \deqn{(y - X\beta)^2/(2\sigma^2) + \lambda ||\beta||^q_q} \cr
#' for fixed \eqn{y}, \eqn{X}, \eqn{\sigma^2 > 0}, \eqn{\lambda > 0 0} and \eqn{q > 0} based on pathwise coordinate descent. \cr \cr
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
#' @usage \code{from.zero()}
#'
#' @param \code{X} design matrix
#' @param \code{y} response vector
#' @param \code{sigma.sq} scalar value \eqn{\sigma^2}
#' @param \code{tau.sq} scalar value \eqn{\tau^2}
#' @param \code{q} scalar value \eqn{q}
#' @param \code{max.iter} maximum number of iterations for coordinate descent
#' @param \code{print.iter} logical value indicating whether iteration count for coordinate descent should be printed
#' @param \code{tol} scalar tolerance value for assessing convergence of objective function
#' @param \code{order} Vector indicating order of variables in coordinate descent, defaults to order 1:p
#' @param \code{return.final} TRUE/FALSE value indicating whether or the order of variables in coordinate descent should be randomly assigned, defaults to FALSE with order 1:p
#'
#' @return Returns a vector of optimal values for \eqn{\beta} based on pathwise coordinate descent. If the coordinate descent algorithm does not meet the optimality conditions given in Marjanovic and Solo (2014), a vector of \code{NA}'s is returned.
#'
#' Note that a non-\code{NA} solution for \eqn{\beta} for any value of \eqn{q} guarantees the global minimum has been attained when \eqn{X} is full rank and orthogonal. Otherwise, when \eqn{q \le 1} the solution will only correspond to a global minimum when the conditions on \eqn{X} given in Marjanovic and Solo (2014) are satisfied.
#'
#' @source For \eqn{q\le 1}, uses coordinate descent algorithm given by Marjanovic and Solo (2014), modified to accomodate X that do not have standardized columns. \cr
#'
#' Box, G. E. P., and G. C. Tiao. "Bayesian inference in statistical inference." Adison-Wesley, Reading, Mass (1973). \cr
#'
#' Marjanovic, Goran, and Victor Solo. "\eqn{l_q} Sparsity Penalized Linear Regression With Cyclic Descent." IEEE Transactions on Signal Processing 62.6 (2014): 1464-1475. \cr
#'
#' @export
from.zero <- function(X, y, sig.sq, tau.sq, q,
                      num.seq = 100,
                      print.iter = FALSE,
                      order = 1:p, estimate.only = TRUE, max.iter = 10000,
                      tol = 10^(-7), warm.start = TRUE) {

  if (num.seq == 1) {
    cat("Need to provide num.seq > 1\n")
    break;
  }

  p <- ncol(X)
  b <- numeric(p)
  for (i in 1:length(b)) {
    b[i] <- crossprod(X[, i], y)
  }

  tau.sq.min <- (sqrt(gamma(3))/max(abs(b)/sig.sq))^2

  tau.sq.seq <- exp(seq(log(tau.sq.min), log(tau.sq), length.out = num.seq))

  betas <- matrix(nrow = num.seq, ncol = p)
  objs <- numeric(num.seq)
  times <- iters <- numeric(num.seq)

  for (i in 1:num.seq) {
    if(print.iter) {cat("i = ", i, "\n")}
    if (warm.start) {
      start <- rep(0, p)
      if (i > 1) {
        start <- betas[i - 1, ]
      }
    } else {
      start <- rnorm(p)
    }
    ti <- system.time(
    CD <- powCD(X, y, sigma.sq = sig.sq, lambda = (gamma(3/q)/gamma(1/q))^(1/2)/sqrt(tau.sq.seq[i])^q,
                q = q, start = start,
                rand.restart = 0, return.obj.iter = TRUE, order = order, max.iter = max.iter, tol = tol)
    )
    betas[i, ] <- CD[["opt.b"]]
    objs[i] <- CD[["obj"]]
    iters[i] <- CD[["iter"]]
    times[i] <- ti[3]
  }
  obj.bb <- objs[num.seq]
  betas.return <- betas
  objs.return <- objs
  iters.return <- iters

  if (estimate.only) {
    res <- list("betas" = betas.return[num.seq, ])
  } else {
    res <- list("betas" = betas.return, "objs" = objs.return, "iters" = iters.return,
                "tau.sqs" = tau.sq.seq, "times" = times.return)
  }

  return(res)

}

from.two <- function(X, y, sig.sq, tau.sq, q,
                     num.seq = 100,
                     print.iter = FALSE,
                     order = 1:p, estimate.only = TRUE, max.iter = 10000,
                     tol = 10^(-7), warm.start = TRUE) {

  if (num.seq == 1) {
    cat("Need to provide num.seq > 1\n")
    break;
  }

  p <- ncol(X)

  q.seq <- exp(seq(log(2), log(q), length.out = num.seq))


  betas <- matrix(nrow = num.seq, ncol = p)
  objs <- numeric(num.seq)
  times <- iters <- numeric(num.seq)

  for (i in 1:num.seq) {
    if(print.iter) {cat("i = ", i, "\n")}
    if ((i > 1 & warm.start) | (!warm.start)) {
      if (warm.start) {
        start <- betas[i - 1, ]
      } else {
        start <- rnorm(p)
      }
      ti <- system.time(
      CD <- powCD(X, y, sigma.sq = sig.sq, lambda = (gamma(3/q.seq[i])/gamma(1/q.seq[i]))^(q.seq[i]/2)/sqrt(tau.sq)^q.seq[i],
                  q = q.seq[i], start = start,
                  rand.restart = 0, return.obj.iter = TRUE, order = order, tol = tol, max.iter = max.iter))
      betas[i, ] <- CD[["opt.b"]]
      objs[i] <- CD[["obj"]]
      iters[i] <- CD[["iter"]]
      times[i] <- ti[3]
    } else {
      betas[i, ] <- crossprod(solve(crossprod(X)/sig.sq + diag(p)/tau.sq), crossprod(X, y)/sig.sq)
      iters[i] <- 1
    }
  }


  obj.bb <- objs[num.seq]
  betas.return <- betas
  objs.return <- objs
  iters.return <- iters
  times.return <- times

  if (estimate.only) {
    res <- list("betas" = betas.return[num.seq, ])
  } else {
    res <- list("betas" = betas.return, "objs" = objs.return, "iters" = iters.return,
                "qs" = q.seq, "times" = times.return)
  }

  return(res)

}

reg.surface <- function(X, y, sig.sq, tau.sq, q,
                        num.seq.q = 100,
                        num.seq.tau.sq = 100,
                        print.iter = FALSE,
                        order = 1:p, estimate.only = TRUE, max.iter = 10000,
                        tol = 10^(-7), rand.restart = 0, warm.start.q = TRUE,
                        warm.start.tau.sq = TRUE) {

  if (num.seq.q == 1 | num.seq.tau.sq == 1) {
    cat("Need to provide num.seq.q > 1 and num.seq.tau.sq > 1\n")
    break;
  }

  q.seq <- exp(seq(log(2), log(q), length.out = num.seq.q))

  p <- ncol(X)
  b <- numeric(p)
  for (i in 1:length(b)) {
    b[i] <- crossprod(X[, i], y)
  }

  tau.sq.min <- (sqrt(gamma(3))/max(abs(b)/sig.sq))^2

  tau.sq.seq <- exp(seq(log(tau.sq.min), log(tau.sq), length.out = num.seq.tau.sq))

  betas <- array(dim = c(num.seq.q, num.seq.tau.sq, p))
  objs <- array(Inf, dim = c(num.seq.q, num.seq.tau.sq))
  times <- iters <- array(dim = c(num.seq.q, num.seq.tau.sq))
  if (warm.start.q | !(warm.start.q & warm.start.tau.sq)) {
    for (i in 1:num.seq.tau.sq) {
      ft <- from.two(y = y, X = X, tau.sq = tau.sq.seq[i],
                     sig.sq = sig.sq, q = q,
                     num.seq = num.seq.q, print.iter = FALSE,
                     estimate.only = FALSE, warm.start = warm.start.q,
                     tol = 10^(-14))
      betas[, i, ] <- ft$betas
      objs[, i] <- ft$objs
      iters[, i] <- ft$iters
      times[, i] <- ft$times

    }
  } else if (warm.start.tau.sq) {
    for (i in 1:num.seq.q) {
      ft <- from.zero(y = y, X = X, tau.sq = tau.sq,
                      sig.sq = sig.sq, q = q.seq[i],
                      num.seq = num.seq.q, print.iter = FALSE,
                      estimate.only = FALSE, warm.start = warm.start.tau.sq,
                      tol = 10^(-14))
      betas[i, , ] <- ft$betas
      objs[i, ] <- ft$objs
      iters[i, ] <- ft$iters
      times[i, ] <- ft$times
    }
  }

  obj.bb <- objs[num.seq.q, num.seq.tau.sq]
  betas.return <- betas
  objs.return <- objs
  iters.return <- iters
  times.return <- times

  if (estimate.only) {
    res <- list("betas" = betas.return[num.seq.q, num.seq.tau.sq, ])
  } else {
    res <- list("betas" = betas.return, "objs" = objs.return, "iters" = iters.return,
                "qs" = q.seq, "tau.sqs" = tau.sq.seq, "times" = times.return)
  }

  return(res)

}

