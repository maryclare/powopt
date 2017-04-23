pow.obj <- function(beta, X, y, sigma.sq.z, lambda, q) {
  log.lik <- -(1/(2*sigma.sq.z))*(crossprod(crossprod(t(X), beta)) - 2*crossprod(beta, crossprod(X, y))) - lambda*sum(abs(beta)^q)
  return(-1*log.lik)
}

# Univariate thresholding function for q < 1
pow.thres <- function(z, lambda, q) {
  
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

pow.cd <- function(X, y, sigma.sq.z, lambda, q, max.iter = 10000,
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
        b.kp.i <- pow.thres(z.k.i, lambda = lambda.ii, q = q)
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
    obj.bb <- pow.obj(beta = bb, X = X, y = y, sigma.sq.z = 1, lambda = lambda, q = q)
    if (obj.bb <= obj.tmp) {
      bb.tmp <- bb
      obj.tmp <- obj.bb
    }
  }
  if (is.infinite(bb.tmp[1])) {bb.tmp <- rep(NA, p)}
  return(bb.tmp)
}

sigma.sq.beta <- 1
zs <- seq(0, 5, by = 0.0001)
plot(zs, pow.thres(zs, lambda = (sigma.sq.beta*gamma(1/q)/gamma(3/q))^(-q/2), q = 0.6), type = "l")
for (i in 1:4) {
  lines(zs, pow.thres(zs, lambda = (sigma.sq.beta*gamma(1/q)/gamma(3/q))^(-q/2), q = 0.6 + i/10), lty = i + 1, col = i + 1)
}
legend("topleft", legend = 0.6 + 0:4/10, col = 1:5, lty = 1:5)

# Replicate Mazumder
qs <- c(0.001, 0.3, 0.7, 1)
zs <- seq(0, 3, by = 0.0001)
plot(zs, pow.thres(zs, lambda = 1, q = qs[1]), type = "l", ylim = c(0, 3))
for (i in 2:4) {
  lines(zs, pow.thres(zs, lambda = 1, q = qs[i]), lty = i, col = i)
}
legend("topleft", col = 1:4, lty = 1:4, legend = qs)

n <- 100
p <- 50
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
# for (i in 1:ncol(X)) {
#   X[, i] <- X[, i]/sqrt(sum(X[, i]^2))
# }
# X <- svd(X)$u%*%diag(svd(X)$d)
# beta <- c(5, 0, 0)
beta <- rnorm(p)
sigma.sq.beta <- 1
y <- X%*%beta + rnorm(n, sd = sqrt(sigma.sq.z))
tol <- 10^(-7)
q <-  100
max.iter = 1000
print.iter = TRUE
ridge.start = FALSE
rand.restart = 10
lambda <- (sigma.sq.beta*gamma(1/q)/gamma(3/q))^(-q/2)

test <- pow.cd(X = X, y = y, sigma.sq.z = sigma.sq.z, lambda = lambda, q = q, 
               max.iter = max.iter, print.iter = FALSE,
               tol = tol, ridge.start = ridge.start, rand.restart = rand.restart)
cbind(test, solve(t(X)%*%X + diag(p)/sigma.sq.beta)%*%t(X)%*%y,
      solve(t(X)%*%X)%*%t(X)%*%y)
