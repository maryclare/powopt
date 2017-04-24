## ---- message = FALSE----------------------------------------------------
library(devtools)
install_github("maryclare/powopt")
library(powopt)
set.seed(100)

## ---- fig.show='hold', fig.align = 'center', fig.width = 6, fig.height = 6----
qs <- c(0.001, 0.3, 0.7, 1)
zs <- seq(0, 3, by = 0.001)
plot(zs, powThresh(zs, lambda = 1, q = qs[1]), type = "l", ylim = c(0, 3),
     main = "Threshold Functions", xlab = expression(beta), ylab = "")
for (i in 2:4) {
  lines(zs, powThresh(zs, lambda = 1, q = qs[i]), lty = i, col = i)
}
legend("topleft", col = 1:4, lty = 1:4, legend = qs, cex = 0.75)

## ---- fig.show='hold', fig.align = 'center', fig.width = 6, fig.height = 6----
qs <- c(0.25, 0.5, 1, 1.25, 2, 3)
sigma.sq.beta <- 1
zs <- seq(0, 5, by = 0.001)
plot(zs, powThresh(zs, lambda = (sigma.sq.beta*gamma(1/qs[1])/gamma(3/qs[1]))^(-qs[1]/2), q = qs[1]), type = "l",
     xlab = "z", ylab = expression(beta))
for (i in 2:length(qs)) {
  q <- qs[i]
  lines(zs, powThresh(zs, lambda = (sigma.sq.beta*gamma(1/q)/gamma(3/q))^(-q/2), q = q), lty = i, col = i)
}
legend("topleft", legend = qs, col = 1:length(qs), lty = 1:length(qs), cex = 0.75)

## ---- fig.show='hold', fig.align = 'center', fig.width = 6, fig.height = 6----
n <- 5
p <- 3
X <- svd(matrix(rnorm(n*p), nrow = n, ncol = p))$u
beta <- rnorm(p)
y <- X%*%beta + rnorm(n)
q <-  2
lambda <- (gamma(1/q)/gamma(3/q))^(-q/2) 
fit <- powCD(X = X, y = y, sigma.sq = 1, lambda = lambda, q = q)
plot(fit, solve(crossprod(X) + diag(p))%*%crossprod(X, y), xlab = "Coordinate Descent", ylab = "Closed Form",
     main = expression(paste("Estimate of ", beta, sep = "")))
abline(a = 0, b = 1)

## ---- fig.show='hold', fig.align = 'center', fig.width = 6, fig.height = 6----
q <-  1000
lambda <- (gamma(1/q)/gamma(3/q))^(-q/2) 
fit <- powCD(X = X, y = y, sigma.sq = 1, lambda = lambda, q = q)
plot(fit, solve(crossprod(X))%*%crossprod(X, y), xlab = "Coordinate Descent", ylab = "Closed Form",
     main = expression(paste("Estimate of ", beta, sep = "")))
abline(a = 0, b = 1)

## ---- fig.show='hold'----------------------------------------------------
qs <-  c(1, 0.5)
lambdas <- seq(0.1, 4, length.out = 100)
fits <- array(dim = c(p, length(lambdas), length(qs)))
for (i in 1:length(lambdas)) {
  for (j in 1:length(qs)) {
    fits[, i, j] <- powCD(X = X, y = y, sigma.sq = 1, lambda = lambdas[i], q = qs[j])
  }
}
for (j in 1:length(qs)) { 
  plot(log(lambdas), fits[1, , j], 
       ylim = range(fits), type = "l", 
       xlab = expression(paste("log(", lambda, ")", sep = "")), 
       ylab = expression(beta), main = paste("q=", q, "\n"))
  for (i in 2:p) {
    lines(log(lambdas), fits[i, , j], col = i)
  }
}
legend("topleft", c(expression(beta[1]),
                        expression(beta[2]),
                        expression(beta[3])), col = 1:3, lty = 1, cex = 0.75)

## ---- fig.show='hold', fig.align = 'center', fig.width = 6, fig.height = 6----
n <- 10
p <- 8
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
beta <- rnorm(p)
y <- X%*%beta + rnorm(n)
q <- 0.5
lambda <- (gamma(1/q)/gamma(3/q))^(-q/2) 
fit1 <- powCD(X = X, y = y, sigma.sq = 1, lambda = lambda, q = q, rand.restart = 1)
fit2 <- powCD(X = X, y = y, sigma.sq = 1, lambda = lambda, q = q, rand.restart = 1)
plot(fit1, fit2, xlab = "Coordinate Descent 1", ylab = "Coordinate Descent 2",
     main = expression(paste("Estimate of ", beta, sep = "")))
abline(a = 0, b = 1)

## ---- fig.show='hold', fig.align = 'center', fig.width = 6, fig.height = 6----
n <- 10
p <- 8
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
beta <- rnorm(p)
y <- X%*%beta + rnorm(n)
q <- 0.5
lambda <- (gamma(1/q)/gamma(3/q))^(-q/2) 
rand.restarts <- seq(1, 10, by = 1)
fits <- matrix(nrow = p, ncol = length(rand.restarts)) 
diff <- numeric(length(rand.restarts) - 1)
for (i in 1:length(rand.restarts)) {
  fits[, i] <- powCD(X = X, y = y, sigma.sq = 1, lambda = lambda, q = q, rand.restart = rand.restarts[i])
  if (i > 1) {
    diff[i - 1] <- mean((fits[, i] - fits[, i - 1])^2)
  }
}
plot(rand.restarts[2:length(rand.restarts)], 
     diff, xlab = "# Random Restarts", ylab = "MSE Compared to Previous Estimate")

