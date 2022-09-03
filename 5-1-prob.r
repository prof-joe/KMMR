# Ch.1 Positive Definite Kernel（Problem 1～15）

## 2
k <- function(x, y, lambda) {
  return(D(abs(x-y) / lambda))
}

n <- 250
x <- 2 * rnorm(n)
y <- sin(2*pi*x) + rnorm(n) / 4
plot(seq(-3, 3, length = 10), seq(-2, 3, length = 10), type = "n",
     xlab = "x", ylab = "y")
points(x, y)
xx <- seq(-3, 3, 0.1)
yy <- NULL
for (zz in xx)
  yy <- c(yy, f(zz, 0.05))
lines(xx, yy, col = "green")
yy <- NULL
for (zz in xx)
  yy <- c(yy, f(zz, 0.50))
lines(xx, yy, col = "red")
title("Nadaraya-Watson Estimate")
legend("topleft", legend = paste0("lambda = ", c(0.05, 0.35, 0.50)),
       lwd = 1, col = c("green", "blue", "red"))

## 8

k <- function(x, y, sigma2) {
  return(exp(-(x-y)^2 / 2 / sigma2))
}

# Data Generation
n <- 100
x <- 2 * rnorm(n)
y <- sin(2*pi*x) + rnorm(n) / 4

m <- n / 10
sigma2.seq <- seq(0.001, 0.01, 0.001)
SS.min <- Inf
for (sigma2 in sigma2.seq) {
  SS <- 0
  for (h in 1:10) {
    test <- ((h-1)*m + 1):(h*m)
    train <- setdiff(1:n, test)
    for (j in test) {
      u <- 0
      v <- 0
      for (i in train) {
        kk <- k(x[i], x[j], sigma2)
        u <- u + kk * y[i]
        v <- v + kk
      }
      if (v != 0) {
        z <- u / v
        SS <- SS + (y[j]-z)^2
      }
    }
  }
  if (SS < SS.min) {
    SS.min <- SS
    sigma2.best <- sigma2
  }
}
paste0("Best sigma2 = ", sigma2.best)
plot(seq(-3, 3, length = 10), seq(-2, 3, length = 10), type = "n",
     xlab = "x", ylab = "y")
points(x, y)
xx <- seq(-3, 3, 0.1)
yy <- NULL
for (zz in xx)
  yy <- c(yy, f(zz, sigma2.best))
lines(xx, yy, col = "red")
title("Nadaraya-Watson estimate")


## 12

string.kernel <- function(x, y) {
  m <- nchar(x)
  n <- nchar(y)
  S <- 0
  for (i in 1:m)
    for (j in i:m)
      for (k in 1:n)
        if (substring(x, i, j) == substring(y, k, k+j-i))
          S <- S + 1
  return(S)
}


## 15

k <- function(s, p) {
  return(rob(s, p) / length(node))
}

prob <- function(s, p) {
  if (length(node[s[1]]) == 0)
    return(0)
  if (length(s) == 1)
    return(p)
  m <- length(s)
  S <- (1 - p) / length(node[s[1]]) * prob(s[2:m], p)
  return(S)
}
