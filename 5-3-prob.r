# 第3章 再生核Hilbert空間（問題31～45）


## 38

H <- function(j, x) {
  if (j == 0) {
    return(1)
  } else if (j == 1) {
    return(2*x)
  } else if (j == 2) {
    return(-2 + 4*x^2)
  } else {
    return(4*x - 8*x^3)
  }
}

cc <- sqrt(5) / 4                                ##
a <- 1 / 4                                       ##

phi <- function(j, x) {                          ###
  return(exp(-(cc-a)*x^2) * H(j, sqrt(2*cc)*x))  ###
}

curve(phi(0, x), -2, 2, ylim = c(-2, 8), col = 1, ylab = "phi")
for (i in 1:3)
  curve(phi(i, x), -2, 2, ylim = c(-2, 8), add = TRUE, ann = FALSE, col = i + 1)
legend("topright", legend = paste("j = ", 0:3), lwd = 1, col = 1:4)
title("Gauss カーネルの固有関数")


## 42

K <- matrix(0, m, m)
for (i in 1:m)
  for (j in 1:m)
    K[i, j] <- k(x[i], x[j])
eig <- eigen(K)
lam.m <- eig$values
lam <- lam.m / m
U <- eig$vector
alpha <- array(0, dim = c(m, m))
for (i in 1:m)
  alpha[, i] <- U[, i] * sqrt(m) / lam.m[i]

F <- function(y, i) {
  S <- 0
  for (j in 1:m)
    S <- S + alpha[j, i] * k(x[j], y)
  return(S)
}
