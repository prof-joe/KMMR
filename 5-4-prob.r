# 第4章 カーネル計算の実際（問題46～64）

## 48

kernel.pca.train <- function(x, k) {
  # データ x とカーネル k からグラム行列を求める。
  res <- eigen(K)
  alpha <- matrix(0, n, n)
  for (i in 1:n)
    alpha[, i] <- res$vector[, i] / res$value[i] ^ 0.5
  return(alpha)
}

kernel.pca.test <- function(x, k, alpha, m, z) {
  # x, k, alpha, m, z から m 次までのスコア pca を求める
  return(pca)
}

sigma.2 <- 0.01

k <- function(x, y) {
  return(exp(-norm(x-y, "2")^2 / 2 / sigma.2))
}

x <- as.matrix(USArrests)
n <- nrow(x)
p <- ncol(x)
alpha <- kernel.pca.train(x, k)
z <- array(dim = c(n, 2))
for (i in 1:n)
  z[i, ] <- kernel.pca.test(x, k, alpha, 2, x[i, ])
min.1 <- min(z[, 1])
min.2 <- min(z[, 2])
max.1 <- max(z[, 1])
max.2 <- max(z[, 2]) 
plot(0, xlim = c(min.1, max.1), ylim = c(min.2, max.2),
     xlab = "First", ylab = "Second", cex.lab = 0.75, cex.axis = 0.75,
     main = "Kernel PCA (Gauss 0.01)")
for (i in 1:n)
  if (i != 5)
    text(z[i, 1], z[i, 2], labels = i, cex = 0.5)
text(z[5, 1], z[5, 2], 5, col = "red")


## 54

sigma <- 10
sigma2 <- sigma ^ 2

z <- function(x) {
  return(sqrt(2/m) * cos(w*x + b))
}

zz = function(x, y) {
  return(sum(z(x)*z(y)))
}


## 61

alpha.m <- function(k, x, y, m) {
  n <- length(x)
  K <- matrix(0, n, n)
  for (i in 1:n)
    for (j in 1:n)
      K[i, j] <- k(x[i], x[j])
  A <- svd(K[1:m, 1:m])
  u <- array(dim = c(n, m))
  for (i in 1:m)
    for (j in 1:n)
      u[j, i] <- sqrt(m/n) * sum(K[j, 1:m] * A$u[1:m, i]) / A$d[i]
  mu <- A$d * n / m
  R <- sqrt(mu[1]) * u[, 1]
  for (i in 2:m)
    R <- cbind(R, sqrt(mu[i]) * u[, i])
  alpha <- (diag(n) - R %*% solve(t(R) %*% R + lambda * diag(m)) %*% t(R)) %*%
    y / lambda
  return(as.vector(alpha))
}
