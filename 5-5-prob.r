# Ch.5　MMD and HSIC（Problem 65～83）


## 76

cc <- function(x, y) {
  return(sum(x*y) / length(x))
}

f <- function(u, v) {
  return(u - cc(u, v) / cc(v, v) * v)
}

## Estimate the top
cc <- function(x, y) {
  return(sum(x*y) / length(x))
}

f <- function(u, v) {
  return(u - cc(u, v) / cc(v, v) * v)
}

x.y <- f(x, y)
y.z <- f(y, z)
z.x <- f(z, x)
x.z <- f(x, z)
z.y <- f(z, y)
y.x <- f(y, x)
v1 <- HSIC.2(x, y.x, z.x, k.x, k.y, k.z)
v2 <- HSIC.2(y, z.y, x.y, k.y, k.z, k.x)
v3 <- HSIC.2(z, x.z, y.z, k.z, k.x, k.y)
if (v1 < v2) {
  if (v1 < v3) {
    top <- 1
  } else {
    top <- 3
  }
} else {
  if (v2 < v3) {
    top <- 2
  } else {
    top <- 3
  }
}

## Estimate the bottom
x.yz <- f(x.y, z.y)
y.zx <- f(y.z, x.z)
z.xy <- f(z.x, y.x)
if (top == 1) {
  v1 <- ## blank (1) ##
  v2 <- ## blank (2) ##
  if (v1 < v2) {
    middle <- 2
    bottom <- 3
  } else {
    middle <- 3
    bottom <- 2
  }
}
if (top == 2) {
  v1 <- ## blank (3) ##
  v2 <- ## blank (4) ##
  if (v1 < v2) {
    middle <- 3
    bottom <- 1
  } else {
    middle <- 1
    bottom <- 3
  }
}
if (top == 3) {
  v1 <- ## blank (5) ##
  v2 <- ## blank (6) ##
  if (v1 < v2) {
    middle <- 1
    bottom <- 2
  } else {
    middle <- 2
    bottom <- 1
  }
}

## Output results
print(paste("top = ", top))
print(paste("middle = ", middle))
print(paste("bottom = ", bottom))


## 77

x <- rnorm(n)
y <- rnorm(n)
u <- HSIC.1(x, y, k.x, k.y)
m <- 100
w <- NULL
for (i in 1:m) {
  x <- x[sample(n, n)]
  w <- c(w, HSIC.1(x, y, k.x, k.y))
}
v <- quantile(w, 0.95)
plot(density(w), xlim = c(min(w, v, u), max(w, v, u)))
abline(v = v, col = "red", lty = 2, lwd = 2)
abline(v = u, col = "blue", lty = 1, lwd = 2)

