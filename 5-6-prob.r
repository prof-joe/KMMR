# 第6章 Gauss過程と関数データ解析（問題84～100）


## 86
gp.2 <- function(x.pred) {
  h <- array(dim = n)
  for(i in 1:n)
    h[i] <- k(x.pred, x[i])
  L <- chol(K + sigma.2 * diag(n))
  alpha <- solve(L, solve(t(L), y - mu(x)))
  mm <- mu(x.pred) + sum(t(h) * alpha)
  gamma <- solve(t(L), h)
  ss <- k(x.pred, x.pred) - sum(gamma ^ 2)
  return(list(mm = mm, ss = ss))
}
