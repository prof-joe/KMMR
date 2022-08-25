# 第6章 Gauss過程と関数データ解析

## 6.1 回帰

### 例83

## (m,k)の定義
m=function(x) 0; k=function(x,y) exp(-(x-y)^2/2)
## 関数 gp.sample の定義
gp.sample=function(x,m,k){
  n=length(x)
  m.x=m(x)
  k.xx=matrix(0,n,n); for(i in 1:n)for(j in 1:n)k.xx[i,j]=k(x[i],x[j])
  R=t(chol(k.xx))
  u=rnorm(n)
  return(as.vector(R%*%u+m.x))
}
## 乱数を発生して、共分散行列を生成して k.xx と比較
x=seq(-2,2,1); n=length(x)
r=100; z=matrix(0,r,n); for(i in 1:r)z[i,]=gp.sample(x,m,k)
k.xx=matrix(0,n,n); for(i in 1:n)for(j in 1:n)k.xx[i,j]=k(x[i],x[j])

cov(z)
k.xx


### 例84

## (m,k)の定義
m=function(x) x[,1]-x[,2]
k=function(x,y) exp(-sum((x-y)^2)/2)
## 関数 gp.sample の定義
gp.sample=function(x,m,k){
  n=nrow(x)
  m.x=m(x)
  k.xx=matrix(0,n,n); for(i in 1:n)for(j in 1:n)k.xx[i,j]=k(x[i,],x[j,])
  R=t(chol(k.xx))
  u=rnorm(n)
  return(R%*%u+m.x)
}
## 乱数を発生して、共分散行列を生成して k.xx と比較
n=5; x=matrix(rnorm(n*2),n,n)
r=100; z=matrix(0,r,n); for(i in 1:r)z[i,]=gp.sample(x,m,k)
k.xx=matrix(0,n,n); for(i in 1:n)for(j in 1:n)k.xx[i,j]=k(x[i],x[j])

cov(z)
k.xx

gp.1=function(x.pred){
  h=array(dim=n); for(i in 1:n)h[i]=k(x.pred,x[i])
  R=solve(K+sigma.2*diag(n))              ## O(n^3)の計算
  mm=mu(x.pred)+t(h)%*%R%*%(y-mu(x))
  ss=k(x.pred,x.pred)-t(h)%*%R%*%h
  return(list(mm=mm,ss=ss))
}
gp.2=function(x.pred){
  h=array(dim=n); for(i in 1:n)h[i]=k(x.pred,x[i])
  L=chol(K+sigma.2*diag(n))               ## O(n^3/3) の計算
  alpha=solve(L,solve(t(L),y-mu(x)))      ## O(n^2) の計算
  mm=mu(x.pred)+sum(t(h)*alpha)
  gamma=solve(t(L),h)                         ## O(n^2) の計算
  ss=k(x.pred,x.pred)-sum(gamma^2)
  return(list(mm=mm,ss=ss))
}


### 例85

sigma.2=0.2     
k=function(x,y)exp(-(x-y)^2/2/sigma.2)  # 共分散関数
mu=function(x) x                        # 平均関数
n=1000; x=runif(n)*6-3; y=sin(x/2)+rnorm(n)  # データ生成
K=array(dim=c(n,n)); for(i in 1:n)for(j in 1:n)K[i,j]=k(x[i],x[j])
## 実行時間を測定
library(tictoc)
tic(); gp.1(0); toc()
tic(); gp.2(0); toc()
# 平均の前後で 3 sigma の幅も記載
u.seq=seq(-3,3,0.1); v.seq=NULL; w.seq=NULL;
for(u in u.seq){res=gp.1(u); v.seq=c(v.seq,res$mm); w.seq=c(w.seq,sqrt(res$ss))}
plot(u.seq,v.seq,xlim=c(-3,3),ylim=c(-3,3),type="l")
lines(u.seq,v.seq+3*w.seq,col="blue"); lines(u.seq,v.seq-3*w.seq,col="blue")
points(x,y)
## サンプルを変えて 5 回
plot(0,xlim=c(-3,3),ylim=c(-3,3),type="n")
n=100
for(h in 1:5){
  x=runif(n)*6-3; y=sin(pi*x/2)+rnorm(n)
  sigma2=0.2
  K=array(dim=c(n,n)); for(i in 1:n)for(j in 1:n)K[i,j]=k(x[i],x[j])
  u.seq=seq(-3,3,0.1); v.seq=NULL
  for(u in u.seq){res=gp.1(u); v.seq=c(v.seq,res$mm)}
  lines(u.seq,v.seq,col=h+1)
}


## 6.2 分類

### 例86

## Iris データ
df=iris
x=df[1:100,1:4]
y=c(rep(1,50),rep(-1,50))
n=length(y)
## 4 個の共変量でカーネルを計算
k=function(x,y) exp(sum(-(x-y)^2)/2)
K=matrix(0,n,n)
for(i in 1:n)for(j in 1:n)K[i,j]=k(x[i,],x[j,])
eps=0.00001
f=rep(0,n)
g=rep(0.1,n)
while(sum((f-g)^2)>eps){
  g=f     ## 比較のため、更新前の値を保存する
  v=exp(-y*f)
  u=y*v/(1+v)
  w=as.vector(v/(1+v)^2)
  W=diag(w); W.p=diag(w^0.5); W.m=diag(w^(-0.5))
  L=chol(diag(n)+W.p%*%K%*%W.p)
  L=t(L)  ## R 言語の chol 関数は転置した行列を出力する
  gamma=W%*%f+u
  beta=solve(L,W.p%*%K%*%gamma)
  alpha=solve(t(L)%*%W.m,beta)
  f=K%*%(gamma-alpha)
}
as.vector(f)

pred=function(z){
    kk=array(0,dim=n); for (i in 1:n)kk[i]=k(z,x[i,])
    mu=sum(kk*as.vector(u))      ## 平均
    alpha=solve(L,W.p%*%kk); sigma2=k(z,z)-sum(alpha^2)    ## 分散
    m=1000; b=rnorm(m,mu,sigma2); pi=sum((1+exp(-b))^(-1))/m   ## 予測値
    return(pi)
}


### 例87

z=array(0,dim=4)
for(j in 1:4)z[j]=mean(x[1:50,j])
pred(z)

for(j in 1:4)z[j]=mean(x[51:100,j])
pred(z)


## 6.3 補助変数法

### 例88

sigma.2=0.05     #本来は推定すべき
k=function(x,y)exp(-(x-y)^2/2/sigma.2)  # 共分散関数
mu=function(x) x                        # 平均関数
n=200; x=runif(n)*6-3; y=sin(x/2)+rnorm(n)  # データ生成
eps=10^(-6)

m=100
index=sample(1:n, m, replace=FALSE)
z=x[index]
m.x=0
m.z=0
K.zz=array(dim=c(m,m)); for(i in 1:m)for(j in 1:m)K.zz[i,j]=k(z[i],z[j])
K.xz=array(dim=c(n,m)); for(i in 1:n)for(j in 1:m)K.xz[i,j]=k(x[i],z[j])
K.zz.inv=solve(K.zz+diag(rep(10^eps,m)))
lambda=array(dim=n)
for(i in 1:n)lambda[i]=k(x[i],x[i])-K.xz[i,1:m]%*%K.zz.inv%*%K.xz[i,1:m]
Lambda.0.inv=diag(1/(lambda+sigma.2))
Q=K.zz+t(K.xz)%*%Lambda.0.inv%*%K.xz    ## Q の計算は、O(n^3) を要求しない
Q.inv=solve(Q+diag(rep(eps,m)))
muu=Q.inv%*%t(K.xz)%*%Lambda.0.inv%*%(y-m.x)
dif=K.zz.inv-Q.inv
K=array(dim=c(n,n)); for(i in 1:n)for(j in 1:n)K[i,j]=k(x[i],x[j])
R=solve(K+sigma.2*diag(n))              ## O(n^3) の計算が必要

gp.ind=function(x.pred){
  h=array(dim=m); for(i in 1:m)h[i]=k(x.pred,z[i])
  mm=mu(x.pred)+h%*%muu
  ss=k(x.pred,x.pred)-h%*%dif%*%h
  return(list(mm=mm,ss=ss))
}                                       ## 補助変数法を用いる

gp.1=function(x.pred){
  h=array(dim=n); for(i in 1:n)h[i]=k(x.pred,x[i])
  mm=mu(x.pred)+t(h)%*%R%*%(y-mu(x))
  ss=k(x.pred,x.pred)-t(h)%*%R%*%h
  return(list(mm=mm,ss=ss))
}　　　　　　　　　　　　　　　　　　　　## 補助変数法を用いない

x.seq=seq(-2,2,0.1)
mmv=NULL; ssv=NULL
for(u in x.seq){
   mmv=c(mmv,gp.ind(u)$mm)
   ssv=c(ssv,gp.ind(u)$ss)
}
plot(0, xlim=c(-2,2),ylim=c(min(mmv),max(mmv)),type="n")
lines(x.seq,mmv,col="red")
lines(x.seq,mmv+3*sqrt(ssv),lty=3,col="red")
lines(x.seq,mmv-3*sqrt(ssv),lty=3,col="red")

x.seq=seq(-2,2,0.1)
mmv=NULL; ssv=NULL
for(u in x.seq){
   mmv=c(mmv,gp.1(u)$mm)
   ssv=c(ssv,gp.1(u)$ss)
}

lines(x.seq,mmv,col="blue")
lines(x.seq,mmv+3*sqrt(ssv),lty=3,col="blue")
lines(x.seq,mmv-3*sqrt(ssv),lty=3,col="blue")
points(x,y)


## 6.4 Karhunen-Loeve展開

### 例89

lambda=function(j) 4/((2*j-1)*pi)^2           ## 固有値
ee=function(j,x) sqrt(2)*sin((2*j-1)*pi/2*x)  ## 固有関数の定義
n=10; m=7
f=function(z,x){                              ## ガウス過程の定義
  n=length(z)
  S=0; for(i in 1:n)S=S+z[i]*ee(i,x)*sqrt(lambda(i))
  return(S)
}
plot(0,xlim=c(-3,3),ylim=c(-2,2),type="n",xlab="x",ylab="f(omega,x)")
for(j in 1:m){
  z=rnorm(n)
  x.seq=seq(-3,3,0.001)
  y.seq=NULL; for(x in x.seq)y.seq=c(y.seq,f(z,x))
  lines(x.seq,y.seq,col=j)
}
title("Brown Motion")

matern=function(nu,l,r){
  p=nu-1/2
  S=0
  for(i in 0:p)S=S+gamma(p+i+1)/gamma(i+1)/gamma(p-i+1)*(sqrt(8*nu)*r/l)^(p-i)
  S=S*gamma(p+2)/gamma(2*p+1)*exp(-sqrt(2*nu)*r/l)
  return(S)
}


### 例90

m=10
l=0.1
for(i in 1:1)curve(matern(i-1/2,l,x),0,0.5,ylim=c(0,10),col=i+1)
for(i in 2:m)curve(matern(i-1/2,l,x),0,0.5,ylim=c(0,10),ann=FALSE,add=TRUE,col=i+1)
legend("topright",legend=paste("nu=",(1:m)+0.5),lwd=2,col=1:m,)
title("Matern Kernel (l=1)")


### 例91

rand.100=function(Sigma){
  L=t(chol(Sigma))     ## 共分散行列 を Cholesky 分解  
  u=rnorm(100)
  y=as.vector(L%*%u)   ## 平均 0 共分散行列 の乱数を 1 組生成
}

x = seq(0,1,length=100)
z = abs(outer(x,x,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
l=0.1

Sigma_OU = exp(-z/l)      ## OU:   matern(1/2,l,z) では遅い
y=rand.100(Sigma_OU)

plot(x,y,type="l",ylim=c(-3,3))
for(i in 1:5){
  y = rand.100(Sigma_OU)
  lines(x,y,col=i+1)
}
title("OU process (nu=1/2,l=0.1)")

Sigma_M=matern(3/2,l,z)   ## Matern 
y = rand.100(Sigma_M)
plot(x,y,type="l",ylim=c(-3,3))
for(i in 1:5){
  y = rand.100(Sigma_M)
  lines(x,y,col=i+1)
}
title("Matern process (nu=3/2,l=0.1)")

## 6.5 関数データ解析

### 例92

library(fda)
g=function(j,x){        ##  基底を p 個用意する
  if(j==1) return(1/sqrt(2*pi)) 
  if(j%%2==0) return(cos((j%/%2)*x)/sqrt(pi))
  else return(sin((j%/%2)*x)/sqrt(pi))
}
beta=function(x,y){     ##  関数の p 個の基底の前の係数を計算する
    X=matrix(0,N,p)
    for(i in 1:N)for(j in 1:p)X[i,j]=g(j,x[i])
    beta=solve(t(X)%*%X+0.0001*diag(p))%*%t(X)%*%y
    return(drop(beta))
}
N=365; n=35; m=5; p=100; df=daily
C=matrix(0,n,p)
for(i in 1:n){x=(1:N)*(2*pi/N)-pi; y=as.vector(df[[2]][,i]); C[i,]=beta(x,y)}
res=prcomp(C)
B=res$rotation
xx=res$x

z=function(i,m,x){   ##  p 個の基底のうち、m 主成分で近似したもとの関数
  S=0
  for(j in 1:p)for(k in 1:m)for(r in 1:p)S=S+C[i,j]*B[j,k]*B[r,k]*g(r,x)
  return(S)
}
x.seq=seq(-pi,pi,2*pi/100)
plot(0,xlim=c(-pi,pi),ylim=c(-15,25),type="n",xlab="Days",ylab="Temp(C)", 
main="Reconstruction for each m")
lines(x,df[[2]][,14],lwd=2)
for(m in 2:6){
    lines(x.seq,z(14,m,x.seq),col=m,lty=1)
}
legend("bottom",legend=c("Original",paste("m=",2:6)), lwd=c(2,rep(1,5)), col=1:6,ncol=2)

lambda=res$sdev^2
ratio=lambda/sum(lambda)
plot(1:5,ratio[1:5],xlab="PC1 through PC5",ylab="Ratio",type="l",main="Ratio ")

h=function(coef,x){     ##  係数を用いて関数を定義する
  S=0
  for(j in 1:p)S=S+coef[j]*g(j,x)
  return(S)
}
plot(0,xlim=c(-pi,pi),ylim=c(-1,1),type="n")
for(j in 1:3)lines(x.seq,h(B[,j],x.seq),col=j)

index=c(10,12,13,14,17,24,26,27)
others=setdiff(1:35,index)
first=substring(df[[1]][index],1,1)
plot(0,xlim=c(-25,35),ylim=c(-15,10),type="n",xlab="PC1",ylab="PC2",main="Canadian Weather")
points(xx[others,1],xx[others,2],pch=4)
points(xx[index,1], xx[index,2], pch = first, cex = 1, col =rainbow(8))
legend("bottom",legend=df[[1]][index], pch=first, col=rainbow(8),ncol=2)
