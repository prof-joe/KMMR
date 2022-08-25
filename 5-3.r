# 第3章 再生核Hilbert空間

## 3.1 RKHS

## 3.2 Sobolev空間

## 3.3 Mercerの定理

### 例59

Hermite=function(j){   ##  R言語は添字が1から
  if(j==0)return(1)
  a=rep(0,j+2); b=rep(0,j+2)
  a[1]=1
  for(i in 1:j){
    b[1]=-a[2]
    for(k in 1:(i+1))b[k+1]=2*a[k]-(k+1)*a[k+2]
    a=b
  }
  return(b[1:(j+1)])   ## Hermite多項式の係数を出力
}

Hermite(2)             ## 2次のHermite多項式
Hermite(3)             ## 3次のHermite多項式
Hermite(4)             ## 4次のHermite多項式

H=function(j,x){
  coef=Hermite(j)
  S=0
  for(i in 0:j) S=S+coef[i+1]*x^i
  return(S)
}

cc=sqrt(5)/4; a=1/4
phi=function(j,x) exp(-(cc-a)*x^2)*H(j,sqrt(2*cc)*x)
curve(phi(0,x),-2,2, ylim=c(-2,8),col=1,ylab="phi")
for(i in 1:3)curve(phi(i,x),-2,2, ylim=c(-2,8), add=TRUE, ann=FALSE, col=i+1)
legend("topright",legend=paste("j=",0:3),lwd=1, col=1:4)
title("Gaussカーネルの固有関数")


### 例62

## カーネルの定義
sigma=1; k=function(x,y)exp(-(x-y)^2/sigma^2)

## mサンプルの発生とグラム行列の設定
m=300; x=rnorm(m)-2*rnorm(m)^2+3*rnorm(m)^3

## 固有値・固有ベクトルの計算
K=matrix(0,m,m)
for(i in 1:m)for(j in 1:m)K[i,j]=k(x[i],x[j])
  eig=eigen(K)
  lam.m=eig$values
  lam=lam.m/m
  U=eig$vector
  alpha=array(0,dim=c(m,m))
  for(i in 1:m)alpha[,i]=U[,i]*sqrt(m)/lam.m[i]
  
## グラフの表示
F=function(y,i){
  S=0; for(j in 1:m)S=S+alpha[j,i]*k(x[j],y)
  return(S)
}
i=1  ## i の値を変えて実行する。
G=function(y)F(y,i)
plot(G,xlim=c(-2,2))
title("Eigen Values and their Eigen Functions")
