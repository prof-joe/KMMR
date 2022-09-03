# Ch.3　Reproducing Kernel Hilbert Space

## 3.1 RKHS

## 3.2 Sobolev Space

## 3.3 Mercer's theorem

### Example 59

Hermite=function(j){   ##  The R language has subscripts from 1
  if(j==0)return(1)
  a=rep(0,j+2); b=rep(0,j+2)
  a[1]=1
  for(i in 1:j){
    b[1]=-a[2]
    for(k in 1:(i+1))b[k+1]=2*a[k]-(k+1)*a[k+2]
    a=b
  }
  return(b[1:(j+1)])   ## Output coefficients of Hermite polynomial
}

Hermite(2)             ## Hermite polynomial of the second degree
Hermite(3)             ## Hermite polynomial of the third degree
Hermite(4)             ## Hermite polynomial of the fourth degree

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
title("Characteristic function of Gauss Kernel")


### Example 62

## Definition of Kernel
sigma=1; k=function(x,y)exp(-(x-y)^2/sigma^2)

## Sample generation and gram matrix setup
m=300; x=rnorm(m)-2*rnorm(m)^2+3*rnorm(m)^3

## Calculation of eigenvalues and eigenvectors
K=matrix(0,m,m)
for(i in 1:m)for(j in 1:m)K[i,j]=k(x[i],x[j])
  eig=eigen(K)
  lam.m=eig$values
  lam=lam.m/m
  U=eig$vector
  alpha=array(0,dim=c(m,m))
  for(i in 1:m)alpha[,i]=U[,i]*sqrt(m)/lam.m[i]
  
## Display graph
F=function(y,i){
  S=0; for(j in 1:m)S=S+alpha[j,i]*k(x[j],y)
  return(S)
}
i=1  ## Run with different i's
G=function(y)F(y,i)
plot(G,xlim=c(-2,2))
title("Eigen Values and their Eigen Functions")
