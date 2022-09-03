# Ch.4　Kernel Computation in Practice

## 4.1 Kernel Ridge regression

alpha=function(k,x,y){
  n=length(x); K=matrix(0,n,n)
  for(i in 1:n)for(j in 1:n)K[i,j]=k(x[i],x[j])
  return(solve(K+10^(-5)*diag(n))%*%y)   ## add 10^(-5) I to K to make it regular
}

### Example 63

k.p=function(x,y) (sum(x*y)+1)^3     ## Definition of Kernel
k.g=function(x,y) exp(-(x-y)^2/2)    ## Definition of Kernel
lambda=0.1
n=50; x=rnorm(n); y=1+x+x^2+rnorm(n)       ## Data Generation
alpha.p=alpha(k.p,x,y); alpha.g=alpha(k.g,x,y)
z=sort(x); u=array(n); v=array(n)
for(j in 1:n){
  S=0;for(i in 1:n)S=S+alpha.p[i]*k.p(x[i],z[j]); u[j]=S
  S=0;for(i in 1:n)S=S+alpha.g[i]*k.g(x[i],z[j]); v[j]=S
}
plot(z,u,type="l",xlim=c(-1,1),xlab="x", ylab="y", ylim=c(-1,5),
  col="red",main="Kernel Regression")
lines(z,v,col="blue"); points(x,y)
legend("topleft", legend = c("Polynomial Kernel","Gauss Kernel"), 
  col = c("red","blue"), lty = 1)


alpha=function(k,x,y){
    n=length(x); K=matrix(0,n,n)
    for(i in 1:n)for(j in 1:n)K[i,j]=k(x[i],x[j])
    return(solve(K+lambda*diag(n))%*%y)
}


### Example 64

k.p=function(x,y) (sum(x*y)+1)^3    ## Definition of Kernel
k.g=function(x,y) exp(-(x-y)^2/2)    ## Definition of Kernel
lambda=0.1
n=50; x=rnorm(n); y=1+x+x^2+rnorm(n)       ## Data Generation
alpha.p=alpha(k.p,x,y); alpha.g=alpha(k.g,x,y)
z=sort(x); u=array(n); v=array(n)
for(j in 1:n){
    S=0;for(i in 1:n)S=S+alpha.p[i]*k.p(x[i],z[j]); u[j]=S
    S=0;for(i in 1:n)S=S+alpha.g[i]*k.g(x[i],z[j]); v[j]=S
}
plot(z,u,type="l",xlim=c(-1,1),xlab="x", ylab="y", ylim=c(-1,5),col="red",main="Kernel Ridge")
lines(z,v,col="blue"); points(x,y)


## 4.2 Kernel Principal Component Analysis


kernel.pca.train=function(x,k){
  n=nrow(x); K=matrix(0,n,n); S=rep(0,n); T=rep(0,n)
  for(i in 1:n)for(j in 1:n)K[i,j]=k(x[i,],x[j,])
  for(i in 1:n)S[i]=sum(K[i,])
  for(j in 1:n)T[j]=sum(K[,j])
  U=sum(K)
  for(i in 1:n)for(j in 1:n)K[i,j]=K[i,j]-S[i]/n-T[j]/n+U/n^2
  res=eigen(K)
  alpha=matrix(0,n,n)
  for (i in 1:n)alpha[,i]=res$vector[,i]/res$value[i]^0.5
  return(alpha)
}

kernel.pca.test=function(x,k,alpha,m,z){
  n=nrow(x)
  pca=array(0,dim=m)
  for(i in 1:n)pca=pca+alpha[i,1:m]*k(x[i,],z)
  return(pca) 
}


### Example 65

#k=function(x,y) sum(x*y)
sigma.2=0.01; k=function(x,y) exp(-norm(x-y,"2")^2/2/sigma.2)
x=as.matrix(USArrests); n=nrow(x); p=ncol(x)
alpha=kernel.pca.train(x,k)
z=array(dim=c(n,2)); for(i in 1:n)z[i,]=kernel.pca.test(x,k,alpha,2,x[i,])
min.1=min(z[,1]); min.2=min(z[,2]); max.1=max(z[,1]); max.2=max(z[,2]) 
plot(0, xlim=c(min.1,max.1),ylim=c(min.2,max.2),xlab="First",ylab="Second",cex.lab=0.75,cex.axis = 0.75, main="Kernel PCA (Gauss 0.01)")
for(i in 1:n)if(i!=5)text(z[i,1],z[i,2],labels=i,cex = 0.5)
text(z[5,1],z[5,2],5,col="red")
## For a normal PCA, the score can be obtained with the following one line.
z=prcomp(x)$x[,1:2]


## 4.3 Kernel SVM

### Example 66

library(quadprog)
K.linear <-function(x,y) {return(t(x)%*%y)}
K.poly <-function(x,y) {return((1+t(x)%*%y)^2)}
svm.2=function(X,y,C, K){ ## Function name is set to svm.2
 eps=0.0001; n=nrow(X); Dmat=matrix(nrow=n,ncol=n);Kmat = matrix(nrow=n,ncol = n)
 for(i in 1:n)for(j in 1:n) {Dmat[i,j]=K(X[i,],X[j,])*y[i]*y[j];Kmat = K(X[i,],X[j,])}
 Dmat=Dmat+eps*diag(n); dvec=rep(1,n)
 Amat=matrix(nrow=(2*n+1),ncol=n); Amat[1,]=y; Amat[2:(n+1),1:n]=-diag(n);
 Amat[(n+2):(2*n+1),1:n]=diag(n) ; Amat=t(Amat)
 bvec=c(0,-C*rep(1,n),rep(0,n)); meq=1
 alpha=solve.QP(Dmat,dvec,Amat,bvec=bvec,meq=1)$solution
 index=(1:n)[0<alpha&alpha<C]
 beta=drop(Kmat%*%(alpha*y)); beta.0=mean(y[index]-beta[index])
 return(list(alpha=alpha,beta.0=beta.0))
}
# Definition of Function 
plot.kernel=function(K, lty){ ## Specify the line type with the parameter lty.
 qq=svm.2(X,y,0.1,K); alpha=qq$alpha; beta.0=qq$beta.0
 f=function(u,v){x=c(u,v); S=beta.0; for(i in 1:n)S=S+alpha[i]*y[i]*K(X[i,], x); return(S)}
 ## is the height at f(x,y). From this we can find the contour.
 u=seq(-2,2,.1);v=seq(-2,2,.1);w=array(dim=c(41,41))
 for(i in 1:41)for(j in 1:41)w[i,j]=f(u[i],v[j])
 contour(u,v,w,level=0,add=TRUE,lty=lty)
}
 # Execution
a=rnorm(1); b=rnorm(1)
n=100; X=matrix(rnorm(n*2),ncol=2,nrow=n); y=sign(a*X[,1]+b*X[,2]+0.3*rnorm(n))
plot(-3:3,-3:3,xlab="X[,1]",ylab="X[,2]", type="n")
for(i in 1:n){
 if(y[i]==1)points(X[i,1],X[i,2],col="red") else points(X[i,1],X[i,2],col= "blue")
}
plot.kernel(K.linear,1); plot.kernel(K.poly,2)


## 4.4 Spline

## Construct function d and h to find the basis
d=function(j,x,knots){
    K=length(knots);
    (max((x-knots[j])^3,0)-max((x-knots[K])^3,0))/(knots[K]-knots[j])
}
h=function(j,x,knots){
    K=length(knots);
    if(j==1) return(1) else if(j==2)return(x) else return(d(j-2,x,knots)-d(K-1,x,knots))
}
## G is the value obtained by integrating a function differentiated twice.
G=function(x){        ## Assuming that each value of x is in ascending order
   n=length(x); g=matrix(0, nrow=n,ncol=n)
   for(i in 3:(n-1))for(j in i:n){
      g[i,j]=12*(x[n]-x[n-1])*(x[n-1]-x[j-2])*(x[n-1]-x[i-2])/(x[n]-x[i-2])/(x[n]-x[j-2])+
      (12*x[n-1]+6*x[j-2]-18*x[i-2])*(x[n-1]-x[j-2])^2/(x[n]-x[i-2])/(x[n]-x[j-2])
      g[j,i]=g[i,j]
   }
return(g)
}
## Main Process
n=100; x=runif(n,-5,5); y=x+sin(x)*2+rnorm(n)  ## Data Generation
index=order(x); x=x[index];y=y[index]
X=matrix(nrow=n,ncol=n); X[,1]=1
for(j in 2:n)for(i in 1:n)X[i,j]=h(j,x[i],x)  ## Generate matrix X
GG=G(x);                                    ## Generate matrix G
lambda.set=c(1,30,80); col.set=c("red","blue","green")
for(i in 1:3){
    lambda=lambda.set[i]
    gamma=solve(t(X)%*%X+lambda*GG)%*%t(X)%*%y
    g=function(u){S=gamma[1]; for(j in 2:n)S=S+gamma[j]*h(j,u,x); return(S)}
    u.seq=seq(-8,8,0.02); v.seq=NULL; for(u in u.seq)v.seq=c(v.seq,g(u))
    plot(u.seq,v.seq,type="l",yaxt="n", xlab="x",ylab="g(x)",ylim=c(-8,8), col=col.set[i])
    par(new=TRUE)
}
points(x,y); legend("topleft", paste0("lambda=",lambda.set), col=col.set, lty=1)
title("Smooth Spline (n=100)")


## 4.5 Random Fourier Features

### Example 68

sigma=10; sigma2=sigma^2
k=function(x,y) exp(-(x-y)^2/2/sigma2)
z=function(x) sqrt(2/m)*cos(w*x+b)
zz=function(x,y) sum(z(x)*z(y))
u=matrix(0,1000,3)
m_seq=c(20,100,400)
for(i in 1:1000){
  x=rnorm(1); y=rnorm(1)
  for(j in 1:3){
    m=m_seq[j]; w=rnorm(m)/sigma; b=runif(m)*2*pi
    u[i,j]=zz(x,y)-k(x,y)
  }
}
boxplot(u[,1],u[,2],u[,3], ylab="Difference from k(x,y)",names=paste0("m=",m_seq),
        col=c("red","blue","green"), main="Kernel approximation with RFF")


### Example 69

sigma=10; sigma2=sigma^2

## Funtion z
m=20; w=rnorm(m)/sigma; b=runif(m)*2*pi
z=function(u,m) sqrt(2/m)*cos(w*u+b)

## Gauss kernel
k.g=function(x,y)exp(-(x-y)^2/2/sigma2)

## Data Generation
n=200; x=rnorm(n)/2; y=1+5*sin(x/10)+5*x^2+rnorm(n)  
x.min=min(x); x.max=max(x); y.min=min(y); y.max=max(y)
lambda=0.001  ## lambda=0.9 also

# Function of low-rank approximation
alpha.rff=function(x,y,m){
  n=length(x)
  Z=array(dim=c(n,m))
  for(i in 1:n)Z[i,]=z(x[i],m)
  beta=solve(t(Z)%*%Z+lambda*diag(m))%*%t(Z)%*%y
  return(as.vector(beta))
}
# Normal Function
alpha=function(k,x,y){
  n=length(x); K=matrix(0,n,n); for(i in 1:n)for(j in 1:n)K[i,j]=k(x[i],x[j])
  alpha=solve(K+lambda*diag(n))%*%y
  return(as.vector(alpha))
}
# Numerical Comparison
alpha.hat=alpha(k.g,x,y)
beta.hat=alpha.rff(x,y,m)
r=sort(x); u=array(n); v=array(n)
for(j in 1:n){
  S=0;for(i in 1:n)S=S+alpha.hat[i]*k.g(x[i],r[j]); u[j]=S
  v[j]=sum(beta.hat*z(r[j],m))
}
plot(r,u,type="l",xlim=c(x.min,x.max),ylim=c(y.min,y.max),xlab="x", ylab="y",col="red",
     main="lambda=10^{-4},m=20,n=200")
lines(r,v,col="blue"); points(x,y)
legend("topleft",lwd=1,c("no approximation","approximation"), col=c("red","blue"))


## 4.6 Nystrom Approximation

### Example 70

sigma2=1; k.g=function(x,y)exp(-(x-y)^2/2/sigma2)
n=300; x=rnorm(n)/2; y=3-2*x^2+3*x^3+2*rnorm(n)       ## Data Generation
lambda=10^(-5)  ## lambda=0.9 also
m=10
# Function of low-rank approximation
alpha.m=function(k,x,y,m){
  n=length(x); K=matrix(0,n,n); for(i in 1:n)for(j in 1:n)K[i,j]=k(x[i],x[j])
  A=svd(K[1:m,1:m])
  u=array(dim=c(n,m)); 
  for(i in 1:m)for(j in 1:n)u[j,i]=sqrt(m/n)*sum(K[j,1:m]*A$u[1:m,i])/A$d[i]
  mu=A$d*n/m
  R=sqrt(mu[1])*u[,1]; for(i in 2:m)R=cbind(R,sqrt(mu[i])*u[,i])
  alpha=(diag(n)-R%*%solve(t(R)%*%R+lambda*diag(m))%*%t(R))%*%y/lambda
  return(as.vector(alpha))
}
# Normal Function
alpha=function(k,x,y){
  n=length(x); K=matrix(0,n,n); for(i in 1:n)for(j in 1:n)K[i,j]=k(x[i],x[j])
  alpha=solve(K+lambda*diag(n))%*%y
  return(as.vector(alpha))
}
# Numerical Comparison
alpha.1=alpha(k.g,x,y); alpha.2=alpha.m(k.g,x,y,m)
z=sort(x); u=array(n); v=array(n)
for(j in 1:n){
  S=0;for(i in 1:n)S=S+alpha.1[i]*k.g(x[i],z[j]); u[j]=S
  S=0;for(i in 1:n)S=S+alpha.2[i]*k.g(x[i],z[j]); v[j]=S
}
plot(z,u,type="l",xlim=c(-1,1),xlab="x", ylab="y", ylim=c(-1,5),col="red",main="Kernel Ridge")
lines(z,v,col="blue"); points(x,y)
legend("topleft",lwd=1,c("no approximation","approximation"), col=c("red","blue"))


## 4.7 Incomplete Cholesky decomposition

im.ch=function(A,m=ncol(A)){
    n=ncol(A); R=matrix(0,n,n); P=diag(n)
	for(i in 1:n)R[i,i]=sqrt(A[i,i])
	max.R=0;for(i in 1:n)if(R[i,i]>max.R){k=i; max.R=R[i,i]}
	R[1,1]=max.R
	if(k != 1){
        w=A[,k]; A[,k]=A[,1]; A[,1]=w
		w=A[k,]; A[k,]=A[1,]; A[1,]=w
		P[1,1]=0; P[k,k]=0; P[1,k]=1; P[k,1]=1
	}
    for(i in 2:n)R[i,1]=A[i,1]/R[1,1]
	if(m>1)for(i in 2:m){
	    for(j in i:n)R[j,j]=sqrt(A[j,j]-sum(R[j,1:(i-1)]^2))
		max.R=0;for(j in i:n)if(R[j,j]>max.R){k=j; max.R=R[j,j]}
		R[i,i]=max.R
		if(k!=i){
		    w=R[i,1:(i-1)]; R[i,1:(i-1)]=R[k,1:(i-1)]; R[k,1:(i-1)]=w
			w=A[,k]; A[,k]=A[,i]; A[,i]=w
			w=A[k,]; A[k,]=A[i,]; A[i,]=w
                        Q=diag(n); Q[i,i]=0; Q[k,k]=0; Q[i,k]=1; Q[k,i]=1; P=P%*%Q
		}
		if(i<n)for(j in (i+1):n)R[j,i]=(A[j,i]-sum(R[i,1:(i-1)]*R[j,1:(i-1)]))/R[i,i]
	}
	if(m<n)for(i in (m+1):n)R[,i]=0
	return(list(P=P,R=R))
}

# Data Generation
n=4; range=-5:5
D=matrix(sample(range,n*n,replace=TRUE),n,n)
A=t(D)%*%D

# Execution example
im.ch(A)
L=im.ch(A)$R; L%*%t(L)
P=im.ch(A)$P; t(P)%*%A%*%P
P%*%(L%*%t(L))%*%t(P)  ## Match A

im.ch(A,2)
L=im.ch(A,2)$R; L%*%t(L)
P=im.ch(A,2)$P; t(P)%*%A%*%P
P%*%(L%*%t(L))%*%t(P)   ## Low rank approximation of A
