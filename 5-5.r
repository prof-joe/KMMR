# Ch.5　MMD and HSIC

## 5.1 Random Variables in RKHS

## 5.2 MMD and 2-sample problem

## Example 71

sigma=1; k=function(x,y)exp(-(x-y)^2/sigma^2)
## Data Generation
n=100
xx=rnorm(n)
yy=rnorm(n)       ## When the distributions are equal
# yy=rnorm(n)*2   ## When the distributions are not equal
x=xx;y=yy 
## Calculation of the null distribution
T=NULL
for(h in 1:100){
  index1=sample(n,n/2)
  index2=setdiff(1:n,index1)
  x=c(xx[index2],yy[index1])
  y=c(xx[index1],yy[index2])
  S=0
  for(i in 1:n)for(j in 1:n)if(i!=j)S=S+k(x[i],x[j])+k(y[i],y[j])-k(x[i],y[j])-k(x[j],y[i])
  T=c(S/n/(n-1),T)
}
v=quantile(T,0.95)
## Calculation of statistics
S=0; for(i in 1:n)for(j in 1:n)if(i!=j)S=S+k(x[i],x[j])+k(y[i],y[j])-k(x[i],y[j])-k(x[j],y[i])
u=S/n/(n-1)
## Illustration of graphs
plot(density(T),xlim=c(min(T,v,u),max(T,v,u)))
abline(v=v,col="red",lty=2,lwd=2)
abline(v=u,col="blue",lty=1,lwd=2)


### Example 73

sigma=1; k=function(x,y)exp(-(x-y)^2/sigma^2)
## Data Generation
n=100
x=rnorm(n)
y=rnorm(n)        ## When the distributions are equal
# y=rnorm(n)*2    ## When the distributions are not equal
## Calculation of the null distribution
K=matrix(0,n,n)
for(i in 1:n)for(j in 1:n)K[i,j]=k(x[i],x[j])+k(y[i],y[j])-k(x[i],y[j])-k(x[j],y[i])
lambda=eigen(K)$values/n
r=20
z=NULL
for(h in 1:10000)z=c(z,1/n*(sum(lambda[1:r]*(rchisq(1:r, df=1)-1))))
v=quantile(z,0.95)
## Calculation of statistics
S=0
for(i in 1:(n-1))for(j in (i+1):n)S=S+k(x[i],x[j])+k(y[i],y[j])-k(x[i],y[j])-k(x[j],y[i])
u=S/n/(n-1)
## Illustration of graphs
plot(density(z),xlim=c(min(z,v,u),max(z,v,u)))
abline(v=v,col="red",lty=2,lwd=2)
abline(v=u,col="blue",lty=1,lwd=2)

## 5.3 HSIC and Independence Tests

HSIC.1=function(x,y,k.x,k.y){
     n=length(x)
     S=0;for(i in 1:n)for(j in 1:n)S=S+k.x(x[i],x[j])*k.y(y[i],y[j])
     T=0;
     for(i in 1:n){
          T.1=0; for(j in 1:n)T.1=T.1+k.x(x[i],x[j]);
          T.2=0; for(l in 1:n)T.2=T.2+k.y(y[i],y[l]);
          T=T+T.1*T.2
          }
     U=0; for(i in 1:n)for(j in 1:n)U=U+k.x(x[i],x[j])
     V=0; for(i in 1:n)for(j in 1:n)V=V+k.y(y[i],y[j])
     return(S/n^2-2*T/n^3+U*V/n^4)
}

HSIC.1=function(x,y,k.x,k.y){
     n=length(x)
     K.x=matrix(0,n,n)
     for(i in 1:n)for(j in 1:n)K.x[i,j]=k.x(x[i],x[j])
     K.y=matrix(0,n,n)
     for(i in 1:n)for(j in 1:n)K.y[i,j]=k.y(y[i],y[j])
     E=matrix(1,n,n)
     H=diag(n)-E/n
     return(sum(diag(K.x%*%H%*%K.y%*%H))/n^2)
}

### Example 76

k.x=function(x,y)exp(-norm(x-y,"2")^2/2); k.y=k.x
n=100
for(a in c(0,0.1,0.2,0.4,0.6,0.8)){      ## a is the correlation coefficient
     x=rnorm(n); z=rnorm(n); y=a*x+sqrt(1-a^2)*z
     print(HSIC.1(x,y,k.x,k.y))
}

HSIC.2=function(x,y,z,k.x,k.y,k.z){
     n=length(x)
     S=0;for(i in 1:n)for(j in 1:n)S=S+k.x(x[i],x[j])*k.y(y[i],y[j])*k.z(z[i],z[j])
     T=0;
     for(i in 1:n){
          T.1=0; for(j in 1:n)T.1=T.1+k.x(x[i],x[j]);
          T.2=0; for(l in 1:n)T.2=T.2+k.y(y[i],y[l])*k.z(z[i],z[l]);
          T=T+T.1*T.2
          }
     U=0; for(i in 1:n)for(j in 1:n)U=U+k.x(x[i],x[j])
     V=0; for(i in 1:n)for(j in 1:n)V=V+k.y(y[i],y[j])*k.z(z[i],z[j])
     return(S/n^2-2*T/n^3+U*V/n^4)
}

### Example 77

cc=function(x,y)sum(x*y)/length(x)       ## Partial correlation coefficient at cc(x,y)/cc(x,x)
f=function(u,v)u-cc(u,v)/cc(v,v)*v       ## residual

## Data Generation ##
n=30
x=rnorm(n)^2-rnorm(n)^2; y=2*x+rnorm(n)^2-rnorm(n)^2; z=x+y+rnorm(n)^2-rnorm(n)^2
x=x-mean(x); y=y-mean(y); z=z-mean(z)
## Estimate the top ##
cc=function(x,y)sum(x*y)/length(x)
f=function(u,v)u-cc(u,v)/cc(v,v)*v
x.y=f(x,y); y.z=f(y,z); z.x=f(z,x); x.z=f(x,z); z.y=f(z,y); y.x=f(y,x)
v1=HSIC.2(x,y.x,z.x,k.x,k.y,k.z); v2=HSIC.2(y,z.y,x.y,k.y,k.z,k.x)
     v3=HSIC.2(z,x.z,y.z,k.z,k.x,k.y)
if(v1<v2){if(v1<v3)top=1 else top=3} else {if(v2<v3)top=2 else top=3} ##

## Estimate the bottom ##
x.yz=f(x.y,z.y); y.zx=f(y.z,x.z); z.xy=f(z.x,y.x)
if(top==1){
     v1=HSIC.1(y.x,z.xy,k.y,k.z); v2=HSIC.1(z.x,y.zx,k.z,k.y)
     if(v1<v2){middle=2; bottom=3} else {middle=3; bottom=2}
}
if(top==2){
     v1=HSIC.1(z.y,x.yz,k.z,k.x); v2=HSIC.1(x.y,z.xy,k.x,k.z)
     if(v1<v2){middle=3; bottom=1} else {middle=1; bottom=3}
}
if(top==3){
     v1=HSIC.1(z.y,x.yz,k.z,k.x); v2=HSIC.1(x.y,z.xy,k.x,k.z)
     if(v1<v2){middle=1; bottom=2} else {middle=2; bottom=1}
}
## Output results ##
print(paste("top=",top))
print(paste("middle=",middle))
print(paste("bottom=",bottom))


### Example 78

## Sort x to display a histogram of HSIC distribution ##
## Data Generation ##
x=rnorm(n); y=rnorm(n); u=HSIC.1(x,y,k.x,k.y)
## Sort x to construct the null distribution
m=100; w=NULL; 
for(i in 1:m){x=x[sample(n,n)]; w=c(w,HSIC.1(x,y,k.x,k.y))}
## Set rejection area
v=quantile(w,0.95)
## Display in graph
plot(density(w),xlim=c(min(w,v,u),max(w,v,u)))
abline(v=v,col="red",lty=2,lwd=2)
abline(v=u,col="blue",lty=1,lwd=2)

h=function(i,j,q,r,x,y,k.x,k.y){ 
  M=combn(c(i,j,q,r),m=4)
  m=ncol(M)
  S=0
  for(j in 1:m){
    t=M[1,j]; u=M[2,j]; v=M[3,j]; w=M[4,j]
    S=S+k.x(x[t],x[u])*k.y(y[t],y[u])+k.x(x[t],x[u])*k.y(y[v],y[w])-2*k.x(x[t],x[u])*k.y(y[t],y[v])
  }
  return(S/m)
}
HSIC.U=function(x,y,k.x,k.y){
M=combn(1:n,m=4)
m=ncol(M)
S=0
for(j in 1:m)S=S+h(M[1,j],M[2,j],M[3,j],M[4,j],x,y,k.x,k.y)
return(S/choose(n,4))
}


### Example 79

sigma=1; k=function(x,y)exp(-(x-y)^2/sigma^2); k.x=k; k.y=k
## Data Generation
n=100; x=rnorm(n)
a=0        ## When it is independent
#a=0.2       ## Correlation coefficient 0.2
y=a*x+sqrt(1-a**2)*rnorm(n)
# y=rnorm(n)*2    ## When the distributions are not equal
## Calculation of the null distribution
K.x=matrix(0,n,n); for(i in 1:n)for(j in 1:n)K.x[i,j]=k.x(x[i],x[j])
K.y=matrix(0,n,n); for(i in 1:n)for(j in 1:n)K.y[i,j]=k.y(y[i],y[j])
F=array(0,dim=n); for(i in 1:n)F[i]=sum(K.x[i,])/n
G=array(0,dim=n); for(i in 1:n)G[i]=sum(K.y[i,])/n
H=sum(F)/n
I=sum(G)/n
K=matrix(0,n,n)
for(i in 1:n)for(j in 1:n)K[i,j]=(K.x[i,j]-F[i]-F[j]+H)*(K.y[i,j]-G[i]-G[j]+I)/6
r=20
lambda=eigen(K)$values/n
z=NULL; for(s in 1:10000)z=c(z,1/n*(sum(lambda[1:r]*(rchisq(1:r, df=1)-1))))
v=quantile(z,0.95)
## Calculation of statistics
u=HSIC.U(x,y,k.x,k.y)
## Illustration of graphs
plot(density(z),xlim=c(min(z,v,u),max(z,v,u)))
abline(v=v,col="red",lty=2,lwd=2)
abline(v=u,col="blue",lty=1,lwd=2)


## 5.4 Characteristic Kernels and Universal Kernels

## 5.5 Introduction to Empirical Processes
