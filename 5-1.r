# Ch.1 Positive Definite Kernel

## 1.1 Positive Definiteness of Matrices

### Example 1
n=3
B=matrix(rnorm(n^2),3,3)
A=t(B)%*%B
eigen(A)

S=NULL
for(i in 1:10){
  z=rnorm(n)
  y=drop(t(z)%*%A%*%z)
  S=c(S,y)
}
print(S)


## 1.2 Kernel

### Example 2

n=250; x=2*rnorm(n); y=sin(2*pi*x)+rnorm(n)/4  ## Data Generation
D=function(t) max(0.75*(1-t^2),0)              ## definition of function  D
k=function(x,y,lambda) D(abs(x-y)/lambda)      ## definition of function  K
f=function(z,lambda){                          ## definition of function  f
  S=0; T=0;
  for(i in 1:n){S=S+k(x[i],z,lambda)*y[i]; T=T+k(x[i],z,lambda)}
  return(S/T)
}
plot(seq(-3,3,length=10),seq(-2,3,length=10),type="n",xlab="x", ylab="y"); points(x,y)
xx=seq(-3,3,0.1)
yy=NULL;for(zz in xx)yy=c(yy,f(zz,0.05)); lines(xx,yy,col="green")
yy=NULL;for(zz in xx)yy=c(yy,f(zz,0.35)); lines(xx,yy,col="blue")
yy=NULL;for(zz in xx)yy=c(yy,f(zz,0.50)); lines(xx,yy,col="red")
title("Nadaraya-Watson Estimate")
legend("topleft",legend=paste0("lambda=",c(0.05, 0.35, 0.50)),
       lwd=1,col=c("green","blue","red"))

## 1.3 Positive Definite Kernel

### Example 12

K=function(x,y,sigma2)exp(-norm(x-y,"2")^2/2/sigma2)
f=function(z,sigma2){ ## definition of function  f
  S=0; T=0;
  for(i in 1:n){S=S+K(x[i],z,sigma2)*y[i]; T=T+K(x[i],z,sigma2)}
  return(S/T)
}

n=100; x=2*rnorm(n); y=sin(2*pi*x)+rnorm(n)/4 ## Data generation
plot(seq(-3,3,length=10),seq(-2,3,length=10),type="n",xlab="x", ylab="y"); points(x,y)
xx=seq(-3,3,0.1)
yy=NULL;for(zz in xx)yy=c(yy,f(zz,0.001)); lines(xx,yy,col="green")
yy=NULL;for(zz in xx)yy=c(yy,f(zz,0.01)); lines(xx,yy,col="blue")
## So far，the curves for sigma2=0.01 and 0.001 are plotted
m=n/10
sigma2.seq=seq(0.001,0.01,0.001); SS.min=Inf
for(sigma2 in sigma2.seq){
     SS=0
     for(k in 1:10){
          test=((k-1)*m+1):(k*m); train=setdiff(1:n,test)
          for(j in test){
               u=0; v=0;
               for(i in train){
                    kk=K(x[i],x[j],sigma2); u=u+kk*y[i]; v=v+kk ## Applying the Kernel
               }
               if(v!=0)z=u/v; SS=SS+(y[j]-z)^2
          }
     }
     if(SS<SS.min){SS.min=SS;sigma2.best=sigma2}
}
paste0("Best sigma2 = ",sigma2.best)
## So far, the optimal lambda value is computed.
yy=NULL;for(zz in xx)yy=c(yy,f(zz,sigma2.best)); lines(xx,yy,col="red")
title("Nadaraya-Watson Estimate")
legend("topleft",legend=paste0("sigma2=",c(0.01, 0.001, "sigma2.best")),
lwd=1,col=c("green","blue","red"))


## 1.4 Probability

## 1.5 Bochner's theorem

## 1.6 String,Tree,and Graph Kernels

### Example 21

string.kernel=function(x,y){
  m=nchar(x)
  n=nchar(y)
  S=0
  for(i in 1:m)for(j in i:m)for(k in 1:n)
    if(substring(x,i,j)==substring(y,k,k+j-i))S=S+1
  return(S)
}

C=c("a","b","c")
m=10; w=sample(C,m,rep=TRUE)
x=NULL; for(i in 1:m)x=paste0(x,w[i])
n=12; w=sample(C,n,rep=TRUE)
y=NULL; for(i in 1:m)y=paste0(y,w[i])

x
y
string.kernel(x,y)

### Example 22

C=function(i,j){
  S=s[[i]]; T=t[[j]]
## returns 0 if the labels of vertex i in tree s or vertex j in tree t do not match.
  if(S[[1]]!=T[[1]])return(0)
## Return 0 if vertex i of tree s or vertex j of tree t has no descendants.
  if(is.null(S[[2]]))return(0)
  if(is.null(T[[2]]))return(0)
  if(length(S[[2]])!=length(T[[2]]))return(0)
  U=NULL; for(x in S[[2]])U=c(U,s[[x]][[1]]); U1=sort(U)
  V=NULL; for(y in T[[2]])V=c(V,t[[y]][[1]]); V1=sort(V)
  m=length(U)
## If the labels of the descendants do not match, return 0.
  for(h in 1:m)if(U1[h]!=V1[h])return(0)
  U2=S[[2]][order(U)]
  V2=T[[2]][order(V)]
  W=1; for(h in 1:m)W=W*(1+C(U2[h],V2[h]))
  return(W)
}
k=function(s,t){
  m=length(s); n=length(t)
  kernel=0
  for(i in 1:m)for(j in 1:n)if(C(i,j)>0)kernel=kernel+C(i,j)
  return(kernel)
}

## Trees are described by lists. Labels and their descendants (displayed as vectors)
s=list()
s[[1]]=list("G",c(2,4)); s[[2]]=list("T",3);    s[[3]]=list("C",NULL)
s[[4]]=list("A",c(5,6)); s[[5]]=list("C",NULL); s[[6]]=list("T",NULL)
t=list()
t[[1]]=list("G",c(2,5)); t[[2]]=list("A",c(3,4)); t[[3]]=list("C",NULL)
t[[4]]=list("T",NULL); t[[5]]=list("T",c(6,7)); t[[6]]=list("C",NULL)
t[[7]]=list("A",c(8,9)); t[[8]]=list("C",NULL); t[[9]]=list("T",NULL)

for(i in 1:6)for(j in 1:9)if(C(i,j)>0)print(c(i,j,C(i,j)))
k(s,t)


### Example 23

k=function(s,p) prob(s,p)/length(node)
prob=function(s,p){
  if(length(node[s[1]])==0) return(0)
  if(length(s)==1) return (p)
  m=length(s)
  S=(1-p)/length(node[s[1]])*prob(s[2:m],p)
  return(S)
}

node=list()
node[[1]]=c(2,4); node[[2]]=4; node[[3]]=c(1,5); node[[4]]=3; node[[5]]=3

k(c(1,4,3,5,3),1/3)

