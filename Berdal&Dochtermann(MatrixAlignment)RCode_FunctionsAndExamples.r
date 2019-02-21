#Supplemental materials: Adaptive alignment of plasticity with genetic variation and selection

#Monica Anderson Berdal a, b 
#Ned A. Dochtermann a, c
#a Department of Biological Sciences; North Dakota State University
#b m.anderson.berdal@gmail.com 
#c corresponding author: ned.dochtermann@gmail.com

#Running title: Alignment of plasticity, G, and selection 


#Supplemental Text 1. R functions and code for the worked example under “Expectations for the alignment of E”

#### Functions ####
#Extract dominant eigenvector
vec.out<-function(x,traits){y<-matrix(x,traits)
eigen(y)$vectors[,1]}

#Equation 1
#Function for calculation of vector correlation
vec.cor<-function(z1=z1,z2=z2){
  abs(sum(z1*z2) / ( sqrt(sum(z1 * z1)) * sqrt(sum(z2 * z2)) ) )
}

#Equation 2
#Converts vector correlation to degrees
vec.angle<-function(x){acos(x)*(180/pi)

#Equation 3 and corresponding significance test
#Calculate significance of a vector correlation
cor.sig<-function(r,r.null=0.975,n){
  Z=(atanh(r)-atanh(r.null))/(sqrt(2/(n-3)))
  p <- 2*pnorm(abs(Z), lower=F) 
  return(list(Z=Z,p=p))
}

#Equation 4
#Estimates (or calculates) the posterior distribution of psi.theta
Ovaskainen.etal2008<-function(A,B,k,samp=1000){
  #function fails if entire posterior is imported, the solution to this was to delete the first posterior sample
  A<-A[-1,]
  B<-B[-1,]
  
  if (dim(A)[2] != dim(B)[2]) {
    stop("matrices/vectors are of different sizes")
  }
  if (dim(A)[1] != dim(B)[1]) {
    stop("posterior distributions are of different lengths")
  }
  if (dim(A)[2] != k^2) {
    stop("number of traits indicated does not match matrix/vector size")
  }
  if (dim(A)[2] != k^2) {
    stop("number of traits indicated does not match matrix/vector size")
  }
  
  post.length<-dim(A)[1]
  if(samp=="All"){i<-((post.length*(post.length-1))/2)
  }  else{i<-samp}
  
  combs<-combinations(post.length,2,repeats.allowed=F)
  if(samp=="All"){combs<-combs
  } else{combs<-combs[sample(dim(combs)[1],i,FALSE),]}
  
  stor<-matrix(NA,i)
  for(j in 1:i){
    domA.1<-vec.out(A[combs[j,1],],k)
    domB.1<-vec.out(B[combs[j,1],],k)
    domA.2<-vec.out(A[combs[j,2],],k)
    domB.2<-vec.out(B[combs[j,2],],k)
    
    cor.A1A2<-vec.cor(domA.1,domA.2)
    cor.B1B2<-vec.cor(domB.1,domB.2)
    cor.A1B2<-vec.cor(domA.1,domB.2)
    cor.A2B1<-vec.cor(domB.1,domA.2)
    
    stor[j]<-(cor.A1A2+cor.B1B2)-(cor.A1B2+cor.A2B1)
  }
  as.mcmc(stor)
}



#Converts correlation matrices to covariance matrices 
#given known variances (vars) and a correlation matrix (cormat)
cor2cov<-function(vars,cormat){
  sdMat<-diag(sqrt(vars))
  corMat<-cormat
  mat<-sdMat %*% corMat %*% t(sdMat)
  return(mat)
}



#runif(3,.1,10)
P.vars<-c(0.60, 4.27, 4.13)
G.var<-.3*P.vars
E.var<-0.7*P.vars

G.cor<-matrix(c(1,.57,.57,
            .57, 1, .57,
            .57,.57,1),3)

E.cor1<-diag(3) #first scenario of worked example

E.cor2<-matrix(c(1,.19,.19,
                 .19, 1, .19,
                 .19,.19,1),3)


E.cor3<-matrix(c(1,0.57,0.57,
                 0.57, 1, 0.57,
                 0.57,0.57,1),3)

#second scenario of worked example (based on equations 7,8, and 9)

G.cov<-cor2cov(G.var,G.cor)
E.cov1<-cor2cov(E.var,E.cor1)
E.cov2<-cor2cov(E.var,E.cor2)
E.cov3<-cor2cov(E.var,E.cor3)

gmax<-vec.out(G.cov,3)
emax1<-vec.out(E.cov1,3)
emax2<-vec.out(E.cov2,3)

vec.cor(gmax,emax1)
vec.angle(vec.cor(gmax,emax1))

vec.cor(gmax,emax2)
vec.angle(vec.cor(gmax,emax2))

