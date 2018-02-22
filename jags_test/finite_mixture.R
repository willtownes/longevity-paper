library(rjags)
library(ggplot2)
library(MASS)

n<-500
mu1<-c(10,10)
mu2<-c(-10,0)
mu3<-c(0,-10)
mu_mat<-cbind(mu1,mu2,mu3)
mu_grand<-c(20,-50)
probs_true<-probs<-c(.5,.3,.2)
z_true<-z<-rmultinom(n,1,probs)
#each row of ymn a data point, each column a dimension
ymn<-t(mu_grand+mu_mat%*%z) 
Sigma<-matrix(c(5,4,4,5),nrow=2)
y<-t(sapply(1:n,function(n){mvrnorm(1,ymn[n,],Sigma)}))
colnames(y)<-paste0("dim",1:2)
ggplot(as.data.frame(y),aes(x=dim1,y=dim2))+geom_point()+geom_density2d()

model_string<- "model {
  # Likelihood:
  for( i in 1 : n ) {
    y[i,1:d] ~ dmnorm(mu[i,1:d], Omega)
    mu[i,1:d] <- mu_cl[z[i], 1:d]
    z[i] ~ dcat(probs)
  }
  # Prior:
  Omega ~ dwish(Omega0, 3) #Sigma^{-1} prior
  for ( k in 1: Nclust ) {
    mu_cl[k,1:d] ~ dmnorm( c(0,0) , .0001*Omega0 )
  }
  probs ~ ddirch(c(1,1,1))
}"

jmod <- jags.model(textConnection(model_string), data = list(y=y,n=nrow(y),d=ncol(y),Nclust=3,Omega0=diag(2)), n.chains=3, quiet=TRUE)
update(jmod, 100, progress.bar="none")
system.time(SAMPS <- coda.samples(jmod, variable.names=c("mu_cl","z","probs","Omega"), n.iter=5000, progress.bar="none")) #list of length n.chains

samps<-do.call("rbind",SAMPS)
mu_cl<-samps[,grepl("mu_cl",colnames(samps))]
z<-samps[,grepl("z",colnames(samps))]
probs<-samps[,grepl("probs",colnames(samps))]
Omega<-samps[,grepl("Omega",colnames(samps))]
