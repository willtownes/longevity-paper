# adapted from https://www4.stat.ncsu.edu/~reich/ABA/code/SSVS
library(rjags)

beta_true<-c(10,-5,.5,0.001,-0.001,0) #last three are essentially zero
mu_true<-30
sy<-3
N<-500
k<-length(beta_true)
X<-matrix(rnorm(N*k),nrow=N)
y<-drop(mu_true+X%*%beta_true+rnorm(N,sd=sy))
summary(lm(y~X))

#jags part
model_string <- "model{
  # Likelihood
  for(i in 1:N){
    #variance=1/tau_y
    y[i] ~ dnorm(eta[i],tau_y) 
    #eta[i] <- mu + inprod(beta[],X[i,])
  }
  eta<-mu+X%*%beta

  #Priors
  for(j in 1:k){
    gamma[j] ~ dnorm(0,tau_gam)
    delta[j] ~ dbern(omega)
    beta[j]  <- gamma[j]*delta[j]
  }

  omega ~ dunif(0,1)
  tau_gam ~ dgamma(.1,.1)
  tau_y ~ dgamma(.001,.001)
  mu ~ dt(0,0.001,4)
}"

SSVS <- jags.model(textConnection(model_string), data = list(y=y,N=N,X=X,k=k), n.chains=3, quiet=TRUE)
update(SSVS, 10000, progress.bar="none")
SAMPS <- coda.samples(SSVS, variable.names=c("beta","delta"), n.iter=50000, progress.bar="none") #list of length n.chains

samps<-do.call("rbind",SAMPS)
beta   <- samps[,1:k]
delta  <- samps[,1:k+k]
colnames(beta)<-colnames(delta)<- colnames(X)
boxplot(beta,outline=FALSE)
colMeans(delta)
hist(beta[,1])
hist(beta[,4],breaks=100)
hist(beta[,6],breaks=100)
