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
#step 1 compile the jags model from the specification file
jmod <- jags.model("spike_slab_regression.bugs", data = list(y=y,N=N,X=X,k=k), n.chains=3, quiet=TRUE)
#step 2 initialize all variables and run for the burn-in period
update(jmod, 1000, progress.bar="none")
#step 3 sample from the model and store in a convenient coda object
samps0<-coda.samples(jmod, variable.names=c("beta","delta"), n.iter=5000, progress.bar="none") #list of length n.chains
traceplot(samps0)

samps<-do.call("rbind",samps0)
beta   <- samps[,1:k]
delta  <- samps[,(k+1):(2*k)]
colnames(beta)<-colnames(delta)<- colnames(X)
boxplot(beta,outline=FALSE)
colMeans(delta)
hist(beta[,1])
hist(beta[,4],breaks=100)
hist(beta[,6],breaks=100)
