library(rjags)
library(ggplot2)
library(MASS)

n <- 500
mu1 <- c(10, 10)
mu2 <- c(-10, 0)
mu3 <- c(0, -10)
mu_mat <- cbind(mu1, mu2, mu3)
mu_grand <- c(20, -50)
probs_true <- probs <- c(.5, .3, .2)
z_true <- z <- rmultinom(n, 1, probs)
#each row of ymn a data point, each column a dimension
ymn <- t(mu_grand + mu_mat %*% z)
Sigma <- matrix(c(5, 4, 4, 5), nrow = 2)
y <- t(sapply(1:n, function(n) {
  mvrnorm(1, ymn[n, ], Sigma)
}))
colnames(y) <- paste0("dim", 1:2)
(baseplot <- ggplot(as.data.frame(y), aes(x = dim1, y = dim2)) +
  geom_point() +
  geom_density2d())

jmod <- jags.model(
  "finite_mixture.bugs",
  data = list(y = y, n = nrow(y), d = ncol(y), Nclust = 3, Omega0 = diag(2)),
  n.chains = 3,
  quiet = TRUE
)
update(jmod, 1000, progress.bar = "none")
system.time(
  samps0 <- coda.samples(
    jmod,
    variable.names = c("mu_cl", "z", "probs", "Omega"),
    n.iter = 5000,
    progress.bar = "none"
  )
)
#list of length n.chains
#convert to data frame
samps <- do.call("rbind", samps0)
mu_cl <- samps[, grepl("mu_cl", colnames(samps))]
z <- samps[, grepl("z", colnames(samps))]
probs <- samps[, grepl("probs", colnames(samps))]
Omega <- samps[, grepl("Omega", colnames(samps))]

#check results invariant to label switching
#covariance matrix
solve(matrix(colMeans(Omega), nrow = 2)) #5,4; 4,5
#check results that are vulnerable to label switching
#use a single sample from very end of chain
tail(probs, 1) #c(.5,.3,.2) order doesn't matter
tail(mu_cl, 1)
s <- sample.int(nrow(mu_cl), 100)
f <- function(k) {
  m <- as.data.frame(mu_cl[s, c(k, k + 3)])
  colnames(m) <- c("dim1", "dim2")
  m$k <- paste0("clust", k)
  m
}
mu_s <- do.call("rbind", lapply(1:3, f))
baseplot +
  geom_point(
    aes(x = dim1, y = dim2),
    data = mu_s,
    colour = "red",
    shape = 2,
    size = 2
  )
