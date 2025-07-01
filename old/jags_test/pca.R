# PCA with JAGS
library(ggplot2)
library(rjags)

norm <- function(v) {
  sqrt(sum(v^2))
}

colNorms <- function(x) {
  #compute the L2 norms of columns of a matrix
  apply(x, 2, norm)
}

ortho <- function(u, v, ret = c("m", "df")) {
  #convert factors to orthonormal basis
  #u is LxN, v is LxG
  #returns NxL data frame of factors (similar to PCs) suitable for plotting
  #also returns the GxL loading matrix
  #loading matrix describes how important feature g is in dimension l
  #note we assume columns are samples, rows are features
  #this is the transpose of what prcomp expects
  #eg if we did prcomp(t(Y)) would get t(Y) = "x"%*%t("rotation")
  #the orientation of the returned loadings matches "rotation" from prcomp
  #the orientation of the returned factors matches "x" from prcomp
  #default is to return matrices, can also return data.frames if specified
  ret <- match.arg(ret)
  L <- nrow(u)
  #step 2, convert loadings to orthonormal basis
  svd_v <- svd(v)
  A <- svd_v$u #LxL
  D <- if (length(svd_v$d) > 1) diag(svd_v$d) else svd_v$d #LxL
  loadings <- svd_v$v #GxL with orthonormal rows
  factors <- crossprod(u, A %*% D) #NxL
  #step 3, reorder dimensions in decreasing magnitude
  o <- order(colNorms(factors), decreasing = TRUE)
  factors <- factors[, o]
  loadings <- loadings[, o]
  colnames(loadings) <- colnames(factors) <- paste0("dim", 1:L)
  if (ret == "df") {
    loadings <- as.data.frame(loadings)
    factors <- as.data.frame(factors)
  }
  mget(c("factors", "loadings"))
}


# Iris Data
#standard PCA
Y <- scale(iris[, 1:4])
pca_factors <- as.data.frame(prcomp(Y, rank = 2)$x)
colnames(pca_factors) <- paste0("dim", 1:2)
pd1 <- cbind(pca_factors, sp = iris$Species, alg = "pca")
ggplot(pd1, aes(x = dim1, y = dim2, colour = sp)) + geom_point()

#7.59 sec for nested loops (L=3 --> 11.6 sec)
#12.46 sec for multivariate normal prior (L=3 --> 17.6 sec)

jmod <- jags.model(
  "pca.bugs",
  data = list(y = Y, n = nrow(Y), d = ncol(Y), alpha = 1e-4, L = 2),
  n.chains = 3,
  quiet = TRUE
)
update(jmod, 1000, progress.bar = "none")
system.time(
  samps0 <- coda.samples(
    jmod,
    variable.names = c("mu", "u", "v", "tau_y"),
    n.iter = 5000,
    progress.bar = "none"
  )
)
#list of length n.chains
#convert to data frame
samps <- do.call("rbind", samps0)
#mu should be close to zero because already centered the data
mu <- colMeans(samps[, grepl("mu", colnames(samps))])
tau_y <- samps[, grepl("tau_y", colnames(samps))]
#only use last sample for u,v to avoid rotational invariance problem
u <- samps[nrow(samps), grepl("^u\\[", colnames(samps))]
U <- t(matrix(u, ncol = 2)) #LxN
v <- samps[nrow(samps), grepl("^v\\[", colnames(samps))]
V <- t(matrix(v, ncol = 2)) #LxG
#check prediction accuracy
Yhat <- crossprod(U, V)
(rmse_bpca <- sqrt(mean((Y - Yhat)^2))) #should be small
#orthogonalize the loadings
res <- ortho(U, V, ret = "df")
pd2 <- cbind(res$factors, sp = iris$Species, alg = "jags")
pd <- rbind(pd1, pd2)
ggplot(pd, aes(x = dim1, y = dim2, colour = sp)) +
  geom_point() +
  facet_wrap(~alg)
