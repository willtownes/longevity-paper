//Full model for integrating gene expression with survival outcome
#ncell_s integer vector of length Ns, number of cells for each survival pert
#ncell_x integer, number of cells for gene expression (assumed const acros perts)
#W a matrix with G rows, Nx cols of gene expression measurements
#columns of W= perturbations, rows of W= expression of genes
#y a vector of length Ns of average survival outcomes for each perturbation
#dir_wts a vector of length K of dirichlet prior weights for gene clusters

data {
  int<lower=1> Nx; //number of gene expression samples
  int<lower=1> Ns; //number of survival samples
  int<lower=1> G; //number of genes measured
  int<lower=1,upper=G> K; //number of gene clusters
  int<lower=1> L; //number of latent dimensions
  int<lower=1,upper=L> Ls; //number of latent dims relevant to survival
  int<lower=1,upper=G> xx[Nx]; //perturbation indicator for gene expression
  int<lower=1,upper=G> xs[Ns]; //pert indicator for survival
  int<lower=1> ncell_s[Ns];
  int<lower=1> ncell_x[Nx];
  matrix[G,Nx] W;
  vector[Ns] y;
  vector[K] dir_wts;
}
parameters {
  vector[L] phi[K];
  vector[L] psi[K];
  vector[L] a[G];
  vector[L] v[G];
  simplex[K] theta;
  vector[Ls] b;
  real mu_y;
  vector[G] mu_w;
  real<lower=0> sigma_y;
  real<lower=0> sigma_w;
  real<lower=0> sigma_b;
  real<lower=0> sigma_a;
  real<lower=0> sigma_v;
  real<lower=0> sigma_phi;
  real<lower=0> sigma_psi;
}
transformed parameters {
  vector[L] us[Ns];
  vector[L] ux[Nx];
  for(n in 1:Ns){
    us[n] = a[xs[n]];
  }
  for(n in 1:Nx){
    ux[n] = a[xx[n]];
  }
}
model {
  real log_theta[K] = log(theta); //cache
  // finite mixture prior
  theta~dirichlet(dir_wts);
  sigma_phi ~ gamma(2,1);
  sigma_psi ~ gamma(2,1);
  for(k in 1:K){
    phi[k] ~ normal(0,sigma_phi);
    psi[k] ~ normal(0,sigma_psi);
  }
  for (g in 1:G) {
    real lps_a[K] = log_theta;
    real lps_v[K] = log_theta;
    for (k in 1:K) {
      lps_a[k] += normal_lpdf(a[g] | phi[k], sigma_a);
      lps_v[k] += normal_lpdf(v[g] | psi[k], sigma_v);
    }
    target += log_sum_exp(lps_a)+log_sum_exp(lps_v);
  }
  // survival outcome
  mu_y~cauchy(0,5);
  b~cauchy(0,1);
  sigma_y~cauchy(0,1);
  for(n in 1:Ns){
    y[n] ~ normal(mu_y+dot_product(us[n],b),sigma_y);
  }
  //gene expression outcome
  mu_w~cauchy(0,5);
  sigma_w~cauchy(0,1);
  for(n in 1:Nx){
    for(g in 1:G){
      w[g,n]~ normal(mu_w[g] + dot_product(ux[n],v[g]), sigma_w);
    }
  }
}