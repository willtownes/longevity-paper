//Full model for integrating gene expression with survival outcome

data {
  int<lower=1> Nx; //number of gene expression samples
  int<lower=1> Ns; //number of survival samples
  int<lower=1> G; //number of genes measured
  int<lower=1,upper=G> K; //number of gene clusters
  int<lower=1> L; //number of latent dimensions
  int<lower=1,upper=L-1> Ls; //number of latent dims relevant to survival
  int<lower=1,upper=G> xx[Nx]; //perturbation indicator for gene expression
  int<lower=1,upper=G> xs[Ns]; //pert indicator for survival
  //matrix[Ns,G] Xs;
  //matrix[G,Nx] Xxt;
  vector[Ns] ncell_s;
  real ncell_x;
  matrix[Nx,G] W; //transpose of usual genomics orientation
  vector[Ns] y;
  vector[K] dir_wts;
}
transformed data {
  int<lower=1,upper=L-1> Lx;
  Lx = L-Ls;
}
parameters {
  vector[Ls] phi_s[K];
  vector[Lx] phi_x[K];
  vector[Ls] psi_s[K];
  vector[Lx] psi_x[K];
  //vector[Ls] a_s[G];
  //vector[Lx] a_x[G];
  //vector[Ls] v_s[G];
  //vector[Lx] v_x[G];
  matrix[G,Ls] As;
  matrix[G,Lx] Ax;
  matrix[Ls,G] Vs;
  matrix[Lx,G] Vx;
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
//transformed parameters {
  //vector[Ls] us[Ns];
  //vector[L] ux[Nx];
  //for(n in 1:Ns){
  //  us[n] = a[xs[n],1:Ls];
  //}
  //for(n in 1:Nx){
  //  ux[n] = a[xx[n]];
  //}
//}
model {
  vector[K] log_theta = log(theta); //cache
  vector[Ns] yhat;
  // finite mixture prior
  theta~dirichlet(dir_wts);
  sigma_phi ~ gamma(2,1);
  sigma_psi ~ gamma(2,1);
  for(k in 1:K){
    phi_s[k] ~ normal(0,sigma_phi);
    phi_x[k] ~ normal(0,sigma_phi);
    psi_s[k] ~ normal(0,sigma_psi);
    psi_x[k] ~ normal(0,sigma_psi);
  }
  for (g in 1:G) {
    //real lps_a[K] = log_theta;
    //real lps_v[K] = log_theta;
    vector[K] lps_a_s = log_theta;
    vector[K] lps_v_s = log_theta;
    vector[K] lps_a_x = log_theta;
    vector[K] lps_v_x = log_theta;
    for (k in 1:K) {
      lps_a_s[k] += normal_lpdf(a_s[g] | phi_s[k], sigma_a);
      lps_v_s[k] += normal_lpdf(v_s[g] | psi_s[k], sigma_v);
      lps_a_x[k] += normal_lpdf(a_x[g] | phi_x[k], sigma_a);
      lps_v_x[k] += normal_lpdf(v_x[g] | psi_x[k], sigma_v);
    }
    target += log_sum_exp(lps_a_s)+log_sum_exp(lps_v_s);
    target += log_sum_exp(lps_a_x)+log_sum_exp(lps_v_x);
  }
  // survival outcome
  mu_y~cauchy(0,5);
  b~cauchy(0,1);
  sigma_y~cauchy(0,1);
  for(n in 1:Ns){
    yhat[n] = mu_y + dot_product(a_s[xs[n]],b);
  }
  y ~ normal(yhat, sigma_y ./ncell_s);
  //gene expression outcome
  mu_w~cauchy(0,5);
  sigma_w~cauchy(0,1);
  for(n in 1:Nx){
    vector[G] w_hat;
    for(g in 1:G){
      w_hat[g] = dot_product(a_s[xx[n]],v_s[g])+dot_product(a_x[xx[n]],v_x[g]);
    }
    W[:,n]~ normal(mu_w + w_hat, sigma_w/ncell_x);
  }
}
