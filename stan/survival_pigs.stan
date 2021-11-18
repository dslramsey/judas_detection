functions {
 
 real average_p(real x, real xc, real[] theta, real[] x_r, int[] x_i) {     
  real b0 = theta[1];
  real b1 = theta[2];
  real b2 = theta[3];
  return (x * (1-exp(-exp(b0 + b1*x + b2))));
  } 
}

data {
  int<lower=0> N;
  int<lower=0> M;
  int<lower=0> K;
  vector<lower=0>[N] np;
  int id1[N];
  int id2[N];
  int<lower=0> time[N];
  matrix[N, K] X;
  int<lower=0> ev[N];
  real maxD[M];
  int<lower=0,upper=1> sex[M];
}
transformed data {
  real x_r[0];
  int x_i[0];
}
parameters {
  real b1;
  vector[K] beta;
  real<lower=0> sigma_id1;
  real<lower=0> sigma_id2;
  vector[M] id1_raw;
  vector[M] id2_raw;
}

transformed parameters {
  real q[N];
  real<lower=0, upper=1> p[N];
  vector[M] eps_id1;
  vector[M] eps_id2;
  
  eps_id1 = sigma_id1 * id1_raw;
  eps_id2 = sigma_id2 * id2_raw;
  
   for(i in 1:N){
      q[i] = b1 + X[i] * beta + eps_id1[id1[i]] + eps_id2[id2[i]];
      p[i] = inv_cloglog(q[i]);
    }
  
}

model {
  
  b1 ~ normal(0, 5);
  id1_raw ~ std_normal();
  id2_raw ~ std_normal(); 
  beta ~ normal(0, 5);  
  sigma_id1 ~ student_t(4, 0, 2); 
  sigma_id2 ~ student_t(4, 0, 2);
  
  for(i in 1:N) 
    target += binomial_lpmf(ev[i] | time[i], p[i]);
   
}

generated quantities {
  vector[M] beta0;
  int yrep[N];
  real pave[M];
  real intc[M];
  
  beta0 = b1 + eps_id1;
  for(n in 1:N) yrep[n] = binomial_rng(time[n], p[n]);
  for(i in 1:M) {
    intc[i] = integrate_1d(average_p, 0, maxD[i], {beta0[i], beta[1], beta[2]*sex[i]}, x_r, x_i, 1e-8);
    pave[i] = intc[i] * 2/maxD[i]^2;
  }
}
