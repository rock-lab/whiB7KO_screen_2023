functions {  
  /*
  * Alternative to neg_binomial_2_log_rng() that 
  * avoids potential numerical problems during warmup
  */
  int neg_binomial_2_log_safe_rng(real eta, real phi) {
    real gamma_rate = gamma_rng(phi, phi / exp(eta));
    if (gamma_rate >= exp(20.7)) gamma_rate = exp(20.7);
    return poisson_rng(gamma_rate);
  }
  
  // 
  // vector twoline_mean_pdf(int N, int[] C_A, real lambda, int[] guide_idx, vector generations, vector[] alpha_l, vector[] beta_l, vector[] gamma, vector[] beta_e) {
  //                       
  //   vector[N] eta;
  //   for(i in 1:N){
  //     if(generations[i] <= gamma[guide_idx[i]]){
  //       eta[i] = log(C_A[i] + lambda) + alpha_l[guide_idx[i]] + beta_l[guide_idx[i]] * generations[i];
  //     }
  //     else{
  //       eta[i] = log(C_A[i] + lambda) + ((alpha_l[guide_idx[i]] + beta_l[guide_idx[i]] * gamma[guide_idx[i]]) + 
  //                   ((beta_e[guide_idx[i]] *(generations[i] - gamma[guide_idx[i]]))));
  //     }
  //   }
  //   return eta;
  // }


}
data {
  int<lower=1> N;   // Number of data-points
  int<lower=1> J;   // Number of guides
  
  int<lower=1, upper=J> guide_idx[N];   // Guide index
  
  int<lower=0> C_B[N];
  vector<lower=0>[N] C_A;
  vector<lower=0>[N] generations;
  vector[J] s;      // vector with predicted strengths
  

  //
  real<lower=0> inv_phi_mean;
  real<lower=0> inv_phi_sd;
  
  
  // PRIORS
  
  real alpha_l_mean; real alpha_l_std;
  real beta_l_mean; real beta_l_std;
  real gamma_mean; real gamma_std;
  real beta_e_mean; real beta_e_std;
  
  real beta_max_mean; real beta_max_std;
  real beta_min_mean; real beta_min_std;
  real M_alpha; real M_beta;
  real H_mean; real H_std;
  
  // logistic
  
  
}
parameters {
  
  vector[J] alpha_l;
  vector[J] beta_l;
  vector<lower=0.1, upper=25>[J] gamma;
  // vector[J] gamma;
  vector[J] beta_e;
  real<lower=0> inv_phi;
  
  vector<lower=0>[J] lambda;
  
  // real<lower=0> nu_beta_e;
  // vector<lower=0>[J] nu_beta_e;
  real<lower=0> beta_e_sigma;
  // vector<lower=0>[J]  beta_e_sigma;
  
  
  // real<lower=0, upper=1> theta;
  
  // LOGISTIC
  real<lower =-2, upper = 0> beta_max;
  real<lower=0, upper=1> M;
  // real beta_min;
  // real beta_max;
  // real M;
  real<lower=0> H;
  
}
transformed parameters {
  // real beta_min = (-(beta_max/(1+exp(H*M)))) / (1 - (1/((1+exp(H*M)))));
  real<lower=0> phi = 1 / inv_phi;

}
model {
  
  
  real eta;
  vector[J] gmu;
  
  
   for (i in 1:N) {
     
    if(generations[i] <= gamma[guide_idx[i]]){
      eta = log(C_A[i] + lambda[guide_idx[i]]) + ((alpha_l[guide_idx[i]] + beta_l[guide_idx[i]] * generations[i]) * log(2));
    }
    else{
      eta = log(C_A[i] + lambda[guide_idx[i]]) + (((alpha_l[guide_idx[i]] + beta_l[guide_idx[i]] * gamma[guide_idx[i]]) + 
                  ((beta_e[guide_idx[i]] *(generations[i] - gamma[guide_idx[i]])))) * log(2))  ;
    }
     
    // target += neg_binomial_2_log_lpmf(C_B[i] | eta, phi);
    
    C_B[i] ~ neg_binomial_2_log(eta, phi);
     
    // if (C_B[i] == 0)
    //   target += log_sum_exp(bernoulli_lpmf(1 | theta),
    //                         bernoulli_lpmf(0 | theta)
    //                           + neg_binomial_2_log_lpmf(C_B[i] | eta, phi));
    // else
    //   target += bernoulli_lpmf(0 | theta)
    //               + neg_binomial_2_log_lpmf(C_B[i] | eta, phi);
  }

  lambda ~ normal(1, 10);

  alpha_l ~ normal(alpha_l_mean, alpha_l_std);
  beta_l ~ normal(beta_l_mean, beta_l_std);
  gamma ~ normal(gamma_mean, gamma_std);
  // beta_e ~ normal(beta_e_mean, beta_e_std);
  
  // gmu = 0 + (beta_max - 0) ./ (1+exp(-H * (s-M)));
  gmu = 0 + (beta_max - 0) ./ (1+exp(-H * (s-M)));
  beta_e ~ normal(gmu, beta_e_sigma);
  
  beta_e_sigma ~ gamma(0.1, 1);
  // nu_beta_e ~ gamma(2, 0.1);
  // beta_e ~ student_t(nu_beta_e, gmu, beta_e_sigma);
  
  
  beta_max ~ normal(beta_max_mean, beta_max_std);
  // beta_min ~ normal(beta_min_mean, beta_min_std);
  // beta_min ~ normal(0.1, 0.5);
  M ~ beta(M_alpha, M_beta);
  H ~ normal(H_mean, H_std);  // slope
  
  

  
  // theta ~ beta(1,10);
  
  // lambda ~ normal(lambda_mu, lambda_std);

}
// generated quantities {
//   int y_rep[N];
//   for (n in 1:N) {
//     y_rep[n] = neg_binomial_2_log_safe_rng(eta[n], inv_phi);
//   }
// }


