functions {

// credit http://srmart.in/informative-priors-for-correlation-matrices-an-easy-approach/
vector lower_tri(matrix mat) {

int d = rows(mat);
int lower_tri_d = d * (d - 1) / 2;
vector[lower_tri_d] lower;
int count = 1;
for(r in 2:d) {
for(c in 1:(r - 1)) {
lower[count] = mat[r,c];
count += 1;
}
}
return(lower); 
}

// credit http://srmart.in/informative-priors-for-correlation-matrices-an-easy-approach/
real lkj_corr_point_lower_tri_lpdf(matrix rho, vector point_mu_lower, vector point_scale_lower) {

real lpdf = lkj_corr_lpdf(rho | 1) + normal_lpdf(lower_tri(rho) | point_mu_lower, point_scale_lower);
return(lpdf);
}


real lkj_corr_cholesky_point_lower_tri_two_lpdf(matrix cor_L, real point_mu_lower, real point_scale_lower) {
    real lpdf = lkj_corr_cholesky_lpdf(cor_L | 1);
    int d = rows(cor_L);
    matrix[d,d] cor = multiply_lower_tri_self_transpose(cor_L);
    lpdf += normal_lpdf(cor[2,1] | point_mu_lower, point_scale_lower);
    return(lpdf);
  }

real normal_lub_rng(real mu, real sigma, real lb, real ub) {
  real p_lb = normal_cdf(lb, mu, sigma);
  real p_ub = normal_cdf(ub, mu, sigma);
  real u = uniform_rng(p_lb, p_ub);
  real y = mu + sigma * Phi(u);
  return y;
}

real gamma_lb_rng(real a, real b, real lb) {

real count=1;
real epsilon_group = gamma_rng(a, b);

  while (epsilon_group < lb && count < 1000)
  epsilon_group = gamma_rng(a, b);
  count = count + 1;
  
return epsilon_group;

}

}

data{
int nAnalytes;	           // number of analytes
int nObs;		           // number of observations
int npH;                   // npH;
int analyte[nObs];	       // analyte indexes
int pHid[nObs];
int<lower=1> steps[nObs];   // steps for gradient retention time aproimation
vector[11] hplcparam[nObs]; // [tg, td, to, te, fio, fik, mod, pHo, alpha1, alpha2, (temp-25)/10]
int<lower=0> mod[nObs];     // MeOH==1, ACN==2 (repeats hplcparam(:,7))

vector[nAnalytes] logPobs; 

int<lower=0,upper=2> maxR;
int<lower=0,upper=2> R[nAnalytes];
ordered[maxR] pKaslit[nAnalytes];
vector[maxR] pKasliterror[nAnalytes];
vector[maxR] groupsA[nAnalytes];
vector[maxR] groupsB[nAnalytes];
vector[maxR+1] chargesA[nAnalytes];
vector[maxR+1] chargesB[nAnalytes];

vector[nObs] trobs; // observed retention factors 
}

transformed data {
int grainsize = 1;
int ind[nObs] = rep_array(1, nObs);
vector[3] point_mu_lower = [0.75,0.75,0.75]';       // mean priors for rho
vector[3] point_scale_lower = [0.125,0.125,0.125]';  // std priors for rho
}

parameters{
real foo;
corr_matrix[3] rho1;	 
cholesky_factor_corr[2] L2;	
}

model{
foo ~ normal(2.2, 1);
rho1 ~ lkj_corr_point_lower_tri(point_mu_lower, point_scale_lower);
L2 ~ lkj_corr_cholesky_point_lower_tri_two(0.75, 0.125);

}

generated quantities{
real logkwHat;	      
real S1mHat;	    
real S1aHat;         
real dlogkwHat[2];    
real dSmHat[2];      
real dSaHat[2];       
real<lower = 0> S2mHat;
real<lower = 0> S2aHat; 
vector[3] beta;         
real dlogkTHat;      
vector[2] alphaAHat; 
vector[2] alphaBHat;  
vector<lower = 0.>[3] omega;    	     
vector<lower = 0>[3] kappa;  
vector<lower = 0>[2] tau;    	
real<lower = 0> omegadlogkT;  
vector[2] apH;
real<lower = 0.> msigma; 
real<lower = 0> ssigma; 
corr_matrix[2] rho2;
 
rho2 = L2 * L2';

  logkwHat  = normal_rng(2.2, 2);
  S1mHat    = normal_rng(4, 1);
  S1aHat    = normal_rng(5, 1);
  dlogkwHat[1] = normal_rng(-1,0.125);
  dlogkwHat[2] = normal_rng(-1,0.125);
  dSmHat[1] = normal_rng(0,0.5);
  dSmHat[2] = normal_rng(0,0.5);
  dSaHat[1] = normal_rng(0,0.5);
  dSaHat[2] = normal_rng(0,0.5);
  S2mHat = lognormal_rng(-1.6,0.125);
  S2aHat = lognormal_rng(0.69,0.125);
  alphaAHat[1] = normal_rng(2,0.25);
  alphaBHat[1] = normal_rng(-1,0.25);
  alphaAHat[2] = normal_rng(2,0.25);
  alphaBHat[2] = normal_rng(-1,0.25);
  beta[1]  = normal_rng(1,0.125);
  beta[2]  = normal_rng(0.5,0.5);
  beta[3]  = normal_rng(0.5,0.5);
  omega[1] = abs(normal_rng(0,2));
  omega[2] = abs(normal_rng(0,2));
  omega[3] = abs(normal_rng(0,2));
  kappa[1] = abs(normal_rng(0,0.5));
  kappa[2] = abs(normal_rng(0,0.5));
  kappa[3] = abs(normal_rng(0,0.5));

  tau[1]     = abs(normal_rng(0,0.5));
  tau[2]     = abs(normal_rng(0,0.5));
  apH[1]  = normal_rng(0,0.1);
  apH[2]  = normal_rng(0,0.1);

  dlogkTHat   = normal_rng(-0.087,0.022);
  omegadlogkT = abs(normal_rng(0,0.022));
  msigma = abs(normal_rng(0,1));
  ssigma = abs(normal_rng(0,1));

}
