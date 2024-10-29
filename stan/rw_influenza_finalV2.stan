// Copyright 2024 Oliver Eales
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


data {
  int num_data;                       // number of data points
  int num_path;                       // number of pathogens
  int Y[num_data];                    // daily number of 'cases'
  int P1[num_path-1, num_data];       // daily number of lab tests positive for influenza A (1st entry) and all other pathogens
  int P2[2, num_data];                // daily number of influenza A H3N2, and influenza A H1N1
  int week_effect;          // Number of days in day of week effect? 1=none, 2=weekends?, 7=all days
  int DOW[num_data];        // integer of day of the week
  int<lower = 0, upper = 2> cov_structure; //0 is tau[1], 1 is tau[num_path], 2 is Sigma[num_path, num_path]
  int<lower = 0, upper = 1> noise_structure; //0 is only includes observation noise (same between pathogens), 1 includes noise in individual pathogens as well
}

transformed data {
  int cols;
  int rows;
  
  if(noise_structure==1){
    cols = num_path;
    rows = num_data;
    
  }
  else{
    cols = 0;
    rows = 0;
  }
}

parameters {
  matrix[num_path, num_data] a;           //First two rows/pathogens will be influenza A H3N2 and influenza A H1N1
  
  matrix<lower=0>[cols,rows] c;
  real<lower=0> eta[noise_structure==1 ? 1: 0 ];
  
  real<lower=0> phi;
  
  real<lower=0> tau[cov_structure==0 ? 1: cov_structure==1? num_path: 0 ];
  cov_matrix[cov_structure==2? num_path: 0] Sigma;
  
  simplex[week_effect] day_of_week_simplex;
}


model {
  //// Priors
  // Second-order random walk prior
  if(cov_structure==2){
    for(i in 3:num_data)
      a[,i] ~ multi_normal(2*a[,(i-1)] - a[,(i-2)], Sigma);
  } 
  
  else if(cov_structure==1){
    for(i in 1:num_path)
      a[i,3:num_data] ~ normal(2*a[i,2:(num_data-1)] - a[i,1:(num_data-2)], tau[i]);
  } 
  
  else{
    for(i in 1:num_path)
      a[i,3:num_data] ~ normal(2*a[i,2:(num_data-1)] - a[i,1:(num_data-2)], tau[1]);
  }
  
  // Uninformative priors on scale parameters?
  
  
  //// Likelihood
  real total_ILI[num_data];
  real total_A;
  vector[num_path-1] theta;
  
  if(noise_structure==1){
    for(i in 1:num_path){
      c[i,] ~ gamma(exp(a[i,])*eta[1], eta[1]);
    }
    
    for(i in 1:num_data){
      total_ILI[i] = sum(c[,i]);
      
      total_A   = sum(c[1:2,i]);
      
      theta[1]  = total_A;
      theta[2:(num_path-1)] = c[3:num_path,i];
      
      P1[,i] ~ multinomial(theta/total_ILI[i]);
      P2[,i] ~ multinomial(c[1:2,i]/total_A);

    }
  } 
  else{
    for(i in 1:num_data){
      total_ILI[i] = sum(exp(a[,i]));
      
      total_A   = sum(exp(a[1:2,i]));
      
      theta[1]  = total_A;
      theta[2:(num_path-1)] = exp(a[3:num_path,i]);
      
      P1[,i] ~ multinomial(theta/total_ILI[i]);
      P2[,i] ~ multinomial(exp(a[1:2,i])/total_A);
    }
  }
  
  // This assumes there is some noise in the number of symptomatic cases for each pathogen individually (with shared parameter eta)
  
  // Total number of cases (Y[i]) is negative-binomially distributed
  // Proportion of each pathogen (influenza A, and others) is multinomially distributed
  // Proportion of influenza A H3N2 and influenza A H1N1 is multinomially distributed

  
  if(week_effect==1){
    for(i in 1:num_data){
      
      Y[i] ~ neg_binomial(total_ILI[i]*phi, phi);
      
    }
  }
  else{
    for(i in 1:num_data){
      
      Y[i] ~ neg_binomial(total_ILI[i]*phi*week_effect*day_of_week_simplex[DOW[i]], phi);
      
    }
  }
}
