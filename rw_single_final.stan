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
  int num_data;              // number of data points
  int Y[num_data];           // daily number of 'cases'
  int week_effect;          // Number of days in day of week effect? 1=none, 2=weekends?, 7=all days
  int DOW[num_data];        // integer of day of the week
}

parameters {
  row_vector[num_data] a;
  
  real<lower=0> phi;
  
  real<lower=0> tau;
  
  simplex[week_effect] day_of_week_simplex;
}


model {
  //// Priors
  // Second-order random walk prior
  a[3:num_data] ~ normal(2*a[2:(num_data-1)] - a[1:(num_data-2)], tau);

  // Uninformative priors on scale parameters?
  
  //// Likelihood
  if(week_effect==1){
    Y ~ neg_binomial(exp(a)*phi, phi);
  } 
  else{
    for(i in 1:num_data)
      Y[i] ~ neg_binomial(exp(a)*phi*week_effect*day_of_week_simplex[DOW[i]], phi);
  }
  
}
