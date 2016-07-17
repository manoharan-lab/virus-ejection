library(data.table)
library(rstan)

file.names = dir("~/Downloads/ejectiondata160513",pattern=".txt")
for (i in 1:length(file.names)) {
  dat = fread(file.names[i])
  dat[,series:=i]
  if (i == 1) {
    alldat = dat
  } else {
    alldat = rbind(alldat,dat)
  }
}




findstartcode = "
data {
int<lower=0> N; // total measurements
vector[N] t; // time
vector[N] mi; // measured intensity
real sigma_lownoise; //fudge factor to deal with the fact the exponential curve used is an approximation,
  //especially at lower values.
real transitionspeed; //this is a fudge to 'round off the corner' and make things differentiable.
  //The bigger it is, the less 'overshoot'; at 15, overshoot is around 1%.
}

parameters {
real<lower=0,upper=1> start; ejection time as portion of full time

real<lower=0.1> scale; //horizontal scale
real<lower=0> sigma_noise;
} 

transformed parameters {
vector[N] ti; // 'true' (modeled) intensity
for (i in 1:G) {
ti[i] = (Phi((start - i/N) * transitionspeed * scale) +
         Phi((i/N - start) * transitionspeed * scale) * exp((start - i/N) * scale));
}
}

model {
mi ~ normal(ti,sigma + (1-ti)^2 * sigma_lownoise);
} 
"




model1stan = stan_model(model_code=findstartcode)








hcode = "
data {
int<lower=0> N; // total measurements
vector[N] t; // time
vector[N] mi; // measured intensity
}

parameters {
vector[N] ti; // true intensity
real<lower=0,upper=1> start; ejection time portion of full time

real gamma_pow; //the power law for the friction
real k; //combined multiplicative constant on force
} 

transformed parameters {
  
  for (i in 1:N) {
    if (i/N > start)
  }
}

model {
int prev;
transitions_high ~ beta(.9,.1);
transitions_low ~ beta(.1,.9);


for (i in 1:G) {
  start[i] ~ neg_binomial(startalpha, startbeta);
  if (start[i] > last[i+1] - last[i]) {
      increment_log_prob(-mult * log(sum( groupMembers[i] .* lambdas) + 1/gsigma));
  }
}
for (i in 1:N) {
  if (i == start[s[i]]) {
    state[i] ~ bernoulli(0.001) + 1;
    ti[i] = 1;
  } else if (i > start[s[i]]) {
    if (ti > breakpoint) {
      state[i]-1 ~ bernoulli(transitions_high[state[i-1]]) + 1;
    }
  }
  state[]
lambdas[i] <- exp(mu + xw * x[i] + gmu[g[i]]);
y[i] ~ poisson(lambdas[i]);
}

} 
"


model1stan = stan_model(model_code=hcode)

