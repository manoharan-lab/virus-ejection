library(data.table)
library(rstan)
library(ggplot2)
library(reshape2)
library(loo)

onedat = fread("../tests/test_data/result_151119chamber3_5_X7_294,347.txt")

names(onedat) = c("time","intensity")  #,"fitted")
standat = list(N=dim(onedat)[1],
               mi = onedat[,intensity],
               sigma_lownoise = 0.2,
               transitionspeed = 15.0)

findstartcode = "
data {
int<lower=0> N; // total measurements
//vector[N] t; // time - ignored
vector[N] mi; // measured intensity
real sigma_lownoise; //fudge factor to deal with the fact the exponential curve used is an approximation,
  //especially at lower values.
real transitionspeed; //this is a fudge to 'round off the corner' and make things differentiable.
  //The bigger it is, the less 'overshoot'; at 15, overshoot is around 1%.
}

parameters {
real<lower=0,upper=1> start; //ejection time as portion of full time

real<lower=0.1> scale; //horizontal scale
real<lower=0> sigma_noise;
} 

transformed parameters {
  vector[N] ti; // 'true' (modeled) intensity
  for (i in 1:N) {
    
    
    ti[i] <- (Phi((start - i/(N*1.0)) * transitionspeed * scale) + 
              Phi((i/(N*1.0) - start) * transitionspeed * scale) * exp((start - i/(N*1.0)) * scale));
  }
}

model {
mi ~ normal(ti,sigma_noise + (1-ti) * sigma_lownoise);
} 
"
model1stan = stan_model(model_code=findstartcode)

b = optimizing(model1stan, data = standat)

onedat[,fitted:=b$par[4:3753]]

b$par[1:3]

mdat = melt(onedat,"time")
ggplot(dat=mdat,aes(x=time,y=value,color=variable)) + geom_point() 


onedat[,resid:=intensity-fitted]
scatter.smooth(onedat[,resid],lpars=list(col="red"))
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
base = "~/Downloads/ejectiondata160513"
file.names = dir(base,pattern=".txt")
for (i in 1:length(file.names)) {
  name=file.names[i]
  fullname=paste(base,name,sep="/")
  print(fullname)
  dat = fread(fullname)
  dat[,series:=i]
  if (i == 1) {
    alldat = dat
  } else {
    alldat = rbind(alldat,dat)
  }
}
names(alldat) = c("time","intensity","series")
ggplot(dat=alldat,aes(x=time,y=intensity,color=as.factor(series)))+geom_point()
for (i in 1:18) {
  plot(alldat[series==i,intensity])
}
N=dim(alldat)[1]
allstandat = list(N = N,
                  low_pause = 5,
                  hi_pause = 150,
                  G = max(alldat[,series]),
                  mi = alldat[,intensity],
                  gid = alldat[,series],
                  d_min = 2.75, # from Purohit, Kondev, and Phillips (2002)
                  d_max = 4.85, # ditto
                  log_nonpauselambda = 0,
                  transitionspeed_low = .05,
                  scale_mean_hyperprior = log(standat$N / b$par[2]), #log(308)
                  scale_sigma_hyperprior = 10, #this is pretty large/"weakly informative"
                  
                  sameTransitions = 0,
                  noPause = 0,
                  sameVar = 0
                  )
hcode = "
data {
int<lower=0> N; // total measurements
int<lower=0> G; // total groups
vector[N] mi; // measured intensity
int gid[N]; // group id

// fudge factors and constants
real low_pause;
real hi_pause;
real transitionspeed_low;
real scale_mean_hyperprior;
real scale_sigma_hyperprior;
real d_min;
real d_max;
//real log_nonpauselambda; //big, so 'non' pauses are very short
//this is a fudge to 'round off the corner' and make things differentiable.
//The bigger it is, the less 'overshoot'; at 15, overshoot is around 1%.


// tuneable fudge factors, for suppressing certain parameters
int sameTransitions; //boolean; if true, put strong prior that sd_log_transitionspeed -> 0.
int noPause; //boolean; if true, pause_length -> 0
int sameVar; //boolean; if true, sd_log_sigma_noise -> 0
}

transformed data {
  int starts[G];
  int lengths[G];
  
  starts[1] <- 1;
  for (i in 2:N){
    if (gid[i] != gid[i-1]) {
      starts[gid[i]] <- i;
      lengths[gid[i-1]] <- i - starts[gid[i-1]];
    }
  }
  lengths[G] <- N + 1 - starts[G];
}


parameters {


real<lower=0> mean_log_sigma_noise; //should vary by group; do that later
real<lower=0> sd_log_sigma_noise; //should vary by group; do that later

real<lower=0> mean_log_transitionspeed; //should vary by group; do that later
real<lower=0> sd_log_transitionspeed; //should vary by group; do that later

real<lower=2> scale; //log of slope of linear multiplier
//real<lower=0,upper=1> pauseprob;
real<lower=0> log_pauselambda;
//real prepause_alpha;
//real prepause_beta;
//real<lower=0> gamma_pow; //the power law for the friction, when not pausing
real<lower=transitionspeed_low> transitionspeed[3,G];

real<lower=0> sigma_noise[G]; //ejection time portion of full time
real<lower=0,upper=1> start[G]; //ejection time portion of full time
vector<lower=low_pause,upper=hi_pause>[G] pause_start;
vector<lower=low_pause>[G] pause_length;

} 

transformed parameters {
  real elapsed;
  real prev;
  vector[N] sigma_noises; // 'true' (modeled) intensity
  vector[N] log_lik;
  vector[N] ti; // 'true' (modeled) intensity


  for (curg in 1:G) {
    int i;
    real force; 
    real basestart[2];
    //real friction; 
    real pause;
    ti[starts[curg]] <- 1;
    sigma_noises[starts[curg]] <- sigma_noise[curg];
    for (j in 1:(lengths[curg] - 1)) {
      i <- starts[curg] + j;
      sigma_noises[i] <- sigma_noise[curg];
      if (gid[i] == gid[i-1]) {
        basestart[1] <- lengths[curg] * start[curg];
        basestart[2] <- pause_length[curg];
        elapsed <- j - basestart[1] + min(basestart); //subtract pauselen: reparameterize for faster convergence
        prev <- ti[i-1];
        if (is_nan(prev)) {
          print(i);
          prev <- 1;
        }
        force <- Phi(elapsed * transitionspeed[1,curg]) / scale * prev;//exp(prev-1) * (1 - prev * multslope[curg]) ;
        //friction <- ((prev * (d_max - d_min) + d_min) ^ gamma_pow);
        pause <- (1 - Phi((elapsed - pause_start[curg]) * transitionspeed[2,curg])  *  
                      Phi((pause_start[curg] + pause_length[curg] - elapsed) * transitionspeed[3,curg]));
        ti[i] <- prev - force  * pause;// / friction;
        if (ti[i] < 0 && prev > 0) {
          print(prev,\"i=\",i,\" force=\",force);//,\" frictione=\",friction,\" pause=\",pause);
        }
        //if (is_nan(force) || is_nan(pause)) {
          //print(\"i=\",i,\" force=\",force[i],\" frictione=\",friction[i],\" pause=\",pause[i]);
          //print(\"prev=\",prev,\" elapsed=\",elapsed,\" pause_start[gid[i]]=\",pause_start[gid[i]],\" multslope[gid[i]]=\",multslope[gid[i]]);
        //}
      } else {
        ti[i] <- 1;
      }
    }
    //print(\"curg=\",curg,\" force=\",force,\" pause=\",pause, \" pausest=\",pause_start[curg],\" pausel=\",pause_length[curg]);
  }

  for (n in 1:N) log_lik[n] <- normal_log(mi[n], ti[n], sigma_noises[n]);
}

model {
  //hyperprior on scale; default flat is good for restart
  scale ~ lognormal(scale_mean_hyperprior,scale_sigma_hyperprior);

  //not overprecise hyperparameters
  increment_log_prob(log(sd_log_transitionspeed * sd_log_sigma_noise));

  //group level variables from hyperparameters
  sigma_noise ~ lognormal(mean_log_sigma_noise,sd_log_sigma_noise);

  for (curg in 1:G) {
    transitionspeed[1,curg] ~ lognormal(mean_log_transitionspeed,sd_log_transitionspeed);
    transitionspeed[2,curg] ~ lognormal(mean_log_transitionspeed,sd_log_transitionspeed);
    transitionspeed[3,curg] ~ lognormal(mean_log_transitionspeed,sd_log_transitionspeed);
  }

  //scale ~ lognormal(mean_logscale,sigma_logscale);
  //multslope ~ lognormal(mean_logmultslope,sigma_logmultslope);
  pause_start ~ uniform(low_pause,hi_pause);//gamma(exp(prepause_alpha),exp(prepause_beta)*.01);
  pause_length ~ exponential(exp(-log_pauselambda));
  //for (i in 1:G) {
  //  increment_log_prob(log(
  //            pauseprob * exp(log_pauselambda + pause_length[i] * exp(log_pauselambda)) +
  //            (1-pauseprob) * exp(log_nonpauselambda + pause_length[i] * exp(log_nonpauselambda))));
  //}
  //print(\"scale=\",scale,\" pauseprob=\",pauseprob,\" pause_length[1]=\",pause_length[1],\" pause_start[=\",pause_start[1]);


  //print(sigma_noise);
  //print(sigma_noise,\" \",ti);

  mi ~ normal(ti,sigma_noises);

}

"


model2stan = stan_model(model_code=hcode)

fullfit = optimizing(model2stan, data = allstandat)
#fullfit = stan(model_code=hcode, data = allstandat)

l = length(fullfit$par)
alldat[,fitted:=fullfit$par[(l-N+1):l]]
alldat[,index:=1:.N]
malldat = melt(alldat,c("time","series","index"))
ggplot(dat=malldat,aes(x=index,y=value,color=variable)) + geom_point() 


fullfit$par[1:80]

loglik = fullfit$par[gsub("\\[.*","",names(fullfit$par)) == "log_lik"]

nparam = 6 + 7 * allstandat$G
bicfull = - 2 * sum(loglik) + nparam * (log(allstandat$N/2/pi))
bicfull

allstandat$sameTransitions = 1
allstandat$noPause = 0
allstandat$sameVar = 0
sameTransFit = optimizing(model2stan, data = allstandat)
stloglik = sameTransFit$par[gsub("\\[.*","",names(fullfit$par)) == "log_lik"]
stnparam = 5 + 4 * allstandat$G
stbic = - 2 * sum(stloglik) + stnparam * (log(allstandat$N/2/pi))
stbic


2










































#backup version below; ignore.



N=dim(alldat)[1]
allstandat = list(N = N,
                  low_pause = 5,
                  hi_pause = 150,
                  G = max(alldat[,series]),
                  mi = alldat[,intensity],
                  gid = alldat[,series],
                  d_min = 2.75, # from Purohit, Kondev, and Phillips (2002)
                  d_max = 4.85, # ditto
                  log_nonpauselambda = 0,
                  transitionspeed = 1,
                  scale_mean_hyperprior = log(standat$N / b$par[2]), #log(308)
                  scale_sigma_hyperprior = 10 #this is pretty large/"weakly informative"
)


hcode = "
data {
int<lower=0> N; // total measurements
int<lower=0> G; // total groups
vector[N] mi; // measured intensity
int gid[N]; // group id

// fudge factors and constants
real low_pause;
real hi_pause;
real scale_mean_hyperprior;
real scale_sigma_hyperprior;
real d_min; //Has to do with friction; currently ignored
real d_max; //Has to do with friction; currently ignored
//real log_nonpauselambda; //big, so 'non' pauses are very short
real transitionspeed; //this is a fudge to 'round off the corner' and make things differentiable.
//The bigger it is, the less 'overshoot'; at 15, overshoot is around 1%.
}

transformed data {
int starts[G];
int lengths[G];

starts[1] <- 1;
for (i in 2:N){
if (gid[i] != gid[i-1]) {
starts[gid[i]] <- i;
lengths[gid[i-1]] <- i - starts[gid[i-1]];
}
}
lengths[G] <- N + 1 - starts[G];
}


parameters {


real<lower=0> sigma_noise; //should vary by group; do that later
real<lower=2> scale; //log of slope of linear multiplier
//real<lower=0,upper=1> pauseprob;
real<lower=0> log_pauselambda;
//real prepause_alpha;
//real prepause_beta;
//real<lower=0> gamma_pow; //the power law for the friction, when not pausing

real<lower=0,upper=1> start[G]; //ejection time portion of full time
vector<lower=low_pause,upper=hi_pause>[G] pause_start;
vector<lower=low_pause>[G] pause_length;

} 

transformed parameters {
real elapsed;
real prev;
vector[N] ti; // 'true' (modeled) intensity

for (curg in 1:G) {
int i;
real force; 
real basestart[2];
//real friction; 
real pause;
ti[starts[curg]] <- 1;
for (j in 1:(lengths[curg] - 1)) {
i <- starts[curg] + j;
if (gid[i] == gid[i-1]) {
basestart[1] <- lengths[curg] * start[curg];
basestart[2] <- pause_length[curg];
elapsed <- j - basestart[1] + min(basestart); //subtract pauselen: reparameterize for faster convergence
prev <- ti[i-1];
if (is_nan(prev)) {
print(i);
prev <- 1;
}
force <- Phi(elapsed * transitionspeed) / scale * prev;//exp(prev-1) * (1 - prev * multslope[curg]) ;
//friction <- ((prev * (d_max - d_min) + d_min) ^ gamma_pow);
pause <- (1 - Phi((elapsed - pause_start[curg]) * transitionspeed)  *  
Phi((pause_start[curg] + pause_length[curg] - elapsed) * transitionspeed));
ti[i] <- prev - force  * pause;// / friction;
if (ti[i] < 0 && prev > 0) {
print(prev,\"i=\",i,\" force=\",force);//,\" frictione=\",friction,\" pause=\",pause);
}
//if (is_nan(force) || is_nan(pause)) {
//print(\"i=\",i,\" force=\",force[i],\" frictione=\",friction[i],\" pause=\",pause[i]);
//print(\"prev=\",prev,\" elapsed=\",elapsed,\" pause_start[gid[i]]=\",pause_start[gid[i]],\" multslope[gid[i]]=\",multslope[gid[i]]);
//}
} else {
ti[i] <- 1;
}
}
//print(\"curg=\",curg,\" force=\",force,\" pause=\",pause, \" pausest=\",pause_start[curg],\" pausel=\",pause_length[curg]);
}
}

model {
//hyperprior on scale; default flat is good for restart
scale ~ lognormal(scale_mean_hyperprior,scale_sigma_hyperprior);

//not overprecise hyperparameters
//increment_log_prob(log(sigma_logscale * sigma_logmultslope));

//group level variables from hyperparameters
//scale ~ lognormal(mean_logscale,sigma_logscale);
//multslope ~ lognormal(mean_logmultslope,sigma_logmultslope);
pause_start ~ uniform(low_pause,hi_pause);//gamma(exp(prepause_alpha),exp(prepause_beta)*.01);
pause_length ~ exponential(exp(-log_pauselambda));
//for (i in 1:G) {
//  increment_log_prob(log(
//            pauseprob * exp(log_pauselambda + pause_length[i] * exp(log_pauselambda)) +
//            (1-pauseprob) * exp(log_nonpauselambda + pause_length[i] * exp(log_nonpauselambda))));
//}
//print(\"scale=\",scale,\" pauseprob=\",pauseprob,\" pause_length[1]=\",pause_length[1],\" pause_start[=\",pause_start[1]);


//print(sigma_noise);
//print(sigma_noise,\" \",ti);

mi ~ normal(ti,sigma_noise);

if (sameTransitions) sd_log_transitionspeed ~ normal(0,0.01);
if (noPause) pause_length  ~ normal(0,0.01);
if (sameVar) sd_log_sigma_noise ~ normal(0,0.01);
}
"


model2stan = stan_model(model_code=hcode)

fullfit = optimizing(model2stan, data = allstandat)
#fullfit = stan(model_code=hcode, data = allstandat)

l = length(fullfit$par)
alldat[,fitted:=fullfit$par[(l-N+1):l]]
alldat[,index:=1:.N]
malldat = melt(alldat,c("time","series","index"))
ggplot(dat=malldat,aes(x=index,y=value,color=variable)) + geom_point() 
















#Here's a model with the kitchen sink (friction term, non-working force function). It fails with 
#"Error evaluating model log probability: Non-finite gradient."


hcode = "
data {
int<lower=0> N; // total measurements
int<lower=0> G; // total groups
vector[N] mi; // measured intensity
int gid[N]; // group id

// fudge factors and constants
real scale_mean_hyperprior;
real scale_sigma_hyperprior;
real d_min;
real d_max;
real log_nonpauselambda; //big, so 'non' pauses are very short
real transitionspeed; //this is a fudge to 'round off the corner' and make things differentiable.
//The bigger it is, the less 'overshoot'; at 15, overshoot is around 1%.
}

transformed data {
vector[G] starts;
vector[G] lengths;

starts[1] <- 1;
for (i in 2:N){
if (gid[i] != gid[i-1]) {
starts[gid[i]] <- i;
lengths[gid[i-1]] <- i - starts[gid[i-1]];
}
}
lengths[G] <- N + 1 - starts[G];
}


parameters {


real<lower=0> sigma_noise;
real mean_logscale; //log of horizontal scale
real<lower=0> sigma_logscale; 
real mean_logmultslope; //log of slope of linear multiplier
real<lower=0> sigma_logmultslope; 
real<lower=0,upper=1> pauseprob;
real log_pauselambda;
real prepause_alpha;
real prepause_beta;
real<lower=0> gamma_pow; //the power law for the friction, when not pausing

real<lower=0,upper=1> start[G]; //ejection time portion of full time
vector<lower=0>[G] pause_start;
vector<lower=0>[G] pause_length;
vector<lower=2>[G] scale;
vector<lower=0,upper=.99>[G] multslope;

} 

transformed parameters {
real elapsed;
real prev;
vector[N] ti; // 'true' (modeled) intensity

ti[1] <- 1;
for (i in 2:N) {
real force; 
real friction; 
real pause;
int curg;
if (gid[i] == gid[i-1]) {
elapsed <- i - starts[gid[i]] - lengths[gid[i]] * start[gid[i]];
prev <- ti[i-1];
if (is_nan(prev)) {
print(i);
prev <- 1;
}
curg <- gid[i];
force <- Phi(elapsed * transitionspeed) / scale[curg] * prev;//exp(prev-1) * (1 - prev * multslope[curg]) ;
friction <- ((prev * (d_max - d_min) + d_min) ^ gamma_pow);
pause <- (1 - Phi((elapsed - pause_start[curg]) * transitionspeed)  *  
Phi((pause_start[curg] + pause_length[curg] - elapsed) * transitionspeed));
ti[i] <- prev - pause * force / friction;
if (ti[i] < 0 && prev > 0) {
print(prev,\"i=\",i,\" force=\",force,\" frictione=\",friction,\" pause=\",pause);
}
if (is_nan(force) || is_nan(friction) || is_nan(pause)) {
//print(\"i=\",i,\" force=\",force[i],\" frictione=\",friction[i],\" pause=\",pause[i]);
//print(\"prev=\",prev,\" elapsed=\",elapsed,\" pause_start[gid[i]]=\",pause_start[gid[i]],\" multslope[gid[i]]=\",multslope[gid[i]]);
}
} else {
ti[i] <- 1;
}
}
}

model {
//hyperprior on scale; default flat is good for restart
mean_logscale ~ normal(scale_mean_hyperprior,scale_sigma_hyperprior);

//not overprecise hyperparameters
increment_log_prob(log(sigma_logscale * sigma_logmultslope));

//group level variables from hyperparameters
scale ~ lognormal(mean_logscale,sigma_logscale);
multslope ~ lognormal(mean_logmultslope,sigma_logmultslope);
pause_start ~ gamma(exp(prepause_alpha),exp(prepause_beta));
for (i in 1:G) {
increment_log_prob(log(
pauseprob * exp(log_pauselambda + pause_length[i] * exp(log_pauselambda)) +
(1-pauseprob) * exp(log_nonpauselambda + pause_length[i] * exp(log_nonpauselambda))));
}


//print(sigma_noise);
//print(sigma_noise,\" \",ti);

mi ~ normal(ti,sigma_noise);
}
"
