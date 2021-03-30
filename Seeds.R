Seeds = function(niter=10^4, init=c(0,0,0,0,10), prop.sd=c(1,1,1,1,1)){   #ajuster prop.sd selon acc.rates
  
  I = 21
  r = c(10, 23, 23, 26, 17, 5, 53, 55, 32, 46, 10, 8, 10, 8, 23, 0, 3, 22, 15, 32, 3)
  n = c(39, 62, 81, 51, 39, 6, 74, 72, 51, 79, 13, 16, 30, 28, 45, 4, 12, 41, 30, 51, 7)
  x1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  x2 = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
  sd = 1000
  shape = 10^-3
  rate = 10^-3
  b = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)  #à updater
  
  chain=matrix(NA, niter+1, 5)
  chain[1,] = init
  
  acc.rates = rep(0, 5)
  
  for (iter in 1:niter){
    current = chain[iter,]
    
    logit_bottom = current[1]+current[2]*x1+current[3]*x2+current[4]*x1*x2+b
    #if(is.na(logit_bottom)) print("NA")
    
    #Mise à jour de alpha_0
    prop = rnorm(1, current[1], prop.sd[1])
    
    logit_top = prop+current[2]*x1+current[3]*x2+current[4]*x1*x2+b
    top = sum(r*logit_top - n*log(1+exp(logit_top)) + log(current[5]^(shape-1))
              - (prop^2+current[2]^2+current[3]^2+current[4]^2+current[5]*rate))  
    bottom = sum(r*logit_bottom - n*log(1+exp(logit_bottom)) + log(current[5]^(shape-1))
                 - (current[1]^2+current[2]^2+current[3]^2+current[4]^2+current[5]*rate))
    acc.prob1 = exp(top - bottom)
    #print(acc.prob1)
    #print(iter)
    print(b)
    
    if (runif(1) < acc.prob1){
      current[1] = prop
      acc.rates[1] = acc.rates[1] + 1
    }
    
    #Mise à jour de alpha_1
    prop = rnorm(1, current[2], prop.sd[2])
    
    logit_top = current[1]+prop*x1+current[3]*x2+current[4]*x1*x2+b
    top = sum(r*logit_top - n*log(1+exp(logit_top)) + log(current[5]^(shape-1))
              - (current[1]^2+prop^2+current[3]^2+current[4]^2+current[5]*rate))  
    bottom = sum(r*logit_bottom - n*log(1+exp(logit_bottom)) + log(current[5]^(shape-1))
                 - (current[1]^2+current[2]^2+current[3]^2+current[4]^2+current[5]*rate))
    acc.prob2 = exp(top - bottom)
    
    if (runif(1) < acc.prob2){
      current[2] = prop
      acc.rates[2] = acc.rates[2] + 1
    }
    
    #Mise à jour de alpha_2
    prop = rnorm(1, current[3], prop.sd[3])
    
    logit_top = current[1]+current[2]*x1+prop*x2+current[4]*x1*x2+b
    top = sum(r*logit_top - n*log(1+exp(logit_top)) + log(current[5]^(shape-1))
              - (current[1]^2+current[2]^2+prop^2+current[4]^2+current[5]*rate))  
    bottom = sum(r*logit_bottom - n*log(1+exp(logit_bottom)) + log(current[5]^(shape-1))
                 - (current[1]^2+current[2]^2+current[3]^2+current[4]^2+current[5]*rate))
    acc.prob3 = exp(top - bottom)
    
    if (runif(1) < acc.prob3){
      current[3] = prop
      acc.rates[3] = acc.rates[3] + 1
    }
    
    #Mise à jour de alpha_12
    prop = rnorm(1, current[4], prop.sd[4])
    
    logit_top = current[1]+current[2]*x1+current[3]*x2+prop*x1*x2+b
    top = sum(r*logit_top - n*log(1+exp(logit_top)) + log(current[5]^(shape-1))
              - (current[1]^2+current[2]^2+current[3]^2+prop^2+current[5]*rate))  
    bottom = sum(r*logit_bottom - n*log(1+exp(logit_bottom)) + log(current[5]^(shape-1))
                 - (current[1]^2+current[2]^2+current[3]^2+current[4]^2+current[5]*rate))
    acc.prob4 = exp(top - bottom)
    
    if (runif(1) < acc.prob4){
      current[4] = prop
      acc.rates[4] = acc.rates[4] + 1
    }
    
    #Mise à jour de tau
    prop = rlnorm(1, log(current[5]), prop.sd[5])
    
    logit_top = current[1]+current[2]*x1+current[3]*x2+current[4]*x1*x2+b
    top = sum(r*logit_top - n*log(1+exp(logit_top)) + log(current[5]^(shape-1))
              - (current[1]^2+current[2]^2+current[3]^2+current[4]^2+prop*rate))  
    bottom = sum(r*logit_bottom - n*log(1+exp(logit_bottom)) + log(current[5]^(shape-1))
                 - (current[1]^2+current[2]^2+current[3]^2+current[4]^2+current[5]*rate))
    acc.prob5 = exp(top - bottom)*(prop/current[5])
    
    if (runif(1) < acc.prob5){
      current[5] = prop
      acc.rates[5] = acc.rates[5] + 1
    }
    
    #Mise à jour de b
    b = rnorm(21,0,1/sqrt(current[5]))
    
    chain[iter+1,] = current
  }
  
  chain[,5]=1/sqrt(chain[,5])
  
  return(list(chain=chain, acc.rates=acc.rates/niter))
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Etude
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++

acc.rates = Seeds()$acc.rates
chain = Seeds()$chain

par(mfrow=c(1,5))
para = c("alpha[0]","alpha[1]","alpha[2]","alpha[12]","tau")
for (i in 1:5){
  plot(chain[,i], type="l", xlab="Iterations", ylab="", main=para[i])
}

colMeans(chain)

#Init
alpha0 <- 0
alpha1 <- 0
alpha2 <- 0
alpha12 <- 0
tau <- 10
b <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0)



data {
  int<lower=0> I;
  int<lower=0> n[I];
  int<lower=0> N[I];
  vector[I] x1; 
  vector[I] x2; 
}

transformed data {
  vector[I] x1x2;
  x1x2 <- x1 .* x2;
} 
parameters {
  real alpha0;
  real alpha1;
  real alpha12;
  real alpha2;
  real<lower=0> tau;
  vector[I] b;
}
transformed parameters {
  real<lower=0> sigma;
  sigma  <- 1.0 / sqrt(tau);
}
model {
  alpha0 ~ normal(0.0,1.0E3);
  alpha1 ~ normal(0.0,1.0E3);
  alpha2 ~ normal(0.0,1.0E3);
  alpha12 ~ normal(0.0,1.0E3);
  tau ~ gamma(1.0E-3,1.0E-3);
  
  b ~ normal(0.0, sigma);
  n ~ binomial_logit(N, alpha0 + alpha1 * x1 + alpha2 * x2 + alpha12 * x1x2 + b);
}
