Seeds = function(niter=10^4, init=c(-0.55,0.1,1.3,-0.8,10000), prop.sd=c(0.4,0.5,0.5,0.7), prop.sd_b=rep(0.05,21)){   #ajuster prop.sd pour approcher 25% de acc.rates
  
  I = 21
  r = c(10, 23, 23, 26, 17, 5, 53, 55, 32, 46, 10, 8, 10, 8, 23, 0, 3, 22, 15, 32, 3)
  n = c(39, 62, 81, 51, 39, 6, 74, 72, 51, 79, 13, 16, 30, 28, 45, 4, 12, 41, 30, 51, 7)
  x1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  x2 = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
  sd = 1000
  shape = 10^-3
  rate = 10^-3
  b = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)  #à updater
  
  chain=matrix(NA, niter+1, 26)
  chain[1,1:5] = init
  chain[1,6:26] = b
  
  acc.rates = rep(0, 25)
  
  for (iter in 1:niter){
    current = chain[iter,]
    
    logit_bottom = current[1]+current[2]*x1+current[3]*x2+current[4]*x1*x2+current[6:26]
    
    #Mise à jour de alpha_0
    prop = rnorm(1, current[1], prop.sd[1])
    
    logit_top = prop+current[2]*x1+current[3]*x2+current[4]*x1*x2+current[6:26]
    top = sum(r*logit_top - n*log(1+exp(logit_top)) 
              - (prop^2+current[2]^2+current[3]^2+current[4]^2)/(2*sd^2))  
    bottom = sum(r*logit_bottom - n*log(1+exp(logit_bottom))
                 - (current[1]^2+current[2]^2+current[3]^2+current[4]^2)/(2*sd^2))
    acc.prob1 = exp(top - bottom)
    
    if (runif(1) < acc.prob1){
      current[1] = prop
      acc.rates[1] = acc.rates[1] + 1
    }
    
    #Mise à jour de alpha_1
    prop = rnorm(1, current[2], prop.sd[2])
    
    logit_top = current[1]+prop*x1+current[3]*x2+current[4]*x1*x2+current[6:26]
    top = sum(r*logit_top - n*log(1+exp(logit_top))
              - (current[1]^2+prop^2+current[3]^2+current[4]^2)/(2*sd^2))  
    bottom = sum(r*logit_bottom - n*log(1+exp(logit_bottom))
                 - (current[1]^2+current[2]^2+current[3]^2+current[4]^2)/(2*sd^2))
    acc.prob2 = exp(top - bottom)
    
    if (runif(1) < acc.prob2){
      current[2] = prop
      acc.rates[2] = acc.rates[2] + 1
    }
    
    #Mise à jour de alpha_2
    prop = rnorm(1, current[3], prop.sd[3])
    
    logit_top = current[1]+current[2]*x1+prop*x2+current[4]*x1*x2+current[6:26]
    top = sum(r*logit_top - n*log(1+exp(logit_top))
              - (current[1]^2+current[2]^2+prop^2+current[4]^2)/(2*sd^2))  
    bottom = sum(r*logit_bottom - n*log(1+exp(logit_bottom))
                 - (current[1]^2+current[2]^2+current[3]^2+current[4]^2)/(2*sd^2))
    acc.prob3 = exp(top - bottom)
    
    if (runif(1) < acc.prob3){
      current[3] = prop
      acc.rates[3] = acc.rates[3] + 1
    }
    
    #Mise à jour de alpha_12
    prop = rnorm(1, current[4], prop.sd[4])
    
    logit_top = current[1]+current[2]*x1+current[3]*x2+prop*x1*x2+current[6:26]
    top = sum(r*logit_top - n*log(1+exp(logit_top))
              - (current[1]^2+current[2]^2+current[3]^2+prop^2)/(2*sd^2))  
    bottom = sum(r*logit_bottom - n*log(1+exp(logit_bottom))
                 - (current[1]^2+current[2]^2+current[3]^2+current[4]^2)/(2*sd^2))
    acc.prob4 = exp(top - bottom)
    
    if (runif(1) < acc.prob4){
      current[4] = prop
      acc.rates[4] = acc.rates[4] + 1
    }
    
    #Mise à jour de tau
    current[5] = rgamma(1, shape=shape+I/2, rate=rate+sum(b^2)/2)
    
    #Mise à jour des b_i
    prop_b = NULL
    vec = current[6:26]
    
    for (i in 1:I){
      prop_b[i] = rnorm(1, current[i+5], prop.sd_b[i])
      vec[i] = prop_b[i]
      
      logit_top = current[1]+current[2]*x1+current[3]*x2+current[4]*x1*x2+vec
      top = -(prop_b[i]^2)/(2*(1/current[5])) + sum(r*logit_top - n*log(1+exp(logit_top)))
      bottom = -(current[i+5]^2)/(2*(1/current[5])) + sum(r*logit_bottom - n*log(1+exp(logit_bottom)))
      acc.prob = exp(top-bottom)
      
      if (runif(1) < acc.prob){
        current[i+5] = prop_b[i]
        acc.rates[i+4] = acc.rates[i+4] + 1
      }
    }
    
    #Mise à jour globale des chaînes
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
para = c("alpha[0]","alpha[1]","alpha[2]","alpha[12]","sigma")
for (i in 1:5){
  plot(chain[,i], type="l", xlab="Iterations", ylab="", main=para[i])
}

colMeans(chain)

acc.rates

par(mfrow=c(3,7))
for (i in 6:26){
  plot(chain[,i], type="l", xlab="Iterations", ylab="")
}

#Burning

chain_burn = chain[-(1:1000),]

par(mfrow=c(1,5))
para_b = c("b[1]","b[2]","b[3]","b[4]","b[5]")
for (i in 6:10){
  plot(chain_burn[,i], type="l", xlab="Iterations", ylab="", main=para_b[i-5])
}

#ACF & Elagage

par(mfrow=c(1,5))
for (i in 1:5){
  acf(chain[,i], lag.max=100, main=para[i])
}

chain_elag = chain[seq(1, nrow(chain), by = 20), 1:5]
par(mfrow=c(2,5))
for (i in 1:5){
  plot(chain_elag[,i], type="l", xlab="Iterations", ylab="", main=para[i])
}
for (i in 1:5){
  acf(chain_elag[,i], lag.max=100, main="")
}

colMeans(chain_elag)

par(mfrow=c(1,5))
for (i in 6:10){
  acf(chain[,i], lag.max=100, main=para_b[i])
}

chain_elag_b = chain_burn[seq(1, nrow(chain_burn), by = 10), 6:10]
par(mfrow=c(2,5))
for (i in 6:10){
  plot(chain_elag_b[,i-5], type="l", xlab="Iterations", ylab="", main=para_b[i-5])
}
for (i in 6:10){
  acf(chain_elag_b[,i-5], lag.max=100, main="")
}
