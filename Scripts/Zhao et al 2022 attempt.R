

library(tidyverse)
library(nimble)
library(AHMbook)


nimbleTest <- nimbleCode({
  
  ### PRIORS ###
  
  for(i in 1:Nspp){
    
    abund.init[i] ~ dnorm(0,0.01)
    sigma.init[i] ~ dnorm(0,0.01)
    b.int[i]      ~ dnorm(0,0.01)
    
    phi.mean[i]  ~ dnorm(0,0.01)
    phi.alpha[i] ~ dnorm(0,0.01)
    phi.int[i]   ~ dnorm(0,0.01)
    
    kappa[i] ~ dnorm(0,0.01)
    
    det.p.mean[i] ~ dnorm(0,0.01) 
    sigma.det[i]  ~ dnorm(0,0.01)
    p.int[i]      ~ dnorm(0,0.01) 
    
    for(s in 1:Nsites){
      
      eta.init[i,s] ~ dnorm(0, sigma.init[i])
      
      for(t in 1:Nyears){
        
        eta.abund[i,s,t] ~ dnorm(0,0.01) 
        
        for(l in 1:Nocc){
          
          eta.obs[i,s,t,l] ~ dnorm(0, sigma.det[i])
          
        } # end l
      } # end t
    } # end s
  } # end i
  
  
  ## SPECIES INTERACTIONS ##
  
  for(i in 1:Nspp){
    for(j in 1:Nspp){
      
      if(i != j){phi.beta[i,j] ~ dnorm(0,0.01)}
      else{phi.beta[i,j] <- 0}
      
    } # end j
  } # end i
  
### PROCESS MODEL ###
  
  ## POPULATION GROWTH RATE ##
  
  for(i in 1:Nspp){
    for(s in 1:Nsites){
      for(t in 2:Nyears){
        
        log(phi[i,s,t]) <- phi.mean[i] + 
                           (phi.alpha[i] * N[i,s,t-1]) + 
                           inprod(phi.beta[i,1:Nspp], N[1:Nspp,s,t-1]) +   # not confident on this...
                           #(phi.beta[i,1:5] %*% N[1:5,s,t-1]) +
                           phi.int[i]
        
      } # end t
    } # end s
  } # end i
  
  
  ## MOVEMENT RATE BETWEEN SITES ##
  
  for(i in 1:Nspp){
    for(s in 1:Nsites){
      
      eta[i,s,1:Nsites] <- exp(-1 * kappa[i] * dist[1:Nsites,s])
      theta[i,s,1:Nsites] <- eta[i,s,1:Nsites]/sum(eta[i,s,1:Nsites])
      
    } # end s
  } # end i
  
  
  ## ABUNDANCE ##
  
  for(i in 1:Nspp){
    for(s in 1:Nsites){
      # Initial Abundance in Year 1
      N[i,s,1] ~ dpois(lambda[i,s,1])
      
      log(lambda[i,s,1]) <- abund.init[i] + b.int[i] + eta.init[i,s]
        
      for(t in 2:Nyears){
        # Abundance in Subsequent Years  
        N[i,s,t] ~ dpois(lambda[i,s,t])
        
        log(lambda[i,s,t]) <- log(N[i,s,t-1] * phi[i,s,t] * theta[i,s,t]) + eta.abund[i,s,t]
          
        for(l in 1:Nocc){
          # Observation Process  
          y[i,s,t,l] ~ dbin(N[i,s,t], det.p[i,s,t,l])
          
          logit(det.p[i,s,t,l]) <- det.p.mean[i] + p.int[i] + eta.obs[i,s,t,l]
            
        } # end l
      } # end t
    } # end s
  } # end i
})


# Generate data for testing
Nspp <- 5;Nsites <- 6;Nyears <- 4;Nocc <- 3

# Site locations & distance matrix
lat <- runif(Nsites, 0, 5);lon <- runif(Nsites, 0, 5)

dist.mat <- matrix(NA, Nsites, Nsites)
for (i in 1:Nsites) {
  for (j in 1:Nsites) {
    dist.mat[i,j] <- sqrt((lat[i] - lat[j]) ^ 2 + (lon[i] - lon[j]) ^ 2)
  } # j
} # i

# Count Data
dat <- expand.grid(1:Nyears, 1:Nsites, 1:Nocc, 1:Nspp) %>% rename("year" = 1, "site" = 2, "survey" = 3, "spp" = 4) 
dat$count <- sample.int(20, nrow(dat), replace=T)
y <- array(unlist(dat[,5]), dim=c(Nspp, Nsites, Nyears, Nocc))

# Bundle data 
model.data <- list(y = y)

model.constants <- list(Nspp = Nspp,
                        Nsites = Nsites,
                        Nyears = Nyears,
                        Nocc = Nocc,
                        dist = dist.mat)

# Initial values
ymax <- apply(y, 1:3, max)

eta.init <- array(NA, dim=c(Nspp, Nsites, Nsites))
theta.init <- array(NA, dim=c(Nspp, Nsites, Nsites))

for(i in 1:Nspp){
  for(s in 1:Nsites){
    
    eta.init[i,s,] <- exp(-1 * 1 * dist.mat[,s])
    theta.init[i,s,] <- eta.init[i,s,]/sum(eta.init[i,s,])
    
  } # end s
} # end i

inits <- list(N = round((ymax + 1) / .5),
              phi.beta = matrix(0,Nspp, Nspp),
              kappa = rep(1, Nspp),
              phi = array(1, dim=c(Nspp, Nsites, Nyears)),
              eta = eta.init,
              theta = theta.init)

# Parameters monitored
params <- c("phi",
            "eta",
            "theta",
            "phi.beta",
            "lambda",
            "N",
            "det.p",
            "y"
) 

# MCMC settings
ni <- 1000; nt <- 3; nb <- 100; nc <- 3

# Call JAGS from R 
test_nim <- nimbleMCMC(code = nimbleTest,
                       constants = model.constants,
                       data = model.data,
                       inits = inits, 
                       monitors = params,  
                       thin  = nt, 
                       niter  = ni, 
                       nburnin  = nb)

 sampSum <- samplesSummary(test_nim)






# Parameters to monitor
params <- c('lambda', 'phi', 'gamma', 'p', 'Nbar')

Rmodel <- nimbleModel(code=code, constants=constants, data=bdata, inits=inits,
                      calculate=FALSE, check=FALSE)

Cmodel <- compileNimble(Rmodel)
mcmcspec<-configureMCMC(Rmodel, monitors=params, thin=1)
mcmcspec$removeSamplers(c("phi","gamma"))
mcmcspec$addSampler(target=c("phi","gamma"), type="AF_slice")

scrMCMC <- buildMCMC(mcmcspec)

CSCRMCMC <- compileNimble(scrMCMC, project = Rmodel)

nb<-1000 
ni<-5000 +nb 
nc<-1
outNim <- runMCMC(CSCRMCMC, niter = ni , nburnin = nb ,
                  nchains = nc,
                  setSeed = FALSE, progressBar = TRUE,
                  samplesAsCodaMCMC = TRUE)

summary(outNim)

# warning: problem initializing stochastic node N[5, 6, 3]: logProb is -Inf.
# warning: logProb of data node y[1, 1, 3, 1]: logProb is NA or NaN.
# Warning: slice sampler reached maximum number of contractions






## build the model
thinModel <- nimbleModel(nimbleTest,
                         data = model.data,
                         constants = model.constants,
                         inits = inits)

## build the MCMC
thinMCMC <- buildMCMC(thinModel)

## run uncompiled to have R errors and debugging available
thinMCMC$run(1)

## I get:
## warning: problem initializing stochastic node var_epsilon[1], logProb is NA or NaN
## Warning message:
## In rnorm(1, mean = (prior_mean * prior_tau + contribution_mean)/(prior_tau +  :
##  NAs produced

## This suggests a problem triggered by something wrong with var_epsilon, so I will look at its value in the model:

thinModel$det.p


## This reveals that var_epsilon contains one Inf and two crazy values, so I suggest putting reasonable values in inits.