

library(tidyverse)
library(nimble)
library(MCMCvis)

# https://groups.google.com/g/nimble-users/c/7y0ez80OYyo
# https://oliviergimenez.github.io/nimble-workshop/#107

# Generate data for testing
Nspp <- 5; Nsites <- 6; Nyears <- 4; Nocc <- 3

# Count Data
dat <- expand.grid(1:Nyears, 1:Nsites, 1:Nocc, 1:Nspp) %>% rename("year" = 1, "site" = 2, "survey" = 3, "spp" = 4) 
dat$count <- sample.int(20, nrow(dat), replace=T)
y <- array(unlist(dat[,5]), dim=c(Nspp, Nsites, Nyears, Nocc))


# SIMPLE MODEL, NO COVS ---------------------------------------------------

nimbleTest1 <- nimbleCode({
  
  ### PRIORS ###
  for(i in 1:Nspp){
    
    b0[i] ~ dnorm(0,0.01)
    p0[i] ~ dlogis(0,1)
  
  }
  
  ### PROCESS MODEL ###
  
  for(i in 1:Nspp){
    for(s in 1:Nsites){
      
      N[i,s,1] ~ dpois(lambda[i,s,1])                  # Initial abundance
      log(lambda[i,s,1]) <- b0[i]
      
      for(t in 2:Nyears){
      
        N[i,s,t] ~ dpois(lambda[i,s,t])                # Abundance in years t > 1
        log(lambda[i,s,t]) <- N[i,s,t-1] + b0[i]
      
      } # end t
      for(t in 1:Nyears){
        for(l in 1:Nocc){
          
          y[i,s,t,l] ~ dbin(det.p[i,s,t,l], N[i,s,t])  # Observation process
          logit(det.p[i,s,t,l]) <- p0[i]
        
        } # end l
      } # end t
    } # end s
  } # end i
}) # END MODEL

# Bundle data 
model.data <- list(y = y)
model.constants <- list(Nspp = Nspp, Nsites = Nsites, Nyears = Nyears, Nocc = Nocc)

# Initial values
ymax <- apply(y, 1:3, max)

inits <- list(N = round((ymax + 1)),
              b0 = rnorm(Nspp, 0, 1),
              p0 = runif(Nspp, 0, 1))

test_nim <- nimbleModel(code = nimbleTest1, constants = model.constants, data = model.data, inits = inits)

test_nim$initializeInfo()

#######################################################################

# MCMC settings
ni <- 1000; nt <- 3; nb <- 100; nc <- 3

# Parameters monitored
params <- c("N",
            "det.p",
            "b0",
            "p0"
) 

# Call JAGS from R 
test_nim <- nimbleMCMC(code = nimbleTest1,
                       constants = model.constants,
                       data = model.data,
                       inits = inits, 
                       monitors = params,  
                       thin  = nt, 
                       niter  = ni, 
                       nburnin  = nb)

sampSum <- samplesSummary(test_nim)




# MODEL V2 (ADD POP. GROWTH RATE) ----------------------------------------------------------------

nimbleTest2 <- nimbleCode({
  
  ### PRIORS ###
  for(i in 1:Nspp){
    
    b0[i] ~ dnorm(0,0.01)  # abundance intercept (unsure if necessary)
    p0[i] ~ dlogis(0,1)    # detection intercept (unsure if necessary)
    
    phi.mean[i] ~ dnorm(0,0.01)    # spp-specific mean pop. growth rate
    phi.alpha[i] ~ dnorm(0,0.01)   # intraspecific density effect
    
    for(j in 1:Nspp){
      
      phi.beta[i,j] ~ dnorm(0,0.01)  # Species interactions (not sure this is right)
      
    } # end j
  } # end i
  
  ### PROCESS MODEL ###
  
  ## POPULATION GROWTH RATE ##
  for(i in 1:Nspp){
    for(s in 1:Nsites){
      for(t in 2:Nyears){
        
        log(phi[i,s,t]) <- phi.mean[i] +                    # getting NAs for t=1, not sure if this is a problem
                           (phi.alpha[i] * N[i,s,t-1]) +
                           inprod(phi.beta[i,1:Nspp], N[1:Nspp,s,t-1])
        
      } # end t
    } # end s
  } # end i
  
  for(i in 1:Nspp){
    for(s in 1:Nsites){
      # Initial abundance
      N[i,s,1] ~ dpois(lambda[i,s,1])                  
      log(lambda[i,s,1]) <- b0[i]
      
      for(t in 2:Nyears){
        # Abundance in years t > 1
        N[i,s,t] ~ dpois(lambda[i,s,t])                
        log(lambda[i,s,t]) <- log(N[i,s,t-1] * phi[i,s,t]) + b0[i]
        
      } # end t
      for(t in 1:Nyears){
        for(l in 1:Nocc){
          # Observation process
          y[i,s,t,l] ~ dbin(det.p[i,s,t,l], N[i,s,t]) 
          logit(det.p[i,s,t,l]) <- p0[i]
          
        } # end l
      } # end t
    } # end s
  } # end i
}) # END MODEL

# Bundle data 
model.data <- list(y = y)
model.constants <- list(Nspp = Nspp, Nsites = Nsites, Nyears = Nyears, Nocc = Nocc)

# Initial values (not entirely confident on these values...)
ymax <- apply(y, 1:3, max)

inits <- list(N = round((ymax + 1)),
              b0 = rnorm(Nspp, 0, 1),
              p0 = runif(Nspp, 0, 1),
              phi.mean = rnorm(Nspp, 0, 1),
              phi.alpha = rnorm(Nspp, 0, 1),
              #phi.int = rnorm(Nspp, 0, 1),
              phi.beta = matrix(rnorm(Nspp*Nspp, 0, 1), nrow = Nspp, ncol = Nspp),
              phi = array(data = rnorm(Nspp*Nsites*Nyears), dim = c(Nspp, Nsites, Nyears)))

# Check model initializing
test_nim <- nimbleModel(code = nimbleTest2, 
                        constants = model.constants, 
                        data = model.data, 
                        inits = inits)

test_nim$initializeInfo()

# Run model 
# Parameters monitored
params <- c("N",
            "det.p",
            "b0",
            "p0",
            "phi",
            "phi.beta"
) 

test_nim2 <- nimbleMCMC(code = nimbleTest2,
                        constants = model.constants,
                        data = model.data,
                        inits = inits, 
                        monitors = params,  
                        thin  = 1, 
                        niter  = 1000, 
                        nburnin  = 100,
                        nchains = 1)

sampSum <- samplesSummary(test_nim2)

str(test_nim2)


# MODEL V3 (ADD POP. GROWTH RATE & MOVEMENT RATE) ----------------------------------------------------------------

nimbleTest3 <- nimbleCode({
  
  ### PRIORS ###
  for(i in 1:Nspp){
    
    b0[i] ~ dnorm(0,0.01)  # abundance intercept (unsure if necessary)
    p0[i] ~ dlogis(0,1)    # detection intercept (unsure if necessary)
    
    phi.mean[i] ~ dnorm(0,0.01)    # spp-specific mean pop. growth rate
    phi.alpha[i] ~ dnorm(0,0.01)   # intraspecific density effect
    
    kappa[i] ~ dnorm(0,0.01)
    
    for(j in 1:Nspp){
      
      phi.beta[i,j] ~ dnorm(0,0.01)  # Species interactions (not sure this is right)
      
    } # end j
    
    for(s in 1:Nsites){
      
      abund.init[i,s] ~ dunif(0,10)
      eta.init[i,s] ~ dnorm(0,0.01)
      
    }
    
  } # end i
  
  ### PROCESS MODEL ###
  
  ## POPULATION GROWTH RATE ##
  for(i in 1:Nspp){
    for(s in 1:Nsites){
      for(t in 2:Nyears){
        
        log(phi[i,s,t]) <- phi.mean[i] +                    # getting NAs for t=1, not sure if this is a problem
                           (phi.alpha[i] * N[i,s,t-1]) +
                           inprod(phi.beta[i,1:Nspp], N[1:Nspp,s,t-1])
        
      } # end t
    } # end s
  } # end i
  
  ## MOVEMENT RATE BETWEEN SITES ##
  # for(i in 1:Nspp){
  #   for(s in 1:Nsites){
  #     
  #     eta[i,s,1:Nsites] <- exp(-1 * kappa[i] * dist.mat[1:Nsites,s])
  #     theta[i,s,1:Nsites] <- eta[i,s,1:Nsites]/sum(eta[i,s,1:Nsites])
  #     
  #   } # end i
  # } # end s
  
  ## ABUNDANCE MODEL ##
  for(i in 1:Nspp){
    for(s in 1:Nsites){
      # Initial abundance
      N[i,s,1] ~ dpois(lambda[i,s,1])                  
      log(lambda[i,s,1]) <- b0[i] + abund.init[i,s] + eta.init[i,s]
      
      for(t in 2:Nyears){
        # Abundance in years t > 1
        N[i,s,t] ~ dpois(lambda[i,s,t])                
        log(lambda[i,s,t]) <- log(N[i,s,t-1] * 
                                  phi[i,s,t]) + 
                                  #theta[i,s,1:Nsites]) + # not sure how to correctly incorporate this... 
                                  b0[i]
        
      } # end t
      for(t in 1:Nyears){
        for(l in 1:Nocc){
          # Observation process
          y[i,s,t,l] ~ dbin(det.p[i,s,t,l], N[i,s,t]) 
          logit(det.p[i,s,t,l]) <- p0[i]
          
        } # end l
      } # end t
    } # end s
  } # end i
}) # END MODEL

# Site locations & distance matrix
lat <- runif(Nsites, 0, 5);lon <- runif(Nsites, 0, 5)
dist.mat <- matrix(NA, Nsites, Nsites)
for (i in 1:Nsites) {
  for (j in 1:Nsites) {
    dist.mat[i,j] <- sqrt((lat[i] - lat[j]) ^ 2 + (lon[i] - lon[j]) ^ 2)
  } # j
} # i

# Bundle data 
model.data <- list(y = y)

model.constants <- list(Nspp = Nspp, 
                        Nsites = Nsites, 
                        Nyears = Nyears, 
                        Nocc = Nocc,
                        dist.mat = dist.mat)

# Initial values (not entirely confident on these values...)
ymax <- apply(y, 1:3, max)

inits <- list(N = round((ymax + 1)),
              b0 = rnorm(Nspp,0,1),
              p0 = runif(Nspp,0,1),
              phi.mean = rnorm(Nspp,0,1),
              phi.alpha = rnorm(Nspp,0,1),
              #phi.int = rnorm(Nspp, 0, 1),
              phi.beta = matrix(rnorm(Nspp*Nspp,0,1), nrow = Nspp, ncol = Nspp),
              phi = array(data = rnorm(Nspp*Nsites*Nyears), dim = c(Nspp, Nsites, Nyears)),
              kappa = rnorm(Nspp,0,1),
              eta = array(data = runif(Nspp*Nsites*Nsites), dim=c(Nspp,Nsites,Nsites)),
              theta = array(data = runif(Nspp*Nsites*Nsites), dim=c(Nspp,Nsites,Nsites)),
              abund.init = matrix(data=runif(Nspp*Nsites,0,10), nrow = Nspp, ncol = Nsites),
              eta.init = matrix(data=rnorm(Nspp*Nsites,0,1), nrow = Nspp, ncol = Nsites)
              )

# Check model initializing
test_nim <- nimbleModel(code = nimbleTest3, 
                        constants = model.constants, 
                        data = model.data, 
                        inits = inits)

test_nim$initializeInfo()

# Run model 
# Parameters monitored
params <- c("N",
            "det.p",
            "b0",
            "p0",
            "phi",
            "phi.beta",
            "abund.init",
            "eta.init"
) 

test_nim3 <- nimbleMCMC(code = nimbleTest3,
                        constants = model.constants,
                        data = model.data,
                        inits = inits, 
                        monitors = params,  
                        thin  = 1, 
                        niter  = 500, 
                        nburnin  = 100,
                        nchains = 2)

sampSum <- samplesSummary(test_nim3)

test_nim$logProb_N

MCMCsummary(test_nim3, round=2, params = "abund.init", ISB=T, exact = F)
MCMCplot(test_nim3, params = "N")     
