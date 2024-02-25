


library(tidyverse)
library(nimble)
library(MCMCvis)



# DATA FORMATTING ---------------------------------------------------------
# Generate some count data for model testing
Nspp  <-  4
Nsites <- 3
Nyears <- 2
Nocc  <-  5

dat <- expand.grid(1:Nspp, 1:Nsites, 1:Nyears, 1:Nocc) %>% rename("spp"=1, "site"=2, "year"=3, "survey"=4)
dat$count <- rpois(nrow(dat),2)


y <- array(unlist(dat[,5]), dim=c(Nspp, Nsites, Nyears, Nocc))
dim(y)

# Generate abundance (abund.var1[s,t]), growth rate(phi.var1[s,t]), and detection (det.var1[i,s,t,m]) covariates
cov.dat <- expand.grid(1:Nsites, 1:Nyears) %>% rename("site" = 1, "year" = 2)
cov.dat$abund.var1 <- runif(nrow(cov.dat), 20, 50)
cov.dat$phi.var1 <- runif(nrow(cov.dat), 2, 8)

abund.var1 <- cov.dat %>% select(-c(phi.var1)) %>% 
  pivot_wider(names_from = "year", values_from = abund.var1) %>% 
  column_to_rownames(var="site") %>% 
  as.matrix()

phi.var1 <- cov.dat %>% select(-c(abund.var1)) %>% 
  pivot_wider(names_from = "year", values_from = phi.var1) %>% 
  column_to_rownames(var="site") %>% 
  as.matrix()

det.dat <- expand.grid(1:Nspp, 1:Nsites, 1:Nyears, 1:Nocc) %>% rename("spp"=1, "site"=2, "year"=3, "survey"=4)
det.dat$det.var1 <- runif(nrow(det.dat),2,10)
det.dat.array <- array(unlist(det.dat[,5]), dim=c(Nspp, Nsites, Nyears, Nocc))

# WRITE THE MODEL ---------------------------------------------------------

nimbleReprex1 <- nimbleCode({
##############  
### PRIORS ###
##############
  
for(i in 1:Nspp){
  
  phi.mu[i] ~ dnorm(0,0.01)
  alpha[i] ~ dnorm(0,0.01)
  phi.beta[i] ~ dnorm(0,0.01)
  
  lambda.mu[i] ~ dunif(0,5)
  
  kappa[i] ~ dnorm(0,0.01)
  
  b.abund[i] ~ dnorm(0,0.01)
  mu.det.p[i] ~ dunif(0,1)
  b.det[i] ~ dlogis(0,1)
  
  for(s in 1:Nsites){
    
    lambda.eta[i,s] ~ dnorm(0,0.01)
    
    for(t in 1:Nyears){
      
      eta.abund[i,s,t] ~ dnorm(0,0.01)
      #lambda[i,s,t] ~ dunif(0,5)
      #phi[i,s,t] ~ dnorm(0,0.01)
      
      for(m in 1:Nocc){
        
        eta.det[i,s,t,m] ~ dlogis(0,1)
        #det.p[i,s,t,m] ~ dunif(0,1)
        
      } # end m
    } # end t
  } # end s
} # end i
  
for(i in 1:Nspp){
  for(j in 1:Nspp){
    
    beta[i,j] ~ dnorm(0,0.01)
    
  } # end j
} # end i  

###################  
## PROCESS MODEL ##
###################

## LOCAL POPULATION GROWTH ##
for(i in 1:Nspp){
  for(s in 1:Nsites){
    
    log(phi[i,s,1]) ~ dnorm(0,0.01)
    
    for(t in 2:Nyears){
      
      log(phi[i,s,t]) <- log(phi.mu[i]) +                            # mean local pop. growth rate for each spp
                         (alpha[i] * N[i,s,t-1]) +                   # effect of intra-spp. density dependence on growth rate
                         inprod(beta[i,1:Nspp], N[1:Nspp,s,t-1]) +   # effect of spp j on pop. growth rate of spp. i
                         phi.beta[i] * phi.var1[s,t]                 # effect of environmental cov(s) on growth rate  
      
    } # end t
  } # end s
} # end i
  
## MOVEMENT RATE BETWEEN SITES ##
for(i in 1:Nspp){
  for(s in 1:Nsites){
    
    eta.theta[i,s,1:Nsites] <- exp(-1 * kappa[i] * dist.mat[1:Nsites,s])
    theta[i,s,1:Nsites] <- eta.theta[i,s,1:Nsites]/sum(eta.theta[i,s,1:Nsites])
    
  } # end i
} # end s
  
  
## ABUNDANCE MODEL ##
for(i in 1:Nspp){
  for(s in 1:Nsites){
    
    N[i,s,1] ~ dpois(lambda[i,s,1])  # abundance in first primary observation occasion
    
    log(lambda[i,s,1]) <- log(lambda.mu[i]) +                # species-specific mean initial abundance
                          b.abund[i] * abund.var1[s,1] +     # effect of environmental cov(s) on initial abundance
                          lambda.eta[i,1]                    # error term 
    
    for(t in 2:Nyears){
      
      N[i,s,t] ~ dpois(lambda[i,s,t])  # abundance in subsequent primary observation occasions
      
      log(lambda[i,s,t]) <- log(sum(N[i,s,t-1] * phi[i,s,t] * theta[i,s,1:Nsites])) + eta.abund[i,s,t]

    } # end t
    
## OBSERVATION MODEL ## 
    for(t in 1:Nyears){
      for(m in 1:Nocc){
      
        y[i,s,t,m] ~ dbinom(det.p[i,s,t,m], N[i,s,t])
        
        logit(det.p[i,s,t,m]) <- logit(mu.det.p[i]) +            # mean local detection prob. for each spp
                                 b.det[i] * det.var1[i,s,t,m] +  # effect of environmental cov(s) on det. prob.
                                 eta.det[i,s,t,m]                # error term
        
      } # end m
    } # end t
  } # end s
} # end i
}) # END MODEL

# COMPILE DATA FOR RUNNING MODEL ------------------------------------------

# Site locations & distance matrix
lat <- runif(Nsites, 0, 5);lon <- runif(Nsites, 0, 5)
dist.mat <- matrix(NA, Nsites, Nsites)
for (i in 1:Nsites) {
  for (j in 1:Nsites) {
    dist.mat[i,j] <- sqrt((lat[i] - lat[j]) ^ 2 + (lon[i] - lon[j]) ^ 2)
  } # j
} # i

# Bundle data 
model.data <- list(y = y,
                   abund.var1 = abund.var1,
                   phi.var1 = phi.var1,
                   det.var1 = det.dat.array,
                   dist.mat = dist.mat)

model.constants <- list(Nspp = Nspp, 
                        Nsites = Nsites, 
                        Nyears = Nyears, 
                        Nocc = Nocc)

# Initial values (not entirely confident on these values...)
ymax <- apply(y, 1:3, max)

inits <- list(N =          round((ymax + 1)),                         # abundance
              lambda.mu =  runif(Nspp,0,5),                           # initial mean abundance  
              lambda =     array(runif(Nspp*Nsites*Nyears,0,10),      # abundance
                                 dim = c(Nspp, Nsites, Nyears)),
              phi.mu =     runif(Nspp,0,5),                           # mean pop. growth rate  
              alpha =      rnorm(Nspp,0,0.01),                        # intra-spp density dependence
              phi.beta =   rnorm(Nspp,0,0.01),                        # effect of growth variable 1
              kappa =      runif(Nspp,0,1),                           # decay param. for site movement (constrain to be positive?)
              b.abund =    rnorm(Nspp,0,0.01),                        # effect of abund variable 1
              mu.det.p =   runif(Nspp,0,1),                           # mean detection prob
              b.det =      rnorm(Nspp,0,0.01),                        # effect of detection variable 1
              lambda.eta = matrix(rnorm(Nspp*Nsites,0,1),             # error for initial abundance
                                  nrow = Nspp),
              eta.abund =  array(rnorm(Nspp*Nsites*Nyears),           # error for abundance
                                 dim = c(Nspp, Nsites, Nyears)),
              det.p   =    array(runif(Nspp*Nsites*Nyears*Nocc,0,1),  # detection probability
                                 dim = c(Nspp, Nsites, Nyears,Nocc)),
              eta.det =    array(rnorm(Nspp*Nsites*Nyears*Nocc),      # error for detection
                                 dim = c(Nspp, Nsites, Nyears,Nocc)),
              beta =       matrix(rnorm(Nspp*Nspp,0,0.01),            # inter-spp interactions
                                  nrow = Nspp),
              phi =        array(rnorm(Nspp*Nsites*Nyears,0,0.01),    # growth rate
                                 dim = c(Nspp, Nsites, Nyears)),
              theta =      array(rnorm(Nspp*Nsites*Nsites,0,0.01),    # movement rate (constrain to be positive?)
                                 dim = c(Nspp, Nsites, Nsites))
              )

# Run model 
# Parameters monitored
params <- c("N",
            "det.p",
            "phi",
            "theta",
            "kappa",
            "beta"
) 

out1 <- nimbleMCMC(code = nimbleReprex1,
                   data = model.data,
                   inits  = inits,
                   constants = model.constants,
                   monitors = params,
                   niter = 1000,
                   nburnin = 200,
                   nchains = 2)

head(out1)
samps <- samplesSummary(out1$chain1)

samps.df <- as.data.frame(samps) %>%  
  rownames_to_column(var="param") %>% 
  filter(str_detect(param, "N"))

samps.df

library(MCMCvis)
out1.sum <- MCMCsummary(out1, round=2)

MCMCplot(out1, params = "N")

MCMCtrace(object = out1,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          Rhat = T,
          n.eff = T,
          params = "beta")





# Check model initializing
befNim1 <- nimbleModel(code = nimbleReprex1, 
                       constants = model.constants, 
                       data = model.data, 
                       inits = inits)

befNim1$initializeInfo()
log(befNim1$calculate("phi.mu"))


# Compile the model
CbefNim1 <- compileNimble(befNim1)

# Configure MCMC
befNim1_Conf <- configureMCMC(befNim1)

# Create MCMC function
befNim1_MCMC <- buildMCMC(befNim1_Conf,
                          monitors = befNim1$getNodeNames(stochOnly = F, includeData = F))
CbefNim1_MCMC <- compileNimble(befNim1_MCMC, 
                               project = befNim1)



samples <- runMCMC(mcmc = CbefNim1_MCMC,
                   niter = 500,
                   nburnin = 100)









Nc1 <- samplesSummary(out1$chain1) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="param") %>% 
  filter(str_detect(param,"N"))

Nc2 <- samplesSummary(out1$chain2) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="param") %>% 
  filter(str_detect(param,"N"))

Nc3 <- samplesSummary(out1$chain3) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="param") %>% 
  filter(str_detect(param,"N"))

rowMeans(cbind(Nc1$Mean, Nc2$Mean, Nc3$Mean))



