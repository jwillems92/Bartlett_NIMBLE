

library(R2jags)
library(tidyverse)


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
model.data <- list(y = y,
                   Nspp = Nspp,
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

inits <- function() list(N = round((ymax + 1) / .5),
                         eta = eta.init,
                         theta = theta.init)

# Parameters monitored
params <- c("phi",
            "eta",
            "theta",
            "phi.beta",
            "lambda",
            "N",
            "det.p"
) 

# MCMC settings
ni <- 1000; nt <- 3; nb <- 100; nc <- 3

jagsTest1 <- jags(model.data, 
                  inits, 
                  params, 
                  "JAGS Version.txt", 
                  n.chains = nc, 
                  n.thin = nt, 
                  n.iter = ni, 
                  n.burnin = nb, 
                  working.directory = getwd())
