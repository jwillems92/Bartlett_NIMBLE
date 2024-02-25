# Simulate data under Ricker + immigration model

# TO MODIFY
# (1) Change simulator... source(SOMETHING ELSE)
# (2) Change stub
# (3) Change simout
# (4) Change simulation parameters
# (5) Change pcountOpen arguments: "dynamics" and "K"
# (6) Results

stub <- "ricki_r0.05K10iota0.5"  # (2)
setwd(paste('/Users/hostetlerj/Documents/BirdClimate/Abundance/simulated', stub, sep='/'))
ls()
source("../sim_ricki.r")  # (1)
ls()

library(unmarked)

nsim <- 1000
simout.mle <- simout.cover <- matrix(NA, nsim, 5)  # (3)
colnames(simout.mle) <- colnames(simout.cover) <-
    c("Lambda", "r", "K", "iota", "p") # (3)

# Simulate
set.seed(345489)
lambda <- 10 # (4)
r <- 0.05    # (4)
K <- 10      # (4)
iota <- 0.5  # (4)
p <- 0.25    # (4)
for(i in 1:nsim) {
    cat("\ndoing", i, "at", format(Sys.time()), "\n")
    data.name <- paste("data_", stub, "_", i, sep="")
    data.file <- paste("data/", data.name, ".gzip", sep="")
    # Create simulated data
    sim.i <- sim.ricki(lambda=lambda, r=r, K=K, iota=iota, p=p, nSites=100, nYears=40) # (1)
    cat("   max(N)=", max(sim.i$N), "\n", sep="")
    assign(data.name, sim.i)
    # Save simulated data
    save(list=data.name, file=data.file)
    # Run pcountOpen model 
    umf.i <- unmarkedFramePCO(y=sim.i$y, numPrimary=ncol(sim.i$y))
    fm.i <- pcountOpen(~1, ~1, ~1, ~1, umf.i, K=200, dynamics="ricker", # (5)
      immigration=TRUE, starts=c(log(lambda), log(r), log(K), qlogis(p), log(iota))) # (5)
    # Save fitted models
    model.name <- paste("fm_", stub, "_", i, sep="")
    model.file <- paste("fm/", model.name, ".gzip", sep="")
    assign(model.name, fm.i)
    save(list=model.name, file=model.file) 
    mle.i <- coef(fm.i)
    simout.mle[i,] <- c(exp(mle.i[1]), exp(mle.i[2]), exp(mle.i[3]), exp(mle.i[5]), plogis(mle.i[4])) # (6)
    cat("   mle=", round(simout.mle[i,], 4), "\n")
    write.csv(simout.mle, paste(stub, "_mle.csv", sep=""),
              row.names=FALSE)
    CI.lam <- confint(fm.i, type="lambda") # (6)
    CI.r <- confint(fm.i, type="gamma")    # (6)
    CI.K <- confint(fm.i, type="omega")    # (6)
    CI.p <- confint(fm.i, type="det")      # (6)
    CI.iota <- confint(fm.i, type="iota")  # (6)
    simout.cover[i,1] <- log(lambda) >= CI.lam[1] & log(lambda) <= CI.lam[2] # (6)
    simout.cover[i,2] <- log(r) >= CI.r[1] & log(r) <= CI.r[2]               # (6)
    simout.cover[i,3] <- log(K) >= CI.K[1] & log(K) <= CI.K[2]               # (6)
    simout.cover[i,4] <- log(iota) >= CI.iota[1] & log(iota) <= CI.iota[2]   # (6)
    simout.cover[i,5] <- qlogis(p) >= CI.p[1] & qlogis(p) <= CI.p[2]         # (6)
    write.csv(simout.cover, paste(stub, "_cover.csv", sep=""),
              row.names=FALSE)
    # Clean up
    rm(list=c(data.name, model.name))
}


