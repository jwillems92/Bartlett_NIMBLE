
# Simulate data under Ricker model
sim.ricker <- function(lambda, r, K, p, nSites=100, nYears=45) {
    y <- N <- matrix(NA, nSites, nYears)
    N[,1] <- rpois(nSites, lambda)
    for(t in 2:nYears) {
        mu <- N[,t-1]*exp(r*(1-N[,t-1]/K))
        N[,t] <- rpois(nSites, mu)
    }
    y[] <- rbinom(nSites*nYears, N, p)
    return(list(N=N, y=y, lambda=lambda, r=r, p=p,
                seed=.Random.seed, seed.kind=RNGkind()))
}
