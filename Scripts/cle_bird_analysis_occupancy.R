
# load libraries
library(tidyverse)
library(cowplot) # to arrange ggplots into a grid
library(lubridate)
library(R2jags)
library(nimble)

# load full bird data
full_bird <- read_csv("BreedingBirdSurvey_Rem_corrected_19-02-2018.csv")

# load full covariate data  
cmp_cov <- read_csv("bird_cov.csv")

# correcting mis-identified species
# combine two versions of flickers
full_bird$species_id[full_bird$species_id == "Colaptes a. auratus"] <- "Colaptes auratus"

# scarlet-browed tanager to scarlet tanager
full_bird$species_id[full_bird$species_id == "Heterospingus xanthopygius"] <- "Piranga olivacea"

# spotted towhee to eastern towhee
full_bird$species_id[full_bird$species_id == "Pipilo maculatus/erythr."] <- "Pipilo erythrophthalmus"

# black-tailed gnatcatcher to blue gray gnatcatcher
full_bird$species_id[full_bird$species_id == "Polioptila melanura"] <- "Polioptila caerulea"

# bahama oriole to baltimore oriole
full_bird$species_id[full_bird$species_id == "Icterus northropi"] <- "Icterus galbula"

# add common names and order
common <- read_csv("bird_common_names_02.csv")
full_bird$common <- common$Common_names[match(full_bird$species_id,common$Scientific_name)]
full_bird$order <- common$Order[match(full_bird$species_id,common$Scientific_name)]

# remove the NAs - these are unidentified woodpeckers
full_bird <- filter(full_bird,!is.na(common ))


# remove irrelevant columns
colnames(full_bird)
# note col 24 and 25 are the same as col 1 - the survey id
full_bird <- full_bird[,-c(2,3,4,5,6,9,11,18,19,20,21,22,24,25,26,27,28,29,31)]

# check comments for anything that should be removed
#write.csv(full_bird,"full_bird.csv")

# one survey didnt get completed bc of rain
removals <- which(full_bird$site == "BE1015" &
          full_bird$date == "6/4/2017")
full_bird <- full_bird[-removals,]


# create detection history ------------------------------------------------

# number of visits per site
test <- full_bird %>%
  group_by(site) %>%
  summarize(visit = n_distinct(survey_fulcrum_id))

#as.factor(date) trick
# need to do this because tidyverse doesnt like posix format
# and doing it with date as character doesnt work bc 11 comes before 2
test_2 <- full_bird
test_2$date <- as.POSIXlt(test_2$date, format = "%m/%d/%Y")
test_2$date2 <- unlist(tapply(test_2$date,test_2$site,function(x){as.numeric(as.factor(x))}))

# check it
d <- data.frame(test_2$date,test_2$site, test_2$date2) # looks good

# this gives each survey a sequential number based upon the date it was surveyed
test <- full_bird %>%
  mutate(survey_no = test_2$date2) 

# change class of date because tidyverse is picky
test$date <- as.character(test$date)

# new df showing whether hoodeds within 50m are detected or not during each survey
hood <- test %>%
  group_by(site,survey_no) %>%
  summarize(det = ifelse(sum(common == "hooded warbler" & distance == "<=50m") == 1,1,0))

hood_wide <- pivot_wider(hood,
                              names_from = survey_no,
                              values_from = det)



# create site + survey covariates ---------------------------------------------

# sky as three level factor
sky <- test %>%
  select(site,survey_no,sky) 

# remove duplicated rows
sky <- sky[!duplicated(sky),]

# widen
sky <-  pivot_wider(sky, names_from = survey_no, 
              values_from = sky)
  
  

# wind as a two-level factor
wind <- test %>%
  select(site,survey_no,wind) 
wind <- wind[!duplicated(wind),]
wind <-  pivot_wider(wind, names_from = survey_no, 
                    values_from = wind)

# date of survey (weeks since start of surveying, as a factor)
surveydate <- test %>%
  select(site,survey_no,date) 
surveydate <- surveydate[!duplicated(surveydate),]
surveydate$date <- as.POSIXlt(surveydate$date, format = "%m/%d/%Y")
surveydate$date <- floor(yday(surveydate$date)/7) 
surveydate$date <- as.factor(surveydate$date - min(surveydate$date) + 1)# weeks since start of sampling
surveydate <-  pivot_wider(surveydate, names_from = survey_no, 
                     values_from = date)

# time of day survey started as hours since sunrise 
# technically hours since earliest survey, so hours since 30 min prior to sunrise
tod <- test %>%
  select(site,survey_no,time_started) 
# a few times were pm instead of am; fix
tod$time_started <- ifelse(hour(tod$time_started) > 11,
                           tod$time_started - 12*60*60,
                           tod$time_started)
tod <- tod[!duplicated(tod),]
tod$time_started <- as.numeric(tod$time_started/60/60) - min(as.numeric(tod$time_started/60/60))

# hours since sunrise as numeric
tod <-  pivot_wider(tod, names_from = survey_no, 
                           values_from = time_started)



# site-level covariates ---------------------------------------------------

# already have most of these in data.frame cmp_cov


# hood null occupancy model ---------------------------------------------------------

sink("hood")
cat("
    model {
    
    #######
    # Species-level priors
    psi ~ dbeta(1,1)  
    p ~ dbeta(1,1)   

    # Likelihood
    # Ecological model for true occurrence
    for (i in 1:R) {  # sites
    
    z[i] ~ dbern(psi)

    for (j in 1:T) { # weeks
    
    y[i,j] ~ dbern(p*z[i])
    
    } #j

    
    } #i
    occ_sites <- sum(z[])
    }
    ",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = hood_wide[,-1], R = nrow(hood_wide[,-1]), 
                 T = ncol(hood_wide[,-1]))

# Initial values
zst <- rep(1,nrow(hood_wide[,-1]))
inits <- function() list(z = zst)

# Parameters monitored
params <- c("psi",
            "p",
            "occ_sites"
) 

# MCMC settings
ni <- 1000; nt <- 3; nb <- 100; nc <- 3

# Call JAGS from R 
hood <- jags(win.data, inits, params, "hood", 
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(hood, dig = 2)


# same as above, but run as dbin
sink("hood")
cat("
    model {
    
    #######
    # Species-level priors
    psi ~ dbeta(1,1)  
    p ~ dbeta(1,1)   
    
    # Likelihood
    # Ecological model for true occurrence
    for (i in 1:R) {  # sites
    
    z[i] ~ dbern(psi)
    y[i] ~ dbin(p * z[i],N[i])
    } #i
    occ_sites <- sum(z[])
    }
    ",fill = TRUE)
sink()

# Bundle data
# create N vector = no. surveys per site
N <- apply(hood_wide[,-1],1,function(x){sum(!is.na(x))})

# create no. detections per row for data
dets <- apply(hood_wide[,-1],1,sum,na.rm=T)

win.data <- list(y = dets, R = nrow(hood_wide[,-1]),
                 N = N)

# Initial values
zst <- rep(1,nrow(hood_wide[,-1]))
inits <- function() list(z = zst)

# Parameters monitored
params <- c("psi",
            "p",
            "occ_sites"
) 

# MCMC settings
ni <- 1000; nt <- 3; nb <- 100; nc <- 3

# Call JAGS from R 
hood <- jags(win.data, inits, params, "hood", 
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(hood, dig = 2)




### JAGS VERSION ###
sink("hood")
cat("
    model {
    
    #######
    # Species-level priors
    psi ~ dbeta(1,1)  
    p ~ dbeta(1,1)   

    # Likelihood
    # Ecological model for true occurrence
    for (i in 1:R) {  # sites
    
    z[i] ~ dbern(psi)

    for (j in 1:T) { # weeks
    
    y[i,j] ~ dbern(p*z[i])
    
    } #j

    
    } #i
    occ_sites <- sum(z[])
    }
    ",fill = TRUE)
sink()

### NIMBLE VERSION ### 
null_hood <- nimbleCode({
  
  #######
  # Species-level priors
  psi ~ dbeta(1,1)  
  p ~ dbeta(1,1)   
  
  # Likelihood
  # Ecological model for true occurrence
  for (i in 1:R) {  # sites
    
    z[i] ~ dbern(psi)
    
    for (j in 1:T) { # weeks
      
      y[i,j] ~ dbern(p*z[i])
      
    } #j
    
  } #i
  occ_sites <- sum(z[1:98])
}
  
)

# Bundle data
win.data <- list(y = hood_wide[,-1], R = nrow(hood_wide[,-1]), 
                 T = ncol(hood_wide[,-1]))

# Initial values
zst <- rep(1,nrow(hood_wide[,-1]))
inits <- function() list(z = zst)

# Parameters monitored
params <- c("psi",
            "p",
            "occ_sites"
) 

# MCMC settings
ni <- 1000; nt <- 3; nb <- 100; nc <- 3

# Call JAGS from R 
hood_nim <- nimbleMCMC(code = null_hood,
                       constants = win.data, inits = inits, 
                       monitors = params,  
                       thin  = nt, niter  = ni, 
             nburnin  = nb, )
print(hood_nim, dig = 2)
summary(hood_nim)


# summary statistics -----------------------------------------------------------

# number of observations
nrow(full_bird) #4122

# number of sites
n_distinct(full_bird$site) # 98

# number of observers
n_distinct(full_bird$observer) # 13

# number of species
n_distinct(full_bird$species_id) # 90

# number of orders
n_distinct(full_bird$order) # 13

# number of surveys
n_distinct(full_bird$survey_fulcrum_id) # 353

# number of surveys per site mean and sd
f <- full_bird %>% 
  group_by(site) %>% 
  summarize(surveys = n_distinct(survey_fulcrum_id))
sd(f$surveys)
summary(f$surveys)

# number of species per survey
f <- full_bird %>% 
  group_by(survey_fulcrum_id) %>% 
  summarize(species = n_distinct(species_id))
sd(f$species)
summary(f$species)

# how many species detected 5 or fewer times?
d <- full_bird %>% 
  count(species_id)

d %>% 
  filter(n < 6) # 21

# how many species detected at least 100 times?
d %>% 
  filter(n > 99) # 15

# summary for distance
table(full_bird$distance)

# change to wider format to get summary for sky and wind conditions
t <- pivot_wider(full_bird, id_col = survey_fulcrum_id,
                 names_from = sky,
                 values_from = sky,
                 values_fn = list(sky = n_distinct))
sum(t$Cloudy, na.rm = T)
sum(t$`Partly Cloudy`, na.rm = T)
sum(t$Clear, na.rm = T)

# another way to do the above
t <- full_bird[,c(1,which(names(full_bird) == "sky"))]
# remove duplicated survey ids - i.e., get one row per survey
t <- t[-which(duplicated(t[,1])),]
table(t$sky)

# repeat for wind
t <- pivot_wider(full_bird, id_col = survey_fulcrum_id,
                 names_from = wind,
                 values_from = wind,
                 values_fn = list(wind = n_distinct))
sum(t$Light, na.rm = T)
sum(t$Calm, na.rm = T)

# summary of urbanization covariate
summary(cmp_cov$olm_imp_500m)
sd(cmp_cov$olm_imp_500m)

# exploratory visualizations --------------------------------------------------------

# number of observations by species for birds within 50 m
d <- full_bird %>%
  filter(distance == "<=50m") %>%
  group_by(common) %>%
  tally()

p1 <- ggplot(d,aes(x = n) ) +
  geom_histogram(aes(y = ..density..), color = "gray30",
                 fill = "white",
                 binwidth =  5) +
  geom_density(alpha = .2, fill = "antiquewhite3") +
  labs(x = "Number of detections per species",
       y = "Density") +
  theme_bw()+
  theme(text = element_text(size=16))
p1

# histogram of number of species recorded by site for birds < 50 m
d <- full_bird %>%
  filter(distance == "<=50m") %>%
  group_by(site) %>%
  summarize(site_richness = n_distinct(species_id)) 

p2 <- ggplot(data = d, aes(x = site_richness)) +
  geom_histogram(aes(y = ..density..), color = "gray30",
                 fill = "white",
                 binwidth =  1) +
  geom_density(alpha = .2, fill = "antiquewhite3") +
  labs(x = "Number of species per site",
       y = "Density") +
  theme_bw() +
  theme(text = element_text(size=16))
p2

# histogram of urban
p3 <- ggplot(cmp_cov, aes(x = olm_imp_500m)) +
  geom_histogram(aes(y = ..density..), color = "gray30",
                 fill = "white",
                 binwidth =  .05) +
  geom_density(alpha = .2, fill = "antiquewhite3") +
  labs(x = "Urbanization (proportion developed within 500m)",
       y = "Density") +
  theme_bw() +
  theme(text = element_text(size=15))
p3

plot_grid(p1, p2, p3, labels = c('A)', 'B)','C)'), label_size = 12)


# scatterplot of no. species by urban index
d

# check to make the bird data aligns with the covariate data
cmp_cov$camera_id == d$site

# add urbanization to the d dataframe
d$urban <-cmp_cov$olm_imp_500m

p4 <- ggplot(d, aes(x = urban, y = site_richness)) +
  geom_point()   +
  geom_smooth() +
  labs(x = "Urbanization \n(proportion developed within 500m)",
       y = "Number of species") +
  theme_bw() +
  theme(text = element_text(size=16))
p4

# urban - rural richness 
urb_rur <- kmeans(cmp_cov$olm_imp_500m,2)
cmp_cov$urb_rur <- recode(urb_rur$cluster,"1" = "urban", "2" = "rural")

d <- mutate(d, urb_rur = cmp_cov$urb_rur)

p5 <- ggplot(d, aes(x = urb_rur, y = site_richness, fill = urb_rur)) +
  geom_boxplot() + 
  labs(x = "Site type",
     y = "Number of species") +
  theme_bw() +
  theme(legend.position = "none")  +
  theme(text = element_text(size=16))
p5  

plot_grid(p4, p5, labels = c('A)', ' B)'), hjust = -0.05, label_size = 12)


# total number of species in urban sites

# first add urb_rur to full_bird
full_bird$urb_rur <- cmp_cov$urb_rur[match(full_bird$site,cmp_cov$camera_id)]

# get total no species by category
full_bird %>%
  filter(distance == "<=50m" &
           order == "Passeriformes" ) %>%
  group_by(urb_rur) %>%
  summarize(no_sp = n_distinct(species_id))


d <- full_bird %>%
  filter(distance == "<=50m" &
           order == "Passeriformes") %>%
  group_by(site,urb_rur) %>%
  summarize(site_richness = n_distinct(species_id)) 

# add gamma diversity as column
d$gamma <- ifelse(d$urb_rur == "urban",46,57)


# boxplot of beta diversity
ggplot(d, aes(x = urb_rur, y = gamma/site_richness)) +
  geom_boxplot(aes(fill = urb_rur)) +
  labs(x = "Site type",
       y = "Beta Diversity") +
  theme_bw()+
  theme(legend.position = "none") +
  theme(text = element_text(size=16)) +
  scale_x_discrete(labels = c('Rural','Urban'))

