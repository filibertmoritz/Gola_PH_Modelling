### first steps in modelling with spOccupancy 
### script written by Filibert Heim, filibert.heim@posteo.de, in Feb 2025

### this script mainly follows the tutorial from Jeff Doser and Marc Kery (2022)
### https://doserlab.com/files/spoccupancy-web/articles/spacetimemodelshtml

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 1. Preparations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load packages 
library(tidyverse)
library(spOccupancy)
library(stars) # for spatial data 
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter

# set seed to be able to generate the same results as in the tutorial 
set.seed(1996)

# get the data and inspect its structure 
data("hbefTrends") # 9 year of presence/absence data for 12 species with 3 visits per year and 373 sites - which gives the dimensions
str(hbefTrends) # structured as a list containing 


# description of the hbfTrends object
## y - detection/non-detection data
hbefTrends$y %>% dim()

## occ.covs - occurence covariates used to explain occupancy (there are different structures in - yearlySiteCovs and siteCovs (in unmarked terminology)) 
hbefTrends$occ.covs$years # yearlySiteCovs in matrix structure 
str(hbefTrends$occ.covs$elev) # siteCovs in vector structure 

## det.covs # detection covariates 
hbefTrends$det.covs

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2. Data preparation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# subset data set to only one species 
revi.data <- hbefTrends
sp.names <- dimnames(hbefTrends$y)[[1]] # provide a name for species column 
revi.data$y <- revi.data$y[sp.names == 'REVI', , , ]
str(revi.data)
# str(hbefTrends) # compare to subsetted data 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 3. Initialise the model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# first have a look at raw occurence 
raw.occ.prob <- apply(revi.data$y, 2, mean, na.rm = T)
plot(x = 2010:2018, raw.occ.prob, type = 'b', 
     xlab = 'Year', ylab = 'Raw Occurence probability', main = 'Raw Occupancy Trajectory',
     ylim = c(0,1) )

# create the formula for the actual occupancy model 
revi.occ.formula <- ~ scale(years) + scale(elev) + I(scale(elev)^2) + (1 | years) + (1 | site.effect)
revi.det.formula <- ~ scale(day) + I(scale(day)^2) + scale(tod)

# specify values for model parameters 
z.inits <- apply(revi.data$y, c(1,2), function(a) as.numeric(sum(a, na.rm = T) > 0))
revi.inits <-list(beta = 0, # occurence coefs
                  alpha = 0, # det coefs
                  sigma.sq.psi = 1, # occ random effect variances 
                  z = z.inits) # latent occurence values 
revi.priors <- list(beta.normal = list(mean = 0, var = 2.72),
                    alpha.normal = list(mean = 0, var = 2.72), 
                    sigma.sq.psi.ig = list(a = 0.1, b = 0.1)) 

# set argument to control MCMK run 
n.chains <- 3
n.batch <- 200
batch.length = 25
n.samples <- n.batch * batch.length
n.burn = 2000
n.thin <- 12

# set ar1 to FALSE to indicate that we will fit model with an AR(1) random effect 
ar1 <- FALSE

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 4. Run model  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set up model to run, set up model with n.report = 50 to report progress after every 50th batch 
out <- tPGOcc(occ.formula = revi.occ.formula, 
              det.formula = revi.det.formula, 
              data = revi.data, 
              n.batch = n.batch, 
              batch.length = batch.length,
              inits = revi.inits,
              priors = revi.priors,
              ar1 = ar1,
              n.burn = n.burn, 
              n.thin = n.thin, 
              n.chains = n.chains, 
              n.report = 50)

# get summary of the model output 
summary(out)

# run the model with temporal autocorrelation parameter with default priors and 
out.ar1 <- tPGOcc(occ.formula = ~ scale(years) + scale(elev) + I(scale(elev)^2) + 
                    (1 | site.effect), 
                  det.formula = revi.det.formula, 
                  data = revi.data, 
                  n.batch = n.batch, 
                  batch.length = batch.length,
                  inits = revi.inits,
                  priors = revi.priors,
                  ar1 = TRUE,
                  n.burn = n.burn, 
                  n.thin = n.thin, 
                  n.chains = n.chains, 
                  n.report = 50)

# get the summary of the temoral autocorrelation model 
summary(out.ar1)

# compare models using WAIC 
waicOcc(out)
waicOcc(out.ar1)
# the WAIC values are fairly similar, thus, go for the simpler model


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 5. Assess GOF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# posterior predictive check on both models - normally not valid for binary response variables, thus group by data by sites (argument group = 1) and use the FT test
ppc.out <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1) 
summary(ppc.out)

# same for temporal random effect model 
ppc.out.ar1 <- ppcOcc(out.ar1, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out.ar1)

# shows appropriate fit since the overall bayesian p-value is close to 0.5 
# for some years the baz p-value is close to 1, suggesting, that our model creates more variability for some years than there is in the actual data for this time period


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 6. Prediction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get data on elevation 
data("hbefElev")
str(hbefElev)

# perform prediction 
predict(object = out.ar1, # model object 
        X.0 = , # design matrix of covariates at the prediction locations (3-dim array, site, prim period and covariate for pred), NB: if intercept is included in model, first cov should include 1s
        t.cols = , # vector of primary periods in design matrix,  we want to predict for, values should indicate cols in data$y
        ignore.RE = , # should we ignore unstructured random effects and just use fixed effects?
        type = 'occupancy') # choose between occupancy and detection

# prepare for prediction 
J.pred <- nrow(hbefElev) # number of prediction sites 
n.years.pred <- 2 # number of prediction years
p.occ <- ncol(out.ar1$beta.samples) # number of predictors 
elev.pred <- (hbefElev$val - mean(revi.data$occ.covs$elev)) / sd(revi.data$occ.covs$elev)
year.pred <- matrix(rep((c(2010, 2018) - mean(revi.data$occ.covs$years)) / 
                          sd(revi.data$occ.covs$years), 
                        length(elev.pred)), J.pred, n.years.pred, byrow = TRUE)
X.0 <- array(1, dim = c(J.pred, n.years.pred, p.occ)) # create 3 dim array
X.0[,,2] <- year.pred # years 
X.0[,,3] <- elev.pred # elev
X.0[,,4] <- elev.pred^2 # elev^2
str(X.0) # check array structure

# Indicate which primary time periods (years) we are predicting for
t.cols <- c(1, 9)
# Approx. run time: < 30 sec
out.pred <- predict(out.ar1, X.0, t.cols = t.cols, ignore.RE = TRUE, type = 'occupancy') # see explanation above!
# Check out the structure
str(out.pred)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 6. Päotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot.dat <- data.frame(x = hbefElev$Easting, 
                       y = hbefElev$Northing, 
                       mean.2009.psi = apply(out.pred$psi.0.samples[, , 1], 2, mean), 
                       mean.2018.psi = apply(out.pred$psi.0.samples[, , 2], 2, mean), 
                       sd.2009.psi = apply(out.pred$psi.0.samples[, , 1], 2, sd), 
                       sd.2018.psi = apply(out.pred$psi.0.samples[, , 2], 2, sd), 
                       stringsAsFactors = FALSE)
# Make a species distribution map showing the point estimates,
# or predictions (posterior means)
dat.stars <- st_as_stars(plot.dat, dims = c('x', 'y'))
# 2009
ggplot() + 
  geom_stars(data = dat.stars, aes(x = x, y = y, fill = mean.2009.psi)) +
  scale_fill_viridis_c(na.value = 'transparent') +
  labs(x = 'Easting', y = 'Northing', fill = '', 
       title = 'Mean REVI occurrence probability 2009') +
  theme_bw()


ggplot() + 
  geom_stars(data = dat.stars, aes(x = x, y = y, fill = mean.2018.psi)) +
  scale_fill_viridis_c(na.value = 'transparent') +
  labs(x = 'Easting', y = 'Northing', fill = '', 
       title = 'Mean REVI occurrence probability 2018') +
  theme_bw()

