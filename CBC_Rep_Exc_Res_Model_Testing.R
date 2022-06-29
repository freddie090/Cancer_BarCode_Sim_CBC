
# Cancer Barcode Simulation
# Freddie Whiting - 2022

# Developing models for the 'CBC-Rep_Exc_Res' (replicate, exclusive resistance) 
# inference. 

# Now I have the theory that also leverages the variance in the distributions, 
# slowly increasing the complexity of the models until I have one thats 
# working for the full 'assay reading' version. 

################################################################################

# Read in source files and simulation outputs. 

rm(list = ls())

source("C:/Users/fwhiting/My Drive/Barcode_Simulations/Cancer_BarCode_Sim_CBC/CBC_Analysis_Custom_Functions.R")
#source("C:/Users/whitin01/Google Drive/Barcode_Simulations/Cancer_BarCode_Sim_CBC/CBC_Analysis_Custom_Functions.R")
#source("C:/Users/Freddie/Google Drive/Barcode_Simulations/Cancer_BarCode_Sim_CBC/CBC_Analysis_Custom_Functions.R")

setwd("C:/Users/fwhiting/My Drive/Barcode_Simulations/Cancer_BarCode_Sim_CBC/Outputs/Rep_Exc_Res/")
#setwd("C:/Users/whitin01/Google Drive/Barcode_Simulations/Cancer_BarCode_Sim_CBC/Outputs/Rep_Exc_Res/")
#setwd("C:/Users/Freddie/Google Drive/Barcode_Simulations/Cancer_BarCode_Sim_CBC/Outputs/Rep_Exc_Res/")

out_dirs <- grep("CBC_Rep_Exc_Res_Exp", list.files(), value = T)

# run through output directories and save outputs. 

res_dfs <- lapply(out_dirs, function(x){
  setwd(x)
  out_df <- read.csv("res_df.csv", stringsAsFactors = F)
  setwd("../")
  return(out_df)
})

exp_res_dfs <- lapply(out_dirs, function(x){
  setwd(x)
  out_df <- read.csv("exp_res_df.csv", stringsAsFactors = F)
  setwd("../")
  return(out_df)
})

res_df <- bind_rows(res_dfs)
exp_res_df <- bind_rows(exp_res_dfs)

res_df$mu <- as.factor(res_df$mu)
res_df$sig <- as.factor(res_df$sig)

exp_res_df$mu <- as.factor(exp_res_df$mu)
exp_res_df$sig <- as.factor(exp_res_df$sig)

################################################################################

# Change to the STAN directory and load STAN libraries. 

library(rstan)
library(bayesplot)

options(mc.cores = parallel::detectCores())

rstan_options(auto_write = TRUE)

setwd("../../STAN_Models/")

################################################################################

# Choose a parameter pair to test the inference on:

res_df_sub <- subset(res_df, mu == 1e-05 & sig == 1e-01 & Nsim == 1)
res_df_sub <- res_df_sub[, !(names(res_df_sub) %in% "ass_read")]

################################################################################

# Read in the theory functions for testing. 


################################
# A function that is the solution to the proportion of resistance
# at time t given the initial condition p = 1.0 (where p = nR / (nr  + nR)) = 
################################

res_prop_t <- function(mu, sig, del, br, dr, t){
  
  if(del == 0){
    
    p_t <- (1 - (mu/(mu + sig)))*exp((-2 * mu * br * t) - (2 * sig * br * t)) + mu/(mu + sig)
    
  } else {
    
    lam_r <- (br - dr)
    lam_R <- lam_r - (lam_r * del)
    bR <- br - (br * del)
    
    a <- (-lam_R + lam_r)
    b <- lam_R - (2*mu*br) - (2*sig*bR) - lam_r
    c <- 2*mu*br
    
    K <- 2/(sqrt(b^2 - (4*a*c))) * atanh((2*a + b)/sqrt(b^2 - (4*a*c)))
    
    p_t <- ((sqrt(b^2 - (4*a*c))) / (2*a)) * tanh(((K - t)/2) * (sqrt(b^2 - (4*a*c)))) - (b/(2*a))
    
  }
  return(p_t)
}

res_prop_t <- Vectorize(res_prop_t)

################################


################################
# Simulate an assay reading. 
################################

ass_read <- function(x1, x2, mes_err, b, d, mu, sig, del, p){
  eqR <- stable_R_pheno(mu = mu, sig = sig, del = del, b = b, d = d)
  m_est <- (x2 - x1)/(1 - eqR)
  c_est <- x2 - (m_est * 1.0)
  x_read <- rnorm(1, mean = ((m_est*p) + c_est), sd = mes_err)
  return(x_read)
}

ass_read <- Vectorize(ass_read)

################################

################################
# A function to return the assay reading expectation: 
################################

ass_exp <- function(x1, x2, mu, sig, del, b, d, exp_pRt){
  
  eq_R <- stable_R_pheno(mu, sig, del, b, d)
  m <- (x2 - x1)/(1 - eq_R)
  c <- x2 - (1 * m)
  exp_x <- (exp_pRt * m) + c
  return(exp_x)
}

ass_exp <- Vectorize(ass_exp)

################################

################################
# A function to return the assay reading variance: 
################################

ass_var <- function(x1, x2, mu, sig, del, b, d, var_pRt){
  
  eq_R <- stable_R_pheno(mu, sig, del, b, d)
  m <- (x2 - x1)/(1 - eq_R)
  c <- x2 - (1 * m)
  var_x <- m^2 * var_pRt
  return(var_x)
}

ass_var <- Vectorize(ass_var)

################################


################################
# A function to return the equilibrium frequency of resistance: 
################################

res_eq_freq<-function(mu, sigma){
  
  eq_R <- mu/(mu + sigma)
  
  return(eq_R)
  
}

res_eq_freq <- Vectorize(res_eq_freq)

################################


################################
# A function to return the expected proportion of resistance at time t: 
################################

exp_res_prop_t<-function(b, d, del_t, mu, sigma){
  
  exp_p_t <- (1 - (mu/(mu + sigma)))*exp((-2 * mu * b * del_t) - (2 * sigma * b * del_t)) + mu/(mu + sigma)
  
  return(exp_p_t)
  
}

exp_res_prop_t <- Vectorize(exp_res_prop_t)

################################


################################
# A function to return the variance of the proportion of resistance at time t: 
################################

var_res_prop_t<-function(b, d, del_t, mu, sigma, Ni){
  
  exp_p_t <- exp_res_prop_t(b, d, del_t, mu, sigma)
  var_p_t <- ((exp_p_t) * (1 - exp_p_t))/Ni
  
  return(var_p_t)
}

var_res_prop_t <- Vectorize(var_res_prop_t)

################################


################################
# A function to return the expected assay reading value at time t: 
################################


exp_x_read_t<-function(x1_est, x2_est, mu, sigma, exp_pR){
  
  eq_R <- res_eq_freq(mu, sigma)
  m_est <- (x2_est - x1_est)/(1 - eq_R)
  c_est <- x2_est - (m_est * 1.0)
  exp_x_t <- (m_est*exp_pR) + c_est
  
  return(exp_x_t)
  
}

exp_x_read_t <- Vectorize(exp_x_read_t)

################################


################################
# A function to return the variance of an assay reading value at time t: 
################################

var_x_read_t<-function(x1_est, x2_est, mu, sigma, var_pR){
  
  eq_R <- res_eq_freq(mu, sigma);
  m_est <- (x2_est - x1_est)/(1 - eq_R);
  c_est <- x2_est - (m_est * 1.0);
  var_x_t <- m_est^2 * var_pR;
  
  return(var_x_t)
  
}

var_x_read_t <- Vectorize(var_x_read_t)

################################


# Just check that the functions written originally, and the ones 'lifted' from 
# the STAN models work identically: 

res_prop_t(1e-05, 1e-01, 0.000, 0.893, 0.200, 6.00)
exp_res_prop_t(0.893, 0.200, 6.00, 1e-05, 1e-01)
# These do. 

ass_exp(1.0, 10.0, 1e-05, 1e-01, 0.0, 0.893, 0.200, 0.3424892)
exp_x_read_t(1.0, 10.0, 1e-05, 1e-01, 0.3424892)
# These do. 

ass_var(1.0, 10.0, 1e-05, 1e-01, 0.0, 0.893, 0.200, 0.3424892)
var_x_read_t(1.0, 10.0, 1e-05, 1e-01, 0.3424892)
# These do. 

################################################################################

# To begin with, we will simulate the assay readings with some given values of
# x1 and x2, and with no additional reading variance. 

# Additional error noise. 
mes_err <- 0.00 

# The number of replicates per reading - same for x1 and x2 as the simulated 
# output readings. 
K <- unique(res_df_sub$K) 

# The number of cells in each replicate recording.
Ni <- unique(res_df_sub$Ni) 

# A vector of recording times, and the number of recording times. 
del_ts <- unique(res_df_sub$delt)
nT <- length(del_ts)

# Birth and death rates from the simulation. 
b <- unique(res_df_sub$b)
d <- unique(res_df_sub$d)

# Choose values for x1 and x2 to normalise assay readings.  
x2 <- 10
x1 <- 1

# Assume that readings ~Normal(mean = true_value, sd = mes_err).
x1_reads <- rnorm(K, mean = x1, sd = mes_err)
x2_reads <- rnorm(K, mean = x2, sd = mes_err)

# Simulate the assay readings 
res_df_sub["ass_read"] <- ass_read(
  rep(x1, nrow(res_df_sub)), 
  rep(x2, nrow(res_df_sub)), 
  rep(mes_err, nrow(res_df_sub)),
  res_df_sub$b, res_df_sub$d, anac(res_df_sub$mu), 
  anac(res_df_sub$sig), anac(res_df_sub$del), 
  res_df_sub$res_prop
)

# Extract the simulated proportions of resistance:

pR_del_t1 <- subset(res_df_sub, delt == del_ts[[1]])$res_prop
pR_del_t2 <- subset(res_df_sub, delt == del_ts[[2]])$res_prop
pR_del_t3 <- subset(res_df_sub, delt == del_ts[[3]])$res_prop

# Extract simulated assay readings 

x_del_t1 <- subset(res_df_sub, delt == del_ts[[1]])$ass_read
x_del_t2 <- subset(res_df_sub, delt == del_ts[[2]])$ass_read
x_del_t3 <- subset(res_df_sub, delt == del_ts[[3]])$ass_read


################################################################################

# Model 1. - Just use the mean of the proportion of resistance to infer mu and 
# sigma (don't include the variance term yet). 
# Fits the total variance for each separately as a free parameter. 

# Model data:
model_dat <- list(b = b, d = d, nT = nT, del_ts = del_ts, K = K,
                  pR_del_t1 = pR_del_t1, pR_del_t2 = pR_del_t2, 
                  pR_del_t3 = pR_del_t3)

# Model fit:
model_fit <- stan(file = "CBC_Rep_Exc_Res_Model1.stan",
                  data = model_dat, chains = 4, iter = 2000)

# Output fits:
mcmc_dens_overlay(model_fit, pars = c("mu", "sigma"))
mcmc_pairs(model_fit, pars = c("mu", "sigma"))
mcmc_trace(model_fit, pars = c("mu", "sigma"))

# Now, we can also look at the fitted values of additional variance (epsi_i) 
# and see how close these are to the theoretical variances (that don't include
# the birth-death variance). 

mcmc_dens_overlay(model_fit, pars = c("epsi_1", "epsi_2", "epsi_3"))
var_res_prop_t(b, d, del_ts, 1e-05, 1e-01, Ni)

# Think i've realised problem before was i wasn't accounting for the fact
# that the STAN normal function takes the s.d., opposed to variance. 

fit_df <- data.frame(as.array(model_fit))
# Get the theoretical expectations: 
epsi_df <- data.frame(t_epsi = sqrt(var_res_prop_t(b, d, del_ts, 1e-05, 1e-01, Ni)),
                      var = c("1", "2", "3"))
#Extract the epsi chains and turn to long format. 
long_fit_df <- fit_df[, grep("epsi", colnames(fit_df))] %>% pivot_longer(
             cols = everything(), 
             names_to = c("chain", "var"),
             names_pattern = "chain\\.(\\d{1})\\.epsi_(\\d{1})",
             values_to = "epsi")

# Plot
ggplot(data = full_epsi_df, aes(x = epsi, group = interaction(var,chain), fill = var)) + 
  geom_density(alpha = 0.6) +
  geom_vline(data = epsi_df, aes(xintercept = t_epsi, colour = var), size = 1) + 
  theme_FW(18)

# So we can see that the model is fitting extra variance that is almost exactly 
# the same as the theoretical expectations. This is promising for the models
# that include the extra variance explicitly. 

################################################################################

# Model 2. - Now continue to use just the proportions of resistance (opposed to
# the assay readings), but now also now explicitly fit the variance. 

# Model data:
model_dat <- list(b = b, d = d, nT = nT, del_ts = del_ts, K = K, Ni = Ni,
                  pR_del_t1 = pR_del_t1, pR_del_t2 = pR_del_t2, 
                  pR_del_t3 = pR_del_t3)

# Model fit:
model_fit <- stan(file = "CBC_Rep_Exc_Res_Model2.stan",
                  data = model_dat, chains = 4, iter = 2000)

# Output fits:
mcmc_dens_overlay(model_fit, pars = c("mu", "sigma"))
mcmc_pairs(model_fit, pars = c("mu", "sigma"))
mcmc_trace(model_fit, pars = c("mu", "sigma"))

# Looking good. 

################################################################################

# Model 3. -  For completeness sake, we'll now try the model we know shouldn't 
# work, where we fit the transformed, assay reading means, but don't also 
# include the known, minimum variance. 
# We'll also just include the 'true' x1 and x2 values for now, opposed to 
# estimates. 

# Model data:
model_dat <- list(b = b, d = d, nT = nT, del_ts = del_ts, K = K, Ni = Ni,
                  x_del_t1 = x_del_t1, x_del_t2 = x_del_t2, 
                  x_del_t3 = x_del_t3,
                  x1_est = x1, x2_est = x2)

# Model fit:
model_fit <- stan(file = "CBC_Rep_Exc_Res_Model3.stan",
                  data = model_dat, chains = 4, iter = 2000)

# Output fits:
mcmc_dens_overlay(model_fit, pars = c("mu", "sigma"))
mcmc_pairs(model_fit, pars = c("mu", "sigma"))
mcmc_trace(model_fit, pars = c("mu", "sigma"))

# As expected, poor fits and strong correlations in the posteriors. 

################################################################################

# Model 4. - Now we can include the theoretical expectations for the variance, 
# both for the proportion of resistance, and subsequently for the assay readings. 
# This should get around the correlations in the posteriors you get in Model 3. 
# due to the unidentifiability of the model. 
# Again, for now, we'll just use the 'true' x1 and x2s for the assay
# normalisation. 

# Model data:
model_dat <- list(b = b, d = d, nT = nT, del_ts = del_ts, K = K, Ni = Ni,
                  x_del_t1 = x_del_t1, x_del_t2 = x_del_t2, 
                  x_del_t3 = x_del_t3,
                  x1_est = x1, x2_est = x2)

# Model fit:

model_fit <- stan(file = "CBC_Rep_Exc_Res_Model4.stan",
                  data = model_dat, chains = 4, iter = 5000)

# Output fits:
mcmc_dens_overlay(model_fit, pars = c("m_neg_exp", "s_neg_exp"))
mcmc_pairs(model_fit, pars = c("m_neg_exp", "s_neg_exp"))
mcmc_trace(model_fit, pars = c("m_neg_exp", "s_neg_exp"))

mcmc_dens_overlay(model_fit, pars = c("mu", "sigma"))
mcmc_pairs(model_fit, pars = c("mu", "sigma"))
mcmc_trace(model_fit, pars = c("mu", "sigma"))

# So it does well with a prior for mu weighted towards 0 (which does 
# make sense), otherwise it struggles...

################################
# Lets try with a lower sigma...
################################+

res_df_sub <- subset(res_df, mu == 1e-05 & sig == 1e-03 & Nsim == 1)
res_df_sub <- res_df_sub[, !(names(res_df_sub) %in% "ass_read")]

# Simulate the assay readings 
res_df_sub["ass_read"] <- ass_read(
  rep(x1, nrow(res_df_sub)), 
  rep(x2, nrow(res_df_sub)), 
  rep(mes_err, nrow(res_df_sub)),
  res_df_sub$b, res_df_sub$d, anac(res_df_sub$mu), 
  anac(res_df_sub$sig), anac(res_df_sub$del), 
  res_df_sub$res_prop
)

# Extract the simulated proportions of resistance:
pR_del_t1 <- subset(res_df_sub, delt == del_ts[[1]])$res_prop
pR_del_t2 <- subset(res_df_sub, delt == del_ts[[2]])$res_prop
pR_del_t3 <- subset(res_df_sub, delt == del_ts[[3]])$res_prop

# Extract simulated assay readings 
x_del_t1 <- subset(res_df_sub, delt == del_ts[[1]])$ass_read
x_del_t2 <- subset(res_df_sub, delt == del_ts[[2]])$ass_read
x_del_t3 <- subset(res_df_sub, delt == del_ts[[3]])$ass_read

# Model data:
model_dat <- list(b = b, d = d, nT = nT, del_ts = del_ts, K = K, Ni = Ni,
                  x_del_t1 = x_del_t1, x_del_t2 = x_del_t2, 
                  x_del_t3 = x_del_t3,
                  x1_est = x1, x2_est = x2)

# Model fit:
init_fun <- function(...) list(mu=1e-08, sigma=1e-08)


model_fit <- stan(file = "CBC_Rep_Exc_Res_Model4.stan",
                  data = model_dat, chains = 4, iter = 2000,
                  init = init_fun)

# Output fits:
mcmc_dens_overlay(model_fit, pars = c("mu", "sigma"))
mcmc_pairs(model_fit, pars = c("mu", "sigma"))
mcmc_trace(model_fit, pars = c("mu", "sigma"))

# So the model does seem to struggle, and we see correlations in the posteriors.

fit_df <- data.frame(as.array(model_fit))

# Seems to also think that mu=~0.50 and sigma=5~5.83e-17 is also a valid 
# solution... have a look manually: 

exp_pRs <- exp_res_prop_t(b, d, 18.0, c(1e-05, 0.50), c(1e-03, 5.83e-17))
var_pRs <- var_res_prop_t(b, d, 18.0, c(1e-05, 0.50), c(1e-03, 5.83e-17), 10000)

exp_pRs; var_pRs

exp_ass <- exp_x_read_t(x1, x2, c(1e-05, 0.50), c(1e-03, 5.83e-17), exp_pRs)
var_ass <- var_x_read_t(x1, x2, c(1e-05, 0.50), c(1e-03, 5.83e-17), var_pRs)

exp_ass; var_ass

# So these don't look similar at all, so not sure how they could be the same...

# Lets look at the posteriors directly... 

mu_fit_df <- fit_df[, grep("mu", colnames(fit_df))] %>% pivot_longer(
  cols = everything(), 
  names_to = c("chain"),
  names_pattern = "chain\\.(\\d{1})\\.mu",
  values_to = "mu")

sig_fit_df <- fit_df[, grep("sig", colnames(fit_df))] %>% pivot_longer(
  cols = everything(), 
  names_to = c("chain"),
  names_pattern = "chain\\.(\\d{1})\\.sig",
  values_to = "sig")

full_fit_df <- mu_fit_df
full_fit_df["sig"] <- sig_fit_df$sig

ggplot(data = full_fit_df, aes(x = mu, y = sig)) + 
  geom_point() + 
  facet_grid(~chain, scales = "free")

ggplot(data = subset(full_fit_df, chain == "2"), aes(x = mu, y = sig)) + 
  geom_point() 

ggplot(data = subset(full_fit_df, chain == "1"), aes(x = mu, y = sig)) + 
  geom_point() + 
  geom_point(data = data.frame(x = 1e-05, y = 1e-03), 
             aes(x = x, y = y), colour = "red")

ggplot(data = subset(full_fit_df, chain == "4"), aes(x = mu, y = sig)) + 
  geom_point() + 
  geom_point(data = data.frame(x = 1e-05, y = 1e-03), 
             aes(x = x, y = y), colour = "red")


# Maybe if i plot the likelihood landscape over a huge range of parameters, I 
# could work out if there are 'valleys' of high likelihood that the sampler
# gets stuck in... 

mus <- seq(0.001, 1.0, 0.001)
sigs <- seq(0.001, 1.0, 0.001)

mu_sig_df <- data.frame(mu = rep(mus, length(sigs)), 
                        sig = rep(sigs, each = length(mus)))

exp_pRs <- exp_res_prop_t(b, d, 18.0, mu_sig_df$mu, mu_sig_df$sig)
var_pRs <- var_res_prop_t(b, d, 18.0, mu_sig_df$mu, mu_sig_df$sig, 10000)

mu_sig_df["exp_pR"] <- exp_pRs
mu_sig_df["var_pR"] <- var_pRs

exp_ass <- exp_x_read_t(x1, x2, mu_sig_df$mu, mu_sig_df$sig, exp_pRs)
var_ass <- var_x_read_t(x1, x2, mu_sig_df$mu, mu_sig_df$sig, var_pRs)

mu_sig_df["exp_ass"] <- exp_ass
mu_sig_df["var_ass"] <- var_ass


p1 <- ggplot(data = mu_sig_df, aes(x = mu, y = sig, fill = exp_ass)) + 
  geom_tile() + 
  scale_fill_viridis_c(trans = "log10")

p2 <- ggplot(data = mu_sig_df, aes(x = mu, y = sig, fill = var_ass)) + 
  geom_tile() + 
  scale_fill_viridis_c(trans = "log10")

cowplot::plot_grid(p1, p2, nrow = 1)


ggplot(data = subset(mu_sig_df, sig <= 0.1),
       aes(x = mu, y = sig, fill = var_ass)) + 
  geom_tile() + 
  scale_fill_viridis_c(trans = "log10") + 
  geom_point(data = data.frame(x = 1e-05, y = 1e-03)) + 
  geom_point(data = data.frame(x = 0.50, y = 5.83e-17))


# Need to find a transform that makes these landscapes 'less homogeneous'... 


mus_nexp <- seq(0, 8, 0.05)
sigs_nexp <- seq(0, 8, 0.05)

mu_sig_nexp_df <- data.frame(mu_nexp = rep(mus_nexp, length(sigs_nexp)),
                             sig_nexp = rep(sigs_nexp, each = length(mus_nexp)))

exp_pRs <- exp_res_prop_t(b, d, 18.0, 10^-(mu_sig_nexp_df$mu_nexp), 
                                      10^-(mu_sig_nexp_df$sig_nexp))
var_pRs <- var_res_prop_t(b, d, 18.0, 10^-(mu_sig_nexp_df$mu_nexp),
                                      10^-(mu_sig_nexp_df$sig_nexp), 10000)

mu_sig_nexp_df["exp_pR"] <- exp_pRs
mu_sig_nexp_df["var_pR"] <- var_pRs

exp_ass <- exp_x_read_t(x1, x2, 10^-(mu_sig_nexp_df$mu_nexp), 
                                10^-(mu_sig_nexp_df$sig_nexp), exp_pRs)
var_ass <- var_x_read_t(x1, x2, 10^-(mu_sig_nexp_df$mu_nexp), 
                                10^-(mu_sig_nexp_df$sig_nexp), var_pRs)

mu_sig_nexp_df["exp_ass"] <- exp_ass
mu_sig_nexp_df["var_ass"] <- var_ass

p1 <- ggplot(data = mu_sig_nexp_df, aes(x = mu_nexp, y = sig_nexp, fill = exp_ass)) + 
  geom_tile() +
  scale_fill_viridis_c(trans = "log10")

p2 <- ggplot(data = mu_sig_nexp_df, aes(x = mu_nexp, y = sig_nexp, fill = var_ass)) + 
  geom_tile() +
  scale_fill_viridis_c(trans = "log10")

cowplot::plot_grid(p1, p2, nrow = 1)


#######
# Trying the sum and the ratio of mu and sig to reparam: 
#######


# lam = mu + sig
# phi = mu / sig

# and therefore:

# sig = lam / (phi + 1)
#  mu = phi * lam / (phi + 1)

# Limits: 
# min, max lam = [0, 2]
# min, max phi = [0, Inf]
# AND: 
# if(lam > 1.0), (phi < 1.0)

abs_scal <- 10^-4

lam_raws <- seq(0, (2/abs_scal), 100)

phi_raws <- seq(0, (1000/abs_scal), )

lam_phi_df <- data.frame(lam = rep(lams, length(phis)), 
                         phi = rep(phis, each = length(lams)))

lam_phi_df <- lam_phi_df[which(lam_phi_df$lam < 1.0 | lam_phi_df$phi < 1.0) ,]

lam_phi_df["mu"] <- (lam_phi_df$phi * lam_phi_df$lam) / (lam_phi_df$phi + 1)
lam_phi_df["sig"] <- lam_phi_df$lam / (lam_phi_df$phi + 1)

#############################
# Manually Calculating the LL
#############################

# Decided its best to actually calculate the LL for some given simulated values
# and then plot the LL density - should then be easy to spot any pathologies in 
# the distribution/try out some different re-parameterisations.
# Also has the added bonus of combining the mean and variance AND multiple 
# time-points. 

# The likelihood = the probability of the data | the model.

# Want to input: b, d, nT, del_ts, K, Ni, mu, sigma, and then the simulated 
#                data - x_del_t1, x_del_t2, x_del_t3, and return the LL.
# Can then compare over a range of mus and sigmas. 

ass_read_LL <- function(b, d, del_ts, Ni, mu, sigma, x1, x2,
                        x_del_t1, x_del_t2, x_del_t3){
  
  # Calculate exp(pR) and var(pR) for each del_t.
  exp_pRs <- exp_res_prop_t(b, d, del_ts, mu, sigma)
  var_pRs <- var_res_prop_t(b, d, del_ts, mu, sigma, Ni)
  # And now repeat for the assay readings. 
  exp_ass <- exp_x_read_t(x1, x2, mu, sigma, exp_pRs)
  var_ass <- var_x_read_t(x1, x2, mu, sigma, var_pRs)
  # Return the probabilities for each simulated assay reading...
  lprob_dt1 <- dnorm(x_del_t1, exp_ass[1], sqrt(var_ass[1]), log = T)
  lprob_dt2 <- dnorm(x_del_t2, exp_ass[2], sqrt(var_ass[2]), log = T)
  lprob_dt3 <- dnorm(x_del_t3, exp_ass[3], sqrt(var_ass[3]), log = T)
  # Return the sum of these values. 
  sum_lprob <- sum(lprob_dt1, lprob_dt2, lprob_dt3)
  return(sum_lprob)    
  
}

# Lets do some tests to see if this is working: 


# First, add the assay readings to all of the res_df. 
res_df["ass_read"] <- ass_read(
  rep(x1, nrow(res_df)), 
  rep(x2, nrow(res_df)), 
  rep(mes_err, nrow(res_df_sub)),
  res_df$b, res_df$d, anac(res_df$mu), 
  anac(res_df$sig), anac(res_df$del), 
  res_df$res_prop
)

# vectors of mus and sigmas to test over: 

mus <- 10^-(seq(0, 10, 0.05))
sigs <- 10^-(seq(0, 10, 0.05))

mu_sig_df <- data.frame(mu = rep(mus, length(sigs)), 
                        sig = rep(sigs, each = length(mus)))

# Now choose some subsets and plot the change in LL over a param...

df_sub1 <- subset(res_df, Nsim == 1 & mu == 1e-04 & sig == 1e-01)

del_ts <- unique(df_sub1$delt)
x_del_t1_1 <- subset(df_sub1, delt == del_ts[1])$ass_read
x_del_t2_1 <- subset(df_sub1, delt == del_ts[2])$ass_read
x_del_t3_1 <- subset(df_sub1, delt == del_ts[3])$ass_read

mu_sig_df["LL"] <- NA 

# ! THIS IS REALLY SLOW - CHANGE TO AN APPLY VERSION...

for(i in 1:nrow(mu_sig_df)){
  
  mu_sig_df[i, "LL"] <- ass_read_LL(b, d, del_ts, Ni, 
                                  mu_sig_df[i ,]$mu, mu_sig_df[i ,]$sig,
                                  x1, x2, 
                                  x_del_t1_1, x_del_t2_1, x_del_t3_1)
  
}

# Disttribution of LLs:

ggplot(data = mu_sig_df, aes(x = -LL)) + 
  geom_histogram(bins = 100) + 
  scale_x_log10() + 
  scale_y_log10()

# to aid plotting...

# Plot the likelihoods for these combinations:

p1 <- ggplot(data = mu_sig_df, aes(x = -log10(mu), y = -log10(sig),
                             fill = -LL+(max(LL)+1))) + 
  geom_tile() + 
  scale_fill_viridis_c(trans = "log10", option = "magma") + 
  geom_point(data = data.frame(x = 4, y = 1,LL=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)


p2 <- ggplot(data = subset(mu_sig_df, -log10(mu) <= 10 & -log10(sig) <= 2), 
       aes(x = -log10(mu), y = -log10(sig),
           fill = -LL+(max(LL)+1))) + 
  geom_tile() + 
  scale_fill_viridis_c(trans = "log10", option = "magma") + 
  geom_point(data = data.frame(x = 4, y = 1,LL=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)

cowplot::plot_grid(p1, p2, nrow = 1)

# And with a smaller mu and smaller sig...

df_sub1 <- subset(res_df, Nsim == 1 & mu == 1e-05 & sig == 1e-04)

del_ts <- unique(df_sub1$delt)
x_del_t1_1 <- subset(df_sub1, delt == del_ts[1])$ass_read
x_del_t2_1 <- subset(df_sub1, delt == del_ts[2])$ass_read
x_del_t3_1 <- subset(df_sub1, delt == del_ts[3])$ass_read

mu_sig_df["LL"] <- NA 

# ! THIS IS REALLY SLOW - CHANGE TO AN APPLY VERSION...

for(i in 1:nrow(mu_sig_df)){
  
  mu_sig_df[i, "LL"] <- ass_read_LL(b, d, del_ts, Ni, 
                                    mu_sig_df[i ,]$mu, mu_sig_df[i ,]$sig,
                                    x1, x2, 
                                    x_del_t1_1, x_del_t2_1, x_del_t3_1)
  
}

# to aid plotting...

# Plot the likelihoods for these combinations:

p1 <- ggplot(data = mu_sig_df, aes(x = -log10(mu), y = -log10(sig),
                                   fill = -LL+(max(LL)+1))) + 
  geom_tile() + 
  scale_fill_viridis_c(trans = "log10", option = "magma") + 
  geom_point(data = data.frame(x = 5, y = 4,LL=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)


p2 <- ggplot(data = subset(mu_sig_df, -log10(mu) > 3 & -log10(sig) <= 5 & 
                             -log10(sig) >= 3), 
             aes(x = -log10(mu), y = -log10(sig),
                 fill = -LL+(max(LL)+1))) + 
  geom_tile() + 
  scale_fill_viridis_c(trans = "log10", option = "magma") + 
  geom_point(data = data.frame(x = 5, y = 4,LL=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)

cowplot::plot_grid(p1, p2, nrow = 1)


# We can now also look at how a re-param changes the LL 'landscape'...

# REPARAM 1.:

# lam = mu + sig         Therefore min, max = [0, 2.0]
# phi = mu / (mu + sig)  Therefore min, max = [0, 1.0]

# sig = lam - (lam * phi)
#  mu = lam * phi

lams <- 2*10^-seq(0, 10, 0.05)
phis <- 10^-(seq(0, 10, 0.05))

#lams <- seq(0, 2.0, 0.02)
#phis <- seq(0, 1.0, 0.01)

lam_phi_df <- data.frame(lam = rep(lams, length(phis)), 
                         phi = rep(phis, each = length(lams)))

lam_phi_df["mu"] <- lam_phi_df$lam * lam_phi_df$phi
lam_phi_df["sig"] <- lam_phi_df$lam - (lam_phi_df$lam * lam_phi_df$phi)

lam_phi_df <- subset(lam_phi_df, mu <= 1.00)
lam_phi_df <- subset(lam_phi_df, sig <= 1.00)

df_sub1 <- subset(res_df, Nsim == 1 & mu == 1e-05 & sig == 1e-04)

del_ts <- unique(df_sub1$delt)
x_del_t1_1 <- subset(df_sub1, delt == del_ts[1])$ass_read
x_del_t2_1 <- subset(df_sub1, delt == del_ts[2])$ass_read
x_del_t3_1 <- subset(df_sub1, delt == del_ts[3])$ass_read

lam_phi_df["LL"] <- NA 

# ! THIS IS REALLY SLOW - CHANGE TO AN APPLY VERSION...

for(i in 1:nrow(lam_phi_df)){
  
  lam_phi_df[i, "LL"] <- ass_read_LL(b, d, del_ts, Ni, 
                                    lam_phi_df[i ,]$mu, lam_phi_df[i ,]$sig,
                                    x1, x2, 
                                    x_del_t1_1, x_del_t2_1, x_del_t3_1)
  
}


# True values: 
true_lam <- 1e-05 + 1e-04
true_phi <- 1e-05 / (1e-05 + 1e-04)

p1 <- ggplot(data = lam_phi_df, aes(x = -log10(lam), y = -log10(phi),
                                   fill = -LL+(max(LL,na.rm = T)+1))) + 
  geom_tile() + 
  scale_fill_viridis_c(trans = "log10", option = "magma", name = "-LL(...)") + 
  geom_point(data = data.frame(x = -log10(true_lam), y = -log10(true_phi),
                               LL=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)


pa <- ggplot(data = lam_phi_df, aes(x = -log10(lam), y = -log10(phi),
                              fill = -log10(mu))) + 
  geom_tile() + 
  scale_fill_viridis_c(option = "magma") + 
  geom_point(data = data.frame(x = -log10(true_lam), y = -log10(true_phi),
                               LL=NA,mu=NA,sig=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)

pb <- ggplot(data = lam_phi_df, aes(x = -log10(lam), y = -log10(phi),
                              fill = -log10(sig))) + 
  geom_tile() + 
  scale_fill_viridis_c(option = "magma") + 
  geom_point(data = data.frame(x = -log10(true_lam), y = -log10(true_phi),
                               LL=NA,mu=NA,sig=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)

cowplot::plot_grid(p1, p1, pa, pb, nrow = 2, ncol = 2)

# SUBSET...

p1 <- ggplot(data = subset(lam_phi_df, -log10(lam) > 3.8 & -log10(lam) < 4.2), 
             aes(x = -log10(lam), y = -log10(phi),
                                    fill = -LL+(max(LL,na.rm = T)+1))) + 
  geom_tile() + 
  scale_fill_viridis_c(trans = "log10", option = "magma", name = "-LL(...)") + 
  geom_point(data = data.frame(x = -log10(true_lam), y = -log10(true_phi),
                               LL=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)


pa <- ggplot(data = subset(lam_phi_df, -log10(lam) > 3.8 & -log10(lam) < 4.2),
             aes(x = -log10(lam), y = -log10(phi),
                                    fill = -log10(mu))) + 
  geom_tile() + 
  scale_fill_viridis_c(option = "magma") + 
  geom_point(data = data.frame(x = -log10(true_lam), y = -log10(true_phi),
                               LL=NA,mu=NA,sig=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)

pb <- ggplot(data = subset(lam_phi_df, -log10(lam) > 3.8 & -log10(lam) < 4.2),
             aes(x = -log10(lam), y = -log10(phi),
                                    fill = -log10(sig))) + 
  geom_tile() + 
  scale_fill_viridis_c(option = "magma") + 
  geom_point(data = data.frame(x = -log10(true_lam), y = -log10(true_phi),
                               LL=NA,mu=NA,sig=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)

cowplot::plot_grid(p1, p1, pa, pb, nrow = 2, ncol = 2)

######################
# DIFFERENT PARAM SET:
######################

df_sub1 <- subset(res_df, Nsim == 1 & mu == 1e-05 & sig == 1e-01)

del_ts <- unique(df_sub1$delt)
x_del_t1_1 <- subset(df_sub1, delt == del_ts[1])$ass_read
x_del_t2_1 <- subset(df_sub1, delt == del_ts[2])$ass_read
x_del_t3_1 <- subset(df_sub1, delt == del_ts[3])$ass_read

lam_phi_df["LL"] <- NA 

# ! THIS IS REALLY SLOW - CHANGE TO AN APPLY VERSION...

for(i in 1:nrow(lam_phi_df)){
  
  lam_phi_df[i, "LL"] <- ass_read_LL(b, d, del_ts, Ni, 
                                     lam_phi_df[i ,]$mu, lam_phi_df[i ,]$sig,
                                     x1, x2, 
                                     x_del_t1_1, x_del_t2_1, x_del_t3_1)
  
}


# True values: 
true_lam <- 1e-05 + 1e-01
true_phi <- 1e-05 / (1e-05 + 1e-01)

p1 <- ggplot(data = lam_phi_df, aes(x = -log10(lam), y = -log10(phi),
                                    fill = -LL+(max(LL,na.rm = T)+1))) + 
  geom_tile() + 
  scale_fill_viridis_c(trans = "log10", option = "magma", name = "-LL(...)") + 
  geom_point(data = data.frame(x = -log10(true_lam), y = -log10(true_phi),
                               LL=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)


pa <- ggplot(data = lam_phi_df, aes(x = -log10(lam), y = -log10(phi),
                                    fill = -log10(mu))) + 
  geom_tile() + 
  scale_fill_viridis_c(option = "magma") + 
  geom_point(data = data.frame(x = -log10(true_lam), y = -log10(true_phi),
                               LL=NA,mu=NA,sig=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)

pb <- ggplot(data = lam_phi_df, aes(x = -log10(lam), y = -log10(phi),
                                    fill = -log10(sig))) + 
  geom_tile() + 
  scale_fill_viridis_c(option = "magma") + 
  geom_point(data = data.frame(x = -log10(true_lam), y = -log10(true_phi),
                               LL=NA,mu=NA,sig=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)

cowplot::plot_grid(p1, p1, pa, pb, nrow = 2, ncol = 2)

# SUBSET...

p1 <- ggplot(data = subset(lam_phi_df, -log10(lam) > 0.8 & -log10(lam) < 1.2), 
             aes(x = -log10(lam), y = -log10(phi),
                 fill = -LL+(max(LL,na.rm = T)+1))) + 
  geom_tile() + 
  scale_fill_viridis_c(trans = "log10", option = "magma", name = "-LL(...)") + 
  geom_point(data = data.frame(x = -log10(true_lam), y = -log10(true_phi),
                               LL=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)


pa <- ggplot(data = subset(lam_phi_df, -log10(lam) > 0.8 & -log10(lam) < 1.2),
             aes(x = -log10(lam), y = -log10(phi),
                 fill = -log10(mu))) + 
  geom_tile() + 
  scale_fill_viridis_c(option = "magma") + 
  geom_point(data = data.frame(x = -log10(true_lam), y = -log10(true_phi),
                               LL=NA,mu=NA,sig=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)

pb <- ggplot(data = subset(lam_phi_df, -log10(lam) > 0.8 & -log10(lam) < 1.2),
             aes(x = -log10(lam), y = -log10(phi),
                 fill = -log10(sig))) + 
  geom_tile() + 
  scale_fill_viridis_c(option = "magma") + 
  geom_point(data = data.frame(x = -log10(true_lam), y = -log10(true_phi),
                               LL=NA,mu=NA,sig=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)

cowplot::plot_grid(p1, p1, pa, pb, nrow = 2, ncol = 2)


########### 
# Seeing if the values of phi it struggles with all lead to similar Reqs. 
          
head(lam_phi_df)

lam_phi_df["R_eq"] <- res_eq_freq(lam_phi_df$mu, lam_phi_df$sig)

p1 <- ggplot(data = lam_phi_df, aes(x = -log10(lam), y = -log10(phi),
                                    fill = -LL+(max(LL,na.rm = T)+1))) + 
  geom_tile() + 
  scale_fill_viridis_c(trans = "log10", option = "magma", name = "-LL(...)") + 
  geom_point(data = data.frame(x = -log10(true_lam), y = -log10(true_phi),
                               LL=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)

p2 <- ggplot(data = lam_phi_df, aes(x = -log10(lam), y = -log10(phi),
                                    fill = R_eq)) + 
  geom_tile() + 
  scale_fill_viridis_c(trans = "log10", option = "magma") + 
  geom_point(data = data.frame(x = -log10(true_lam), y = -log10(true_phi),
                               LL=NA,R_eq=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)

cowplot::plot_grid(p1, p2, nrow = 1)

# NO - NOT THIS. 

# REPARAM 2.:

# lam = mu + sig          Therefore min, max = [0, 2.0]
# phi = sig / (mu + sig)  Therefore min, max = [0, 1.0]

# mu = lam - (lam * phi)
# sig = lam * phi

lams <- 2*10^-seq(0, 10, 0.05)
phis <- 10^-(seq(0, 10, 0.05))

lam_phi_df <- data.frame(lam = rep(lams, length(phis)), 
                         phi = rep(phis, each = length(lams)))

lam_phi_df["sig"] <- lam_phi_df$lam * lam_phi_df$phi
lam_phi_df["mu"] <- lam_phi_df$lam - (lam_phi_df$lam * lam_phi_df$phi)

lam_phi_df <- subset(lam_phi_df, mu <= 1.00)
lam_phi_df <- subset(lam_phi_df, sig <= 1.00)

df_sub1 <- subset(res_df, Nsim == 1 & mu == 1e-05 & sig == 1e-04)

del_ts <- unique(df_sub1$delt)
x_del_t1_1 <- subset(df_sub1, delt == del_ts[1])$ass_read
x_del_t2_1 <- subset(df_sub1, delt == del_ts[2])$ass_read
x_del_t3_1 <- subset(df_sub1, delt == del_ts[3])$ass_read

lam_phi_df["LL"] <- NA 

# ! THIS IS REALLY SLOW - CHANGE TO AN APPLY VERSION...

for(i in 1:nrow(lam_phi_df)){
  
  lam_phi_df[i, "LL"] <- ass_read_LL(b, d, del_ts, Ni, 
                                     lam_phi_df[i ,]$mu, lam_phi_df[i ,]$sig,
                                     x1, x2, 
                                     x_del_t1_1, x_del_t2_1, x_del_t3_1)
  
}


# True values: 
true_lam <- 1e-05 + 1e-04
true_phi <- 1e-04 / (1e-05 + 1e-04)

p1 <- ggplot(data = lam_phi_df, aes(x = -log10(lam), y = -log10(phi),
                                    fill = -LL+(max(LL,na.rm = T)+1))) + 
  geom_tile() + 
  scale_fill_viridis_c(trans = "log10", option = "magma", name = "-LL(...)") + 
  geom_point(data = data.frame(x = -log10(true_lam), y = -log10(true_phi),
                               LL=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)

p2 <- ggplot(data = subset(lam_phi_df, -log10(lam) > 2.5 & -log10(lam) < 5.0 &
                             -log10(phi) < 2.0), 
             aes(x = -log10(lam), y = -log10(phi),
                                    fill = -LL+(max(LL,na.rm = T)+1))) + 
  geom_tile() + 
  scale_fill_viridis_c(trans = "log10", option = "magma", name = "-LL(...)") + 
  geom_point(data = data.frame(x = -log10(true_lam), y = -log10(true_phi),
                               LL=NA),aes(x=x,y=y),
             shape = 8, colour = "white", size = 3) + 
  theme_FW(12)


cowplot::plot_grid(p1, p2, nrow = 1)




################################################################################

# Model 5. - The model is struggling with the assay-transformed version (even
# when also explicitly including the minimum variance expected).  

res_df_sub <- subset(res_df, mu == 1e-05 & sig == 1e-01 & Nsim == 1)
res_df_sub <- res_df_sub[, !(names(res_df_sub) %in% "ass_read")]

# Simulate the assay readings 
res_df_sub["ass_read"] <- ass_read(
  rep(x1, nrow(res_df_sub)), 
  rep(x2, nrow(res_df_sub)), 
  rep(mes_err, nrow(res_df_sub)),
  res_df_sub$b, res_df_sub$d, anac(res_df_sub$mu), 
  anac(res_df_sub$sig), anac(res_df_sub$del), 
  res_df_sub$res_prop
)

# Extract the simulated proportions of resistance:
pR_del_t1 <- subset(res_df_sub, delt == del_ts[[1]])$res_prop
pR_del_t2 <- subset(res_df_sub, delt == del_ts[[2]])$res_prop
pR_del_t3 <- subset(res_df_sub, delt == del_ts[[3]])$res_prop

# Extract simulated assay readings 
x_del_t1 <- subset(res_df_sub, delt == del_ts[[1]])$ass_read
x_del_t2 <- subset(res_df_sub, delt == del_ts[[2]])$ass_read
x_del_t3 <- subset(res_df_sub, delt == del_ts[[3]])$ass_read


# Model data:
model_dat <- list(b = b, d = d, nT = nT, del_ts = del_ts, K = K, Ni = Ni,
                  x_del_t1 = x_del_t1, x_del_t2 = x_del_t2, 
                  x_del_t3 = x_del_t3,
                  x1_est = x1, x2_est = x2)

# Model fit:
model_fit <- stan(file = "CBC_Rep_Exc_Res_Model5.stan",
                  data = model_dat, chains = 4, iter = 5000)


# true lam and phi...

#lam:
true_lam <- anac(unique(res_df_sub$mu)) + anac(unique(res_df_sub$sig))
#phi:
true_phi <- anac(unique(res_df_sub$mu)) / (anac(unique(res_df_sub$mu)) + anac(unique(res_df_sub$sig)))

# Output fits:
mcmc_dens_overlay(model_fit, pars = c("lam_neg_exp", "phi_neg_exp"))
mcmc_pairs(model_fit, pars = c("lam_neg_exp", "phi_neg_exp"))
mcmc_trace(model_fit, pars = c("lam_neg_exp", "phi_neg_exp"))

mcmc_dens_overlay(model_fit, pars = c("lam", "phi"))
mcmc_pairs(model_fit, pars = c("lam", "phi"))
mcmc_trace(model_fit, pars = c("lam", "phi"))

mcmc_dens_overlay(model_fit, pars = c("mu", "sigma"))
mcmc_pairs(model_fit, pars = c("mu", "sigma"))
mcmc_trace(model_fit, pars = c("mu", "sigma"))

# Compare the posterior with the prior for phi to see if there is any 
# shrinkage or its just fitting to the prior....


fit_df <- data.frame(as.array(model_fit))

long_fit_df <- fit_df[, grep("phi_neg_exp", colnames(fit_df))] %>% pivot_longer(
  cols = everything(), 
  names_to = c("chain"),
  names_pattern = "chain\\.(\\d{1})\\.phi_neg_exp",
  values_to = "phi_neg_exp")

ggplot(data = long_fit_df, aes(x = phi_neg_exp, fill = "posterior")) + 
  geom_density(alpha = 0.6) + 
  geom_area(data = data.frame(x = seq(0, 20, 0.01), 
                              y = dgamma(seq(0, 20, 0.01), 5.0, 1.0)),
            aes(x=x,y=y,fill="prior"),alpha=0.6) + 
  geom_vline(xintercept = -log10(true_phi), size = 1, linetype = "dashed") + 
  geom_vline(xintercept = mean(long_fit_df$phi_neg_exp), size = 1, 
             colour = "red", linetype = "dashed")
  
########
# AND WITH A LOWER SIG..
########

res_df_sub <- subset(res_df, mu == 1e-05 & sig == 1e-03 & Nsim == 1)
res_df_sub <- res_df_sub[, !(names(res_df_sub) %in% "ass_read")]

# Simulate the assay readings 
res_df_sub["ass_read"] <- ass_read(
  rep(x1, nrow(res_df_sub)), 
  rep(x2, nrow(res_df_sub)), 
  rep(mes_err, nrow(res_df_sub)),
  res_df_sub$b, res_df_sub$d, anac(res_df_sub$mu), 
  anac(res_df_sub$sig), anac(res_df_sub$del), 
  res_df_sub$res_prop
)

# Extract the simulated proportions of resistance:
pR_del_t1 <- subset(res_df_sub, delt == del_ts[[1]])$res_prop
pR_del_t2 <- subset(res_df_sub, delt == del_ts[[2]])$res_prop
pR_del_t3 <- subset(res_df_sub, delt == del_ts[[3]])$res_prop

# Extract simulated assay readings 
x_del_t1 <- subset(res_df_sub, delt == del_ts[[1]])$ass_read
x_del_t2 <- subset(res_df_sub, delt == del_ts[[2]])$ass_read
x_del_t3 <- subset(res_df_sub, delt == del_ts[[3]])$ass_read

# Model data:
model_dat <- list(b = b, d = d, nT = nT, del_ts = del_ts, K = K, Ni = Ni,
                  x_del_t1 = x_del_t1, x_del_t2 = x_del_t2, 
                  x_del_t3 = x_del_t3,
                  x1_est = x1, x2_est = x2)

# Model fit:
model_fit <- stan(file = "CBC_Rep_Exc_Res_Model5.stan",
                  data = model_dat, chains = 4, iter = 5000)

# true lam and phi...

#lam:
true_lam <- anac(unique(res_df_sub$mu)) + anac(unique(res_df_sub$sig))
#phi:
true_phi <- anac(unique(res_df_sub$sig)) / (anac(unique(res_df_sub$mu)) + anac(unique(res_df_sub$sig)))

# Output fits:
mcmc_dens_overlay(model_fit, pars = c("lam_neg_exp", "phi_neg_exp"))
mcmc_pairs(model_fit, pars = c("lam_neg_exp", "phi_neg_exp"))
mcmc_trace(model_fit, pars = c("lam_neg_exp", "phi_neg_exp"))

mcmc_dens_overlay(model_fit, pars = c("lam", "phi"))
mcmc_pairs(model_fit, pars = c("lam", "phi"))
mcmc_trace(model_fit, pars = c("lam", "phi"))

mcmc_dens_overlay(model_fit, pars = c("mu", "sigma"))
mcmc_pairs(model_fit, pars = c("mu", "sigma"))
mcmc_trace(model_fit, pars = c("mu", "sigma"))


# Compare the posterior with the prior for phi to see if there is any 
# shrinkage or its just fitting to the prior....

fit_df <- data.frame(as.array(model_fit))

long_fit_df <- fit_df[, grep("phi_neg_exp", colnames(fit_df))] %>% pivot_longer(
  cols = everything(), 
  names_to = c("chain"),
  names_pattern = "chain\\.(\\d{1})\\.phi_neg_exp",
  values_to = "phi_neg_exp")

ggplot(data = long_fit_df, aes(x = phi_neg_exp, fill = "posterior")) + 
  geom_density(alpha = 0.6) + 
  geom_area(data = data.frame(x = seq(0, 0.10, 0.0001), 
                              y = dgamma(seq(0, 0.10, 0.0001), 5.0, 1.0)),
            aes(x=x,y=y,fill="prior"),alpha=0.6) + 
  geom_vline(xintercept = -log10(true_phi), size = 1, linetype = "dashed") + 
  geom_vline(xintercept = mean(long_fit_df$phi_neg_exp), size = 1, 
             colour = "red", linetype = "dashed")


################################################################################


