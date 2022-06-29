
# Cancer Barcode Simulation
# Freddie Whiting - 2022

# Looking through the outputs of the 'CBC-Rep_Exc_Res' (replicate, exclusive
# resistance) simulations. 

################################################################################

rm(list = ls())

#source("C:/Users/fwhiting/My Drive/Barcode_Simulations/Cancer_BarCode_Sim_CBC/CBC_Analysis_Custom_Functions.R")
#source("C:/Users/whitin01/Google Drive/Barcode_Simulations/Cancer_BarCode_Sim_CBC/CBC_Analysis_Custom_Functions.R")
source("C:/Users/Freddie/Google Drive/Barcode_Simulations/Cancer_BarCode_Sim_CBC/CBC_Analysis_Custom_Functions.R")

#setwd("C:/Users/fwhiting/My Drive/Barcode_Simulations/Cancer_BarCode_Sim_CBC/Outputs/Rep_Exc_Res/")
#setwd("C:/Users/whitin01/Google Drive/Barcode_Simulations/Cancer_BarCode_Sim_CBC/Outputs/Rep_Exc_Res/")
setwd("C:/Users/Freddie/Google Drive/Barcode_Simulations/Cancer_BarCode_Sim_CBC/Outputs/Rep_Exc_Res/")

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

# Compare the change in the resistant proportion over time in 

# i) the expanded populations (unobservable experimentally):

# ggplot(data = exp_res_df, aes(x = del_t, y = res_prop, fill = sig)) +
#   geom_point(size = 4, shape = 21, alpha = 0.6, colour = "black") +
#   theme_FW(12) +
#   scale_fill_manual(values = pal) +
#   ggtitle("Expanded Populations") +
#   xlab("\U0394t") +
#   ylab("Resistant Proportion (p)") + 
#   facet_wrap(~mu, ncol = 1)

# ii) the replicate, sampled populations (observable experimentally):

# ggplot(data = res_df, aes(x = delt, y = res_prop, fill = sig)) +
#   geom_point(size = 4, shape = 21, alpha = 0.6, colour = "black") +
#   theme_FW(12) +
#   scale_fill_manual(values = pal) +
#   ggtitle("Replicate Sampled Populations") +
#   xlab("\U0394t") +
#   ylab("Resistant Proportion (p)") +
#   scale_x_continuous(limits = c(0, max(res_df$delt))) +
#   facet_wrap(~mu, ncol = 1)

# Update: is there information held in the variance of the resistant 
# proportion?... 

summ_df_rp <- plyr::ddply(res_df, .(mu, sig, delt), summarise, 
                       mean = mean(res_prop), var = var(res_prop))

summ_df_nr <- plyr::ddply(res_df, .(mu, sig, delt), summarise, 
                          mean = mean(res_n), var = var(res_n))


ggplot(data = summ_df_rp, aes(x = mean, y = var, fill = sig)) + 
  geom_point(shape = 21, colour = "black", size = 4) + 
  facet_grid(~mu)

# Can see that the variance peaks when the mean proportion of resistance 
# = 0.5, and then starts to decrease again, with there being more information 
# the closer to 0.5 you get (which happens when sigma is high... because 
# then you have had time in the time windows to move significantly to either be
# approaching or traveled past pR = 0.5...)

# It looks like as long as mu is 'low' (<= 1e-03) then it doesn't influence
# these values much...



########################################################
# Compare simulated outputs to theoretical expectations.
########################################################

# A function that is the solution to the proportion of resistance
# at time t given the initial condition p = 1.0 (where p = nR / (nr  + nR)) = 

# Solution: 

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

# Make a dataframe with the theoretical expectations over the given: 
#  i) sigma values, 
# ii) time windows,

mus <- anac(unique(res_df$mu))
sigs <- anac(unique(res_df$sig))

mus <- rep(mus, length(sigs))
sigs <- rep(sigs, each = length(anac(unique(res_df$mu))))
ts <- unique(res_df$delt)
# expand the times
ts_exp <- seq(0.0, max(ts), 0.1)
otl <- length(ts_exp)
ts_exp <- rep(ts_exp, length(mus))

mus <- rep(mus, each = otl)
sigs <- rep(sigs, each = otl)

res_props <- res_prop_t(mus, sigs, 0.0, unique(res_df$b), unique(res_df$d), 
                        ts_exp)


# theory dataframe: 
tdf <- data.frame(mu = as.factor(mus), sig = as.factor(sigs), 
                  delt = ts_exp, res_prop = res_props)

# Can now compare these theoretical expectations with the simulated values:

# ggplot(data = res_df, aes(x = delt, y = res_prop, fill = sig)) +
#   geom_point(size = 4, shape = 21, alpha = 0.6, colour = "black") +
#   geom_line(data = tdf, aes(x = delt, y = res_prop, colour = sig), size = 1) +
#   facet_wrap(~mu, ncol = 1) +
#   theme_FW(12) +
#   scale_fill_manual(values = pal) +
#   scale_colour_manual(values =  pal) +
#   ggtitle("Replicate Sampled Populations") +
#   xlab("\U0394t") +
#   ylab("Resistant Proportion (p)") +
#   scale_x_continuous(limits = c(0, max(res_df$delt)))

# So we can see that the simulations accurately replicate the expected resistant
# proportions over time, given some parameters that control the evolution of the
# resistant phenotype.
# Here i've looked at mu = 10^-5. The power to detect differences will diminish
# as:
# 1. mu increases - this will lead to a higher equilibrium frequency of
#    resistance and therefore the relative difference between the initial
#    condition (p = 1.0) and the equilibrium frequency of resistance (eqR).
# 2. sig decreases - for the same reasons as 1.

###############################################################################

# Update: explicitly including the expected variance. 

# Q - can i derive an analytical expectation of what the variance should be...
# i'm guessing it will depend on the birth-death process, and the population/
#bottleneck sizes... 

# So think I have worked out an analytical solution to the variance of the 
# proportion of resistance at a given time, t, given the initial condition
# that pR (the proportion of resistance) = 1.0, i.e. pR(t)=1.0|t=0.0. 

# If we consider the two-state continuous state markov chain where we ignore
# population size and the birth-death model, where we just have a single cell
# that can transition from sensitive to resistant with rate mu, and resistant
# to sensitive with rate sigma - we know that the expected 'proportion' of 
# resistance = res_prop_t(mu ,sig, del, b, d, t) (see function below). 

# Can consider this a p.m.f where the probability of resistance at time t 
# is the output of this function, and sensitivity is the complement...
# and can use this to calculate the variance: 
# Var(X) = E[X^2] - (E[X])^2. 
# So if, for example, E[pR(t)] = 0.7, (the probability the single cell, i.e.
#                                      proportion of resistance, pR = 1.0, at 
#                                      time t). 
#                     Var(X) = 0.7  - (0.7)^2 = 0.21

# Write a function that will return the expected variance given a two-state
# markov chain with switching rates mu and sig (ignoring the birth-death)
# process for now.

exp_rp_var <- function(mu, sig, del, b, d, t){
  
  np <- res_prop_t(mu, sig, del, b, d, t)
  
  ev_p <- np - (np^2)
  
  return(ev_p)
  
}

# This gives the variance given a single cell switching between resistant and 
# sensitive, independent of a birth-death process. 
# The mean is given as a proportion, so this holds when looking at the 
# relative frequencies, however need to rescale the variance;
#   E[rel_freq_R] = E[nR]/N
# var[rel_freq_R] = var[nR]/N^2

summ_df_nr["exp_nR"] <- res_prop_t(anac(summ_df_nr$mu), anac(summ_df_nr$sig),
                                   0.0, unique(res_df$b), unique(res_df$d), 
                                   anac(summ_df_nr$delt))*unique(res_df$Ni)

summ_df_nr["var_nR"] <- exp_rp_var(anac(summ_df_nr$mu), anac(summ_df_nr$sig),
                                   0.0, unique(res_df$b), unique(res_df$d),
                                   anac(summ_df_nr$delt))*unique(res_df$Ni)

ggplot(data = summ_df_nr, aes(x = mean, y = var, fill = sig)) + 
  geom_point(shape = 21, colour = "black", size = 4) + 
  geom_point(aes(x = exp_nR, y = var_nR)) +
  facet_grid(~mu)

summ_df_rp["exp_pR"] <- res_prop_t(anac(summ_df_nr$mu), anac(summ_df_nr$sig),
                                   0.0, unique(res_df$b), unique(res_df$d), 
                                   anac(summ_df_nr$delt))

summ_df_rp["var_pR"] <- exp_rp_var(anac(summ_df_nr$mu), anac(summ_df_nr$sig),
                                   0.0, unique(res_df$b), unique(res_df$d),
                                   anac(summ_df_nr$delt))/unique(res_df$Ni)

ggplot(data = summ_df_rp, aes(x = mean, y = var, fill = sig)) + 
  geom_point(shape = 21, colour = "black", size = 4) + 
  geom_point(aes(x = exp_pR, y = var_pR)) +
  facet_grid(~mu)

# So we can see that looking at these estimates alone, we under-estimate the
# variance - think we'll need to also include the birth-death variance to 
# accurately estimate the variances and subsequently use this to fit the model. 

# Actually! - think that the observed variances are actually primarily due to 
# the binomial sampling from the expanded populations -> replicate 
# sub-populations. 

# Can look at the expanded populations pre-sampling... 

summ_df_nr2 <- plyr::ddply(subset(exp_res_df, del_t != 0), 
                           .(mu, sig, del_t), summarise, 
                          mean = mean(res_n), var = var(res_n))

summ_df_nr2["exp_nR"] <- res_prop_t(anac(summ_df_nr2$mu), anac(summ_df_nr2$sig),
                                   0.0, unique(exp_res_df$b), unique(exp_res_df$d), 
                                   anac(summ_df_nr2$del_t))*6.4*10^6

summ_df_nr2["var_nR"] <- exp_rp_var(anac(summ_df_nr2$mu), anac(summ_df_nr2$sig),
                                   0.0, unique(exp_res_df$b), unique(exp_res_df$d), 
                                   anac(summ_df_nr2$del_t))*6.4*10^6

ggplot(data = summ_df_nr2, aes(x = mean, y = var, fill = sig)) + 
  geom_point(shape = 21, colour = "black", size = 4) + 
  geom_point(aes(x = exp_nR, y = var_nR)) + 
  facet_grid(del_t~mu) + 
  scale_x_log10()

# So this isn't calculating the expected proportion of resistance 

# Can it instead be explained by the binomial sampling variance? 

bin_ps <- summ_df_nr$exp_nR/unique(res_df$Ni)
bin_qs <- (1 - bin_ps)
bin_ns <- unique(res_df$Ni)

summ_df_nr["binom_var"] <- bin_ps * bin_qs * bin_ns

ggplot(data = summ_df_nr, aes(x = mean, y = var, fill = sig)) + 
  geom_point(shape = 21, colour = "black", size = 4) + 
  geom_point(aes(x = exp_nR, y = binom_var)) +
  facet_grid(~mu)

# So the binomial variance is identical to that calculated by the 
# Var[pR(t)] = E[pR(t)^2] - (e[pR(t)])^) ...?
# ... and is therefore also still an underestimate. 

# In the un-sampled populations, the estimated variance is a huge underestimate
# - I therefore think that there is some variance in the birth-death process
# that i'm not considering, and that this is highest pre-sampling, following
# the extended bith-death process for each del_t-i.... 

exp_Nt <- function(N0, b, d, t){
  
  expNt <- N0 * exp((b - d)*t)
  
  return(expNt)
  
}

var_Nt <- function(N0, b, d, t){
  
  varn1 <- (b+d)/(b-d) * exp((b - d)*t) * (exp((b - d)*t) - 1)
  
  varNt <- N0*varn1
  
  return(varNt)
  
}

# UPDATE: I think the solution to the two-state markov chain 'super-imposed' 
# onto the birth-death model (also known as a multi-state branching process)
# is too complicated... instead, I think i can make use of the fact that i know
# the minimum variance should be that of the binomial sampling... and then i fit
# unknown variance on top of this. I think this should prevent the re-normalised
# model (following the simulated 'assay readings') being unidentifiable. 


################################################################################



# The next step is to also simulate the process of measuring the proportion of
# resistance in a replicate sub-population (i.e. a 96-well) as part of an 
# in vitro growth assay. 

# Set a value for the measurement error

mes_err <- 0.0

# Choose the number of measurements to take per condition. 

mes_n <- K

# First step is to get readings for what 'exclusive resistance' would look like. 
# If we assume that there is no cost of resistance, this would be the same as
# sensitive cells in a non-treated environment. 
# We will call these readings 'x2'. 
# Set some arbitrary, true value for 'reading x2'.

x2 <- 10

# Now we need some measurements for the equilibrium frequency of resistance 
# given the evolutionary parameters (eqR). This can be found by exposing the 
# naive cells (POT/pre-treatment) to the chosen concentration of drug. We call
# this x1. 
# Again, choose some arbitrary, true value for 'reading x1'. 

x1 <- 1

# If we assume that the readings scale linearly between these two points - 
# as a function of the proportion of resistance in any given population - we can
# derive a linear equation for the reading value given the proportion of 
# resistance. 

# Lets do an example with b = 0.893, d = 0.200, mu = 1e-05, sig = 0.100

eqR <- stable_R_pheno(mu = 1e-05, sig = 0.100, del = 0.00, b = 0.893, d = 0.200)

m <- (x2 - x1)/(1 - eqR)
c <- x2 - (m * 1.0)

# Can now write a function that returns the underlying, 'true' reading (before
# measurement error): 

true_read <- function(p){
  x <- (m*p) + c
  return(x)
}

# For any given resistance proportion, we can now simulate a reading by taking
# this 'true_read' value as the mean, and adding measurement noise with 
# mes_err: 


# In reality our measurements of x1 and x2 are subject to measurement 
# error as well, so need to simulate the repeated sampling and then take the 
# mean as our x1 and x2 estimation: 

# Assume that readings ~Normal(mean = true_value, var = mes_err).

x1_reads <- rnorm(mes_n, mean = x1, sd = mes_err^2)
x2_reads <- rnorm(mes_n, mean = x2, sd = mes_err^2)

x1_est <- mean(x1_reads)
x2_est <- mean(x2_reads)

# And now we can write a function that returns our simulated, estimated 
# 'assay reading', given a proportion of resistance (p) and some measurement
# error, mes_err: 

ass_read <- function(x1_est, x2_est, mes_err, b, d, mu, sig, del, p){
  eqR <- stable_R_pheno(mu = mu, sig = sig, del = del, b = b, d = d)
  m_est <- (x2_est - x1_est)/(1 - eqR)
  c_est <- x2_est - (m_est * 1.0)
  x_read <- rnorm(1, mean = ((m_est*p) + c_est), sd = mes_err^2)
  return(x_read)
}

# Vectorize the function

ass_read <- Vectorize(ass_read)

# Can run this for all the simulated proportions of resistance: 

res_df["ass_read"] <- ass_read(rep(x1_est, nrow(res_df)), 
                               rep(x2_est, nrow(res_df)), 
                               rep(mes_err, nrow(res_df)),
                               res_df$b, res_df$d, anac(res_df$mu), 
                               anac(res_df$sig), anac(res_df$del), 
                               res_df$res_prop)

# Can now look at the distribution of recorded values vs the resistant 
# proportion. 

# ggplot(data = res_df, aes(x = delt, y = ass_read, fill = sig)) +
#   geom_point(size = 4, shape = 21, alpha = 0.6, colour = "black") +
#   facet_wrap(~sig, nrow = 1) +
#   theme_FW(12) +
#   scale_fill_manual(values = pal)
# 
# ggplot(data = res_df, aes(x = ass_read, group = delt, fill = sig)) +
#   geom_density(alpha = 0.6, colour = "black") +
#   facet_wrap(~sig, ncol = 1) +
#   theme_FW(12) +
#   scale_fill_manual(values = pal)

# Now have the output of a full, simulated experiment: we allow exclusively 
# resistant cells to evolve in the absence of treatment for some known time
# windows according to the parameters that control the evolution of the 
# resistance parameter. 

# Compare the relationship between mean vs variance after simulating the assay
# readings...

ass_summ_df <- plyr::ddply(res_df, .(delt, mu, sig, Nsim), summarise, 
                           mean = mean(ass_read), var = var(ass_read))

# NB. THAT THIS IS FOR ALL SIMULATIONS - SO MORE DATA THAN I'D ACTUALLY HAVE...

p1 <- ggplot(data = summ_df_rp, aes(x = mean, y = var, fill = sig)) + 
  geom_point(shape = 21, colour = "black", size = 4) + 
  geom_point(aes(x = exp_pR, y = var_pR)) +
  facet_grid(~mu)

p2 <- ggplot(data = ass_summ_df, aes(x = mean, y = var, fill = sig)) + 
  geom_point(shape = 21, colour = "black", size = 4) + 
  facet_grid(~mu)

cowplot::plot_grid(p1, p2, ncol = 1)



ggplot(data = subset(ass_summ_df, sig != 1e-05), aes(x = mean, y = var, fill = sig)) + 
  geom_point(shape = 21, colour = "black", size = 4) + 
  facet_grid(delt~mu)


###################################

# CALCULATING THE COMPOUND VARIANCE: 

# I think I'm going to set a minimum variance based on these theoretical 
# predictions... and then fit 2x extra components of variance: 
# 1. Some 'unknown' variance that is the (analytically intractable) variance 
#    introduced by the birth-death model. 
# 2. The measurement error noise - this should be constant across all replicates

# Think i have to calculate the compound variance for a reading, given the 
# distribution of the proportion of resistance... 

# Lets pick a value and see if this works. 

res_df_sub <- subset(res_df, mu == 1e-04 & sig == 0.1 & delt == 6)

# So the theoretical expectations - 
# mean proportion of resistance: 


mRp <- res_prop_t(mu = 1e-04, sig = 0.1, 
                  del = 0.0, br = 0.893, dr = 0.200, t = 6)

mRn <- mRp * unique(res_df_sub$Ni)

# And variance:

vRp <- mRp * (1 - mRp) / unique(res_df_sub$Ni)

vRn <- mRp * (1 - mRp) * unique(res_df_sub$Ni)

# Now, let say we're going to take 'mes_n' * readings using K cells per reading.
# Well now we're performing mes_n rounds of sampling 
#....
# i've already done the sampling for the assay... thats what K corresponds to 
# in the dataframe... 

# Lets do an assay reading where the additional measurement error = 0.0...

res_df_sub["ass_read_0e"] <- ass_read(x1, x2, mes_err = 0.00,
                                      b = res_df_sub$b, d = res_df_sub$d, 
                                      mu = 1e-04, sig = 0.1, del = 0.0, 
                                      p = res_df_sub$res_prop)

ggplot(data = res_df_sub, aes(x = res_prop, y = ass_read_0e)) + 
  geom_point()


eq_R <- stable_R_pheno(1e-04, 0.1, 0.0, res_df_sub$b, res_df_sub$d)

# First get the gradient and the intercept of the transform...
m <- (x2 - x1)/(1 - eq_R) # CAN SEE THAT AS eq_R INCREASES, SO DOES THE GRADIENT
c <- x2 - (1 * m)

# So the mean of the transformed variable: 
t_mRp <- m * mRp + c

t_mRp
mean(res_df_sub$ass_read_0e)

# ... and the variance: 
t_vRp <- m^2 * vRp

t_vRp
var(res_df_sub$ass_read_0e)

# So these match up - means I now have a way to translate the minimum variance
# i expect (due to the binomial sampling) into a minimum variance for the assay
# readings... I can then fit this minimum variance, and an additional 'unknown' 
# birth-death variance, and then finally a fixed measurement error variance. 

# I also think that I have additional information, because the transformed
# variance is contingent on the stable equilibrium frequency (because m is 
# a function of eq_R)... so after the transformation due to the assay readings, 
# I think I might even have more information from the variance than before?...

# Come back to this... don't think they all match up in the jupyter notebook: 

########################
# FROM JUPYTER NOTEBOOK:
########################

res_df <- bind_rows(res_dfs)
res_df$mu <- as.factor(res_df$mu)
res_df$sig <- as.factor(res_df$sig)

ass_read <- function(x1, x2, mes_err, b, d, mu, sig, del, p){
  eqR <- stable_R_pheno(mu = mu, sig = sig, del = del, b = b, d = d)
  # Gradient of linear relationship. 
  m <- (x2 - x1)/(1 - eqR)
  # Intercept of linear relationship. 
  c <- x2 - (m * 1.0)
  x_read <- rnorm(1, mean = ((m*p) + c), sd = mes_err^2)
  return(x_read)
}

# Vectorize the function
ass_read <- Vectorize(ass_read)

res_df["ass_read"] <- ass_read(rep(x1, nrow(res_df)), 
                               rep(x2, nrow(res_df)), 
                               rep(mes_err, nrow(res_df)),
                               res_df$b, res_df$d, anac(res_df$mu), 
                               anac(res_df$sig), anac(res_df$del), 
                               res_df$res_prop)

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

# First lets collect the mean and variances from our simulations: 

summ_df <- plyr::ddply(res_df, .(mu, sig, delt, Nsim), summarise, 
                       mean_nR = mean(res_n), var_nR = var(res_n),
                       mean_pR = mean(res_prop), var_pR = var(res_prop))

head(summ_df)

# Add a column that is the expected proportion of resistance given these parameters:

summ_df["th_mean_nR"] <- res_prop_t(anac(summ_df$mu), anac(summ_df$sig), 0.0, unique(res_df$b), unique(res_df$d), summ_df$delt)*unique(res_df$Ni)

summ_df["th_mean_pR"] <- res_prop_t(anac(summ_df$mu), anac(summ_df$sig), 0.0, unique(res_df$b), unique(res_df$d), summ_df$delt)

head(summ_df)

# Use this as the probability of success for the expected variance of resistance given these parameters: 

# binomial p
bin_ps <- summ_df$th_mean_nR/unique(res_df$Ni)
# binomial q
bin_qs <- (1 - bin_ps)
# binomial n 
bin_ns <- unique(res_df$Ni)

summ_df["th_var_nR"] <- bin_ps * bin_qs * bin_ns

summ_df["th_var_pR"] <- bin_ps * bin_qs / unique(res_df$Ni)

head(summ_df)

# A function to return the theoretical exepcation of variance: 

ass_exp <- function(x1, x2, mu, sig, del, b, d, exp_pRt){
  
  eq_R <- stable_R_pheno(mu, sig, del, b, d)
  m <- (x2 - x1)/(1 - eq_R)
  c <- x2 - (1 * m)
  exp_x <- (exp_pRt * m) + c
  return(exp_x)
}

ass_exp <- Vectorize(ass_exp)

ass_var <- function(x1, x2, mu, sig, del, b, d, var_pRt){
  
  eq_R <- stable_R_pheno(mu, sig, del, b, d)
  m <- (x2 - x1)/(1 - eq_R)
  c <- x2 - (1 * m)
  var_x <- m^2 * var_pRt
  return(var_x)
}

ass_var <- Vectorize(ass_var)

# Create a summary dataframe for the simulated assay readings

summ_ass_df <- plyr::ddply(res_df, .(mu, sig, delt, Nsim), summarise, mean_ass = mean(ass_read), var_ass = var(ass_read))

# Include the expectations and variances from the replicate sub-populations,. 

summ_ass_df["th_mean_pR"] <- res_prop_t(anac(summ_ass_df$mu), anac(summ_ass_df$sig), 0.0, unique(res_df$b), unique(res_df$d), summ_df$delt)

bin_ps <- summ_ass_df$th_mean_pR
bin_qs <- (1 - bin_ps)
bin_ns <- unique(res_df$Ni)
summ_ass_df["th_var_pR"] <- bin_ps * bin_qs / unique(res_df$Ni)

head(summ_ass_df)

# Add the theoretical mean and variance of the assay readings: 

summ_ass_df["th_mean_ass"] <- ass_exp(x1, x2, anac(summ_ass_df$mu), anac(summ_ass_df$sig), 0.0, unique(res_df$b), unique(res_df$d), summ_ass_df$th_mean_pR)

summ_ass_df["th_var_ass"] <- ass_var(x1, x2, anac(summ_ass_df$mu), anac(summ_ass_df$sig), 0.0, unique(res_df$b), unique(res_df$d), summ_ass_df$th_var_pR)

head(summ_ass_df)

ggplot(data = summ_ass_df, aes(x = mean_ass, y = var_ass, fill = sig)) + 
  geom_point(shape = 21, alpha = 0.6, colour = "black", size = 4) + 
  geom_point(aes(x = th_mean_ass, y = th_var_ass)) + 
  facet_grid(delt~mu) +
  theme_FW(14) +
  scale_fill_manual(values = pal)


ggplot(data = subset(summ_ass_df, mu == 1e-06), 
       aes(x = mean_ass, y = var_ass, fill = sig)) + 
  geom_point(shape = 21, alpha = 0.6, colour = "black", size = 4) + 
  geom_point(aes(x = th_mean_ass, y = th_var_ass)) + 
  facet_grid(~mu) +
  theme_FW(14) +
  scale_fill_manual(values = pal)


  











################################################################################

# Testing the Bayesian Model:
# Model1. 

library(rstan)
library(bayesplot)

options(mc.cores = parallel::detectCores())

rstan_options(auto_write = TRUE)

setwd("../../STAN_Models/")

# First, trying Model1 - doesn't include the measurement error noise, just
# takes the true, resistant proportions and tries to infer mu and sigma...
# think because of how i've coded up the model - there is no likelihood 
# function for the ODE solution - will just have to set the measurement 
# variance to ~0... 

res_df_sub <- subset(res_df, mu == 1e-05 & sig == 1e-01 & Nsim == 1)

b <- unique(res_df$b)
d <- unique(res_df$d)

nT <- anac(unique(res_df_sub$n_delt))
K <- nrow(res_df_sub)/nT
del_ts <- unique(res_df_sub$delt)

p_del_t1 <- subset(res_df_sub, delt == del_ts[[1]])$res_prop
p_del_t2 <- subset(res_df_sub, delt == del_ts[[2]])$res_prop
p_del_t3 <- subset(res_df_sub, delt == del_ts[[3]])$res_prop 

# TRY FITTING A DIFFERENT ERROR (epsi) FOR EACH del_t... 

plyr::ddply(res_df_sub, .(delt), function(x){sd(x$res_prop)})

model_dat <- list(b = b, d = d, nT = nT, del_ts = del_ts, K = K,
                  p_del_t1 = p_del_t1, p_del_t2 = p_del_t2, p_del_t3 = p_del_t3)


model_fit <- stan(file = "CBC_Rep_Exc_Res_Model1.stan",
                               data = model_dat, chains = 4, iter = 2000)



mcmc_dens_overlay(model_fit, pars = c("mu", "sigma"))

mcmc_pairs(model_fit, pars = c("mu", "sigma"))

mcmc_trace(model_fit, pars = c("mu", "sigma"))

# This version works. Now try implementing the measurement error as per the 
# method outlined above. 

###############################################################################

# Measurement Noise Version: 
# Model2.

# Set a value for the measurement error = the STD.DEV of the normal. 
mes_err <- 0.00
# Set the number of measurements equal to the number of replicates in the 
# simulated values. 
mes_n <- K
# Choose values for x1 and x2 'true' readings. 
x2 <- 10
x1 <- 1
# Assume that readings ~Normal(mean = true_value, sd = mes_err).
x1_reads <- rnorm(mes_n, mean = x1, sd = mes_err)
x2_reads <- rnorm(mes_n, mean = x2, sd = mes_err)
# We assume that x1 corresponds to pR = 1.0, and 
#                x2 corresponds to pR = eqR | b, d, mu, sig
# We also assume that the assay reading is a linear function of the proportion
# of resistance such that x = m(pR) + c. 
# We use this assumption to now simulate assay reading values for the simulated
# proportions of resistance. 
# n.b. use the 'true' x1 and x2 readings to simulated this relationship, whilst
# the estimated 'noisy' values are used to infer the relationship in the model. 
ass_read <- function(x1, x2, mes_err, b, d, mu, sig, del, p){
  eqR <- stable_R_pheno(mu = mu, sig = sig, del = del, b = b, d = d)
  m_est <- (x2 - x1)/(1 - eqR)
  c_est <- x2 - (m_est * 1.0)
  x_read <- rnorm(1, mean = ((m_est*p) + c_est), sd = mes_err)
  return(x_read)
}

ass_read <- Vectorize(ass_read)

res_df_sub <- subset(res_df, mu == 1e-05 & sig == 1e-01 & Nsim == 1)
res_df_sub <- res_df_sub[, !(names(res_df_sub) %in% "ass_read")]

res_df_sub["ass_read"] <- ass_read(
                               rep(x1, nrow(res_df_sub)), 
                               rep(x2, nrow(res_df_sub)), 
                               rep(mes_err, nrow(res_df_sub)),
                               res_df_sub$b, res_df_sub$d, anac(res_df_sub$mu), 
                               anac(res_df_sub$sig), anac(res_df_sub$del), 
                               res_df_sub$res_prop
                               )

# Extract simulated readings. 

x_del_t1 <- subset(res_df_sub, delt == del_ts[[1]])$ass_read
x_del_t2 <- subset(res_df_sub, delt == del_ts[[2]])$ass_read
x_del_t3 <- subset(res_df_sub, delt == del_ts[[3]])$ass_read

x1_m <- mean(x_del_t1)
x2_m <- mean(x_del_t2)
x3_m <- mean(x_del_t3)

plyr::ddply(res_df_sub, .(delt), function(x){sd(x$ass_read)})

ggplot(data = res_df_sub, aes(x = ass_read, fill = as.factor(delt))) + 
  geom_density(alpha = 0.6, colour = "black")

# Now, i think because we've included noise that is >> than the variance 
# introduced by the birth-death process alone, it makes sense to only include
# the (unknown) measurement noise to begin with, and then assess the model fit. 

# model_dat <- list(b = b, d = d, nT = nT, del_ts = del_ts, K = K,
#                   x_del_t1 = x_del_t1, x_del_t2 = x_del_t2, x_del_t3 = x_del_t3,
#                   x1_reads = x1_reads, x2_reads = x2_reads,
#                   x1_m = x1_m, x2_m = x2_m, xest_sd = 1.0)

model_dat <- list(b = b, d = d, nT = nT, del_ts = del_ts, K = K,
                  x_del_t1 = x_del_t1, x_del_t2 = x_del_t2, x_del_t3 = x_del_t3,
                  x1_est = x1, x2_est = x2)

model_fit <- stan(file = "CBC_Rep_Exc_Res_Model2.stan", 
                  data = model_dat, chains = 4, iter = 2000)


mcmc_dens_overlay(model_fit, pars = c("mu", "sigma"))

mcmc_pairs(model_fit, pars = c("mu", "sigma"))

mcmc_trace(model_fit, pars = c("mu", "sigma"))

#mcmc_dens_overlay(model_fit, pars = c("m_neg_exp", "s_neg_exp"))

# mcmc_dens_overlay(model_fit, pars = c("mu", "sigma", 
#                                       "mes_err1", "mes_err2", "mes_err3"))
# 
# 
# mcmc_trace(model_fit, pars = c("mu", "sigma", 
#                                       "mes_err1", "mes_err2", "mes_err3"))


ggplot(data = res_df_sub, aes(x = ass_read, fill = as.factor(delt))) + 
  geom_density(alpha = 0.6, colour = "black")

# TROUBLESHOOTING BY WRITING OUT THE BAYESIAN STEPS....

# Write in the custom functions. 

res_eq_freq <- function(mu, sigma){
  
  eq_R <- mu/(mu + sigma)
  
  return(eq_R)
  
}

res_prop_t <- function(b, d, del_t, mu, sigma){
  
  p_t <- (1 - (mu/(mu + sigma)))*exp((-2 * mu * b * del_t) - (2 * sigma * b * del_t)) + mu/(mu + sigma);
  
  return(p_t)
  
}

x_read_t <- function(x1_est, x2_est, mu, sigma, pR){
  
  eq_R <- res_eq_freq(mu, sigma)
  m_est <- (x2_est - x1_est)/(1 - eq_R)
  c_est <- x2_est - (m_est * 1.0)
  x_t <- (m_est*pR) + c_est
  
  return(x_t)
  
}

# 1. lets draw some values of mu, sig and err_mes...

samp_mu <- 0.075
samp_sig <- 0.02
samp_mes_err <- 0.01

# 2. Calculate the expected proportion of resistance at each time window...

exp_p_del_t1 <- res_prop_t(b, d, del_ts[1], samp_mu, samp_sig)
exp_p_del_t2 <- res_prop_t(b, d, del_ts[2], samp_mu, samp_sig)
exp_p_del_t3 <- res_prop_t(b, d, del_ts[3], samp_mu, samp_sig)

# 3. Calculate the expected assay readings given the linear relationship, 
# and the expected proportions of resistance... use the simulated, 'true' 
# values of 

exp_x_del_t1 <- x_read_t(x1, x2, samp_mu, samp_sig, exp_p_del_t1)
exp_x_del_t2 <- x_read_t(x1, x2, samp_mu, samp_sig, exp_p_del_t2)
exp_x_del_t3 <- x_read_t(x1, x2, samp_mu, samp_sig, exp_p_del_t3)


# 4. Simulate what the distributions would look like given these values and 
# plot them...

sim_t1 <- rnorm(K, exp_x_del_t1, samp_mes_err)
sim_t2 <- rnorm(K, exp_x_del_t2, samp_mes_err)
sim_t3 <- rnorm(K, exp_x_del_t3, samp_mes_err)

sim_df <- data.frame(x = c(sim_t1, sim_t2, sim_t3), 
                     delt = c(rep(6, K), rep(12, K), rep(18, K)))

# 5.Compare the simulated distributions given these values vs the true, 
# simulated assay reading values...

p1 <- ggplot(data = res_df_sub, aes(x = ass_read, fill = as.factor(delt))) + 
  geom_density(alpha = 0.6, colour = "black") + 
  scale_x_continuous(limits = c(0, max(sim_df$x+1))) 

p2 <- ggplot(data = sim_df, aes(x = x, fill = as.factor(delt))) + 
  geom_density(alpha = 0.6, colour = "black") + 
  scale_x_continuous(limits = c(0, max(sim_df$x+1))) + 
  ggtitle(paste0("mu = ", samp_mu,"; sig = ", samp_sig))

cowplot::plot_grid(p1, p2, ncol = 1)

stable_R_pheno(mu = 1e-05, sig = 1e-01, del = 0.0, b, d)
stable_R_pheno(mu = 0.075, sig = 0.02, del = 0.0, b, d)

res_prop_t(b, d, 12.0, mu = 1e-05, sig = 1e-01)
res_prop_t(b, d, 12.0, mu = 0.075, sig = 0.02)

# So these both lead to similar expectations from the assay despite different
# stable resistance frequencies and resistant proportions at time t...
# I think actually the high difference in equilibrium frequency is the reason 
# the model can't distinguish between these two... 

# Because x1 is calibrated against a 'sensitive' population where its assumed
# that the equilibrium frequency is determined by mu and sigma... so the two 
# scenarios the model can't distinguish: 

# i) there is a low equilibrium frequency of resistance and high sigma (for 
# example), and the population quickly approaches this 
# ii) ...

###//

# stable_R_pheno <- Vectorize(stable_R_pheno)
# 
# mus <- 10^-(seq(0, 10))
# mus <- rep(mus, each = length(mus))
# sigs <- 10^-(seq(0, 10))
# sigs <- rep(sigs, length(sigs))
# 
# eqRs <- stable_R_pheno(mus, sigs, 0.0, b, d)
# 
# eq_df <- data.frame(mu = mus, sig = sigs, eqR = eqRs)
# 
# ggplot(data = eq_df, aes(x = as.factor(mu), y = as.factor(sig), 
#                          fill = eqR)) + 
#   geom_tile() + 
#   scale_fill_viridis_c(trans = "log10")

###

###############################################################################

# Measurement Noise Version: 
# Model3.

library(rstan)
library(bayesplot)

options(mc.cores = parallel::detectCores())

rstan_options(auto_write = TRUE)

setwd("../../STAN_Models/")


# New model: this version explicitly includes 3 components of variance: 

#   i) the known, minimum variance introduced due to the phenotypic switching
#      given an initial condition of exclusive resistance - this is converted
#      into the minimum variance expected in the assay (currently assuming that
#      there is a linear relationship between the proportion of resistance and
#      the assay reading). 

#  ii) the unknown variance introduced due to the birth-death process. This is 
#      analytically intractable, but we do know that it must be _additional_ 
#      variance and therefore this is fit on top of the variance from i). 

# iii) the variance introduced when taking the assay readings. This is fixed
#      amongst all experimental readings, and is assumed to be an intrinsic
#      property of the assay. 

# Choose a parameter combination for inference. 

res_df_sub <- subset(res_df, mu == 1e-05 & sig == 1e-01 & Nsim == 1)
res_df_sub <- res_df_sub[, !(names(res_df_sub) %in% "ass_read")]

mes_err <- 0.00
mes_n <- unique(res_df_sub$K)
K <- unique(res_df_sub$K)
Ni <- unique(res_df_sub$Ni)
del_ts <- unique(res_df_sub$delt)
nT <- length(del_ts)
b <- unique(res_df_sub$b)
d <- unique(res_df_sub$d)

# Choose values for x1 and x2 'true' readings. 
x2 <- 10
x1 <- 1
# Assume that readings ~Normal(mean = true_value, sd = mes_err).
x1_reads <- rnorm(mes_n, mean = x1, sd = mes_err)
x2_reads <- rnorm(mes_n, mean = x2, sd = mes_err)

# We assume that x1 corresponds to pR = 1.0, and 
#                x2 corresponds to pR = eqR | b, d, mu, sig
# We also assume that the assay reading is a linear function of the proportion
# of resistance such that x = m(pR) + c. 

ass_read <- function(x1, x2, mes_err, b, d, mu, sig, del, p){
  eqR <- stable_R_pheno(mu = mu, sig = sig, del = del, b = b, d = d)
  m_est <- (x2 - x1)/(1 - eqR)
  c_est <- x2 - (m_est * 1.0)
  x_read <- rnorm(1, mean = ((m_est*p) + c_est), sd = mes_err)
  return(x_read)
}

ass_read <- Vectorize(ass_read)

res_df_sub["ass_read"] <- ass_read(
  rep(x1, nrow(res_df_sub)), 
  rep(x2, nrow(res_df_sub)), 
  rep(mes_err, nrow(res_df_sub)),
  res_df_sub$b, res_df_sub$d, anac(res_df_sub$mu), 
  anac(res_df_sub$sig), anac(res_df_sub$del), 
  res_df_sub$res_prop
)

# Extract simulated readings. 

x_del_t1 <- subset(res_df_sub, delt == del_ts[[1]])$ass_read
x_del_t2 <- subset(res_df_sub, delt == del_ts[[2]])$ass_read
x_del_t3 <- subset(res_df_sub, delt == del_ts[[3]])$ass_read

# FOR TROUBLESHOOTING:

pR_del_t1 <- subset(res_df_sub, delt == del_ts[[1]])$res_prop
pR_del_t2 <- subset(res_df_sub, delt == del_ts[[2]])$res_prop
pR_del_t3 <- subset(res_df_sub, delt == del_ts[[3]])$res_prop

# Create the model data object. 

# model_dat <- list(b = b, d = d, nT = nT, del_ts = del_ts, K = K, Ni = Ni,
#                   x_del_t1 = x_del_t1, x_del_t2 = x_del_t2, x_del_t3 = x_del_t3,
#                   x1s = x1_reads, x2s = x2_reads)
# model_dat <- list(b = b, d = d, nT = nT, del_ts = del_ts, K = K, Ni = Ni,
#                   x_del_t1 = x_del_t1, x_del_t2 = x_del_t2, x_del_t3 = x_del_t3,
#                   x1_est = x1, x2_est = x2)
model_dat <- list(b = b, d = d, nT = nT, del_ts = del_ts, K = K, Ni = Ni,
                  pR_del_t1 = pR_del_t1, pR_del_t2 = pR_del_t2, pR_del_t3 = pR_del_t3,
                  x1_est = x1, x2_est = x2)

model_fit <- stan(file = "CBC_Rep_Exc_Res_Model3.stan", 
                  data = model_dat, chains = 4, iter = 2000)


mcmc_dens_overlay(model_fit, pars = c("mu", "sigma"))

mcmc_pairs(model_fit, pars = c("mu", "sigma"))

mcmc_trace(model_fit, pars = c("mu", "sigma"))


##################

# TROUBLESHOOTING:

# FOLLOWING CAN MANUALLY CHECK MEAN V VARIANCE FOR GIVEN MU AND SIGMA VALUES...

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
ass_exp <- function(x1, x2, mu, sig, del, b, d, exp_pRt){
  
  eq_R <- stable_R_pheno(mu, sig, del, b, d)
  m <- (x2 - x1)/(1 - eq_R)
  c <- x2 - (1 * m)
  exp_x <- (exp_pRt * m) + c
  return(exp_x)
}
ass_exp <- Vectorize(ass_exp)
ass_var <- function(x1, x2, mu, sig, del, b, d, var_pRt){
  
  eq_R <- stable_R_pheno(mu, sig, del, b, d)
  m <- (x2 - x1)/(1 - eq_R)
  c <- x2 - (1 * m)
  var_x <- m^2 * var_pRt
  return(var_x)
}
ass_var <- Vectorize(ass_var)

comp_df <- data.frame(mu = rep(c(1e-05, 0.0830), 3), sig = rep(c(1e-01, 0.0165), 3), t = rep(c(6.0, 12.0, 18.0), each = 2),
                      param_set = rep(c(1, 2),3))
comp_df["pR"] <- res_prop_t(comp_df$mu, comp_df$sig, 0.0, 0.893, 0.200, comp_df$t)
comp_df["vR"] <- comp_df$pR * (1 - comp_df$pR) /unique(res_df$Ni)
comp_df["exp_ass"] <- ass_exp(x1, x2, comp_df$mu, comp_df$sig, 0.0, 0.893, 0.200, comp_df$pR)
comp_df["var_ass"] <- ass_var(x1, x2, comp_df$mu, comp_df$sig, 0.0, 0.893, 0.200, comp_df$vR)
comp_df$param_set <- as.factor(comp_df$param_set)
ggplot(data = comp_df, aes(x = exp_ass, y = var_ass, fill = param_set)) + 
  geom_point(shape = 21, alpha = 0.6, colour = "black", size = 4) + 
  theme_FW(14)


# LOOK MILES APART... MANUALLY WRITE IN THE STAN FUNCTIONS AND CHECK...

res_eq_freq<-function(mu, sigma){
  
  eq_R <- mu/(mu + sigma)
  
  return(eq_R)
  
}

res_eq_freq <- Vectorize(res_eq_freq)

exp_res_prop_t<-function(b, d, del_t, mu, sigma){
    
    exp_p_t <- (1 - (mu/(mu + sigma)))*exp((-2 * mu * b * del_t) - (2 * sigma * b * del_t)) + mu/(mu + sigma)
    
    return(exp_p_t)
    
}

exp_res_prop_t <- Vectorize(exp_res_prop_t)
  
var_res_prop_t<-function(b, d, del_t, mu, sigma, Ni){
    
    exp_p_t <- exp_res_prop_t(b, d, del_t, mu, sigma)
    var_p_t <- ((exp_p_t) * (1 - exp_p_t))/Ni
    
    return(var_p_t)
}

var_res_prop_t <- Vectorize(var_res_prop_t)

exp_x_read_t<-function(x1_est, x2_est, mu, sigma, exp_pR){
  
  eq_R <- res_eq_freq(mu, sigma)
  m_est <- (x2_est - x1_est)/(1 - eq_R)
  c_est <- x2_est - (m_est * 1.0)
  exp_x_t <- (m_est*exp_pR) + c_est
  
  return(exp_x_t)
  
}

exp_x_read_t <- Vectorize(exp_x_read_t)

var_x_read_t<-function(x1_est, x2_est, mu, sigma, var_pR){
  
  eq_R <- res_eq_freq(mu, sigma);
  m_est <- (x2_est - x1_est)/(1 - eq_R);
  c_est <- x2_est - (m_est * 1.0);
  var_x_t <- m_est^2 * var_pR;
  
  return(var_x_t)
  
}

var_x_read_t <- Vectorize(var_x_read_t)

exp_pR <- exp_res_prop_t(0.893, 0.200, c(6.0, 12.0, 18.0), 1e-05, 1e-01)
var_pR <- var_res_prop_t(0.893, 0.200, c(6.0, 12.0, 18.0), 1e-05, 1e-01, 10000)

tp1 <- rnorm(10000, mean = exp_pR[1], sd = sqrt(var_pR[1]))
tp2 <- rnorm(10000, mean = exp_pR[2], sd = sqrt(var_pR[2]))
tp3 <- rnorm(10000, mean = exp_pR[3], sd = sqrt(var_pR[3]))

tpdf_1 <- data.frame(x = c(tp1, tp2, tp3), t = rep(c(1, 2, 3), each = 10000))

exp_pR <- exp_res_prop_t(0.893, 0.200, c(6.0, 12.0, 18.0), 0.85, 0.26)
var_pR <- var_res_prop_t(0.893, 0.200, c(6.0, 12.0, 18.0), 0.85, 0.26, 10000)



exp_xs <- exp_x_read_t(x1, x2, 1e-05, 1e-01, exp_pR)
var_xs <- var_x_read_t(x1, x2, 1e-05, 1e-01, var_pR)

# So looks like the initial variance calculation function wasn't right... STILL
# NEED TO DOUBLE CHECK WHY THE SCALING WORKS LIKE THIS... MAYBE CHECK BY 
# RUNNING A RANDOM Ni... 
