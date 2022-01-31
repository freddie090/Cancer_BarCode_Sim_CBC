
# Cancer Barcode Simulation
# Freddie Whiting - 2021

# Collection of custom functions that are used for extracting
# summary statistics from simulated barcode lineage distributions 
# as well as plotting the outputs. 

# Load required R packages: 

library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(scales)
library(reldist)
library(reshape2)
library(magrittr)
library(stringr)
library(vegan)
library(tidyr)
library(plyr)


#################################################################
# Summary Statistics Functions
#################################################################

# # Convert directory into dataframe of sim parameters. 

dir_to_params <- function(dir){
  param_df <- data.frame(
    sim_N   = sub(".*Sim_(.*)", replacement = "\\1", sim_dir),
    N       = sub(".*_N-([0-9]+)_.*", replacement = "\\1", dir),
    b       = sub(".*_b-([0-9]+\\.[0-9]+)_.*", replacement = "\\1", dir),
    d       = sub(".*_d-([0-9]+\\.[0-9]+)_.*", replacement = "\\1", dir),
    p       = sub(".*_p-(.*)_mu.*", replacement = "\\1", dir),
    mu      = sub(".*_mu-(.*)_sig.*", replacement = "\\1", dir),
    sig     = sub(".*_sig-(.*)_del.*", replacement = "\\1", dir),
    del     = sub(".*_del-(.*)_psi.*", replacement = "\\1", dir),
    psi     = sub(".*_psi-(.*)_al.*", replacement = "\\1", dir),
    R_real  = sub(".*_R_real-(\\D+)_n_pulse.*", replacement = "\\1", dir),
    n_pulse = sub(".*_n_pulse-([0-9]+)_.*", replacement = "\\1", dir),
    Nmax    = sub(".*_Nmax-([0-9]+)_.*", replacement = "\\1", dir),
    N_seed  = sub(".*_N_seed-([0-9]+)_.*", replacement = "\\1", dir),
    t_CO    = sub(".*_t_CO-([0-9]+\\.[0-9]+)_.*", replacement = "\\1", dir),
    t_DT    = sub(".*_t_DT-([0-9]+\\.[0-9]+).*", replacement = "\\1", dir),
    ik      = sub(".*_ik-([0-9]+).*", replacement = "\\1", dir),
    lp      = sub(".*_lp-([0-9]+)", replacement = "\\1", dir),
    stringsAsFactors = F
  )
  
  return(param_df)
}

# A function that adds '1's in a new column - cum_rf_gt - until the cumulative frequency 
# is equal to or greater than a value, set with 'limit' - inside the function.
# N.B. X <= 1.0. 

cum_rf_gt <- function(df, X){
  xvec <- rep(0, nrow(df))
  limit <- X
  for(i in seq_along(xvec)){
    if(df[i, "cum_rf"] >= limit){
      xvec[i] <- 1
      df["cum_rf_gt"] <- xvec
      return(df)
    } else {
      xvec[i] <- 1
    }
  }
}


# Function to return the relative frequency dataframe, but only for barcodes that make
# up the top (X*100)% of the sample in at least one sample. 

sub_gt_rfs <- function(df, X){
  
  # Track which barcodes are counted before the cumulative frequency of the sample > (X*100)%. 
  df["bc"] <- 1:nrow(df)
  df <- gather(df, key = sample, value = rel_freq, grep("bc", colnames(df), invert = T))
  # Order each from highest to lowest rel_freq, by sample. 
  df <- with(df, df[order(sample, rel_freq, decreasing = T) ,])    
  # Add a cumsum column 
  df <- plyr::ddply(df, .(sample), mutate, cum_rf = cumsum(rel_freq))
  # Add a 1 in a new column - cum_rf_gt - until reaching some threshold defined in 
  # the cum_rf_gt function. 
  df <- plyr::ddply(df, .(sample), cum_rf_gt, X=X)
  # remove any with a rel_freq of 0 
  df <- df[df$rel_freq != 0 ,]
  # remove all those that don't contribute to the top X %.
  df <- df[df$cum_rf_gt == 1 ,]
  # turn back to wide format. 
  df <- reshape(df[, c("bc", "sample", "rel_freq")], idvar = "bc", timevar = "sample", direction = "wide")
  # remove NAs 
  df[is.na(df)] <- 0
  # convert names back to original names
  colnames(df) <- sub("rel_freq\\.", "", colnames(df))
  # remove barcode column 
  df <- df[, grep("bc", colnames(df), invert = T)]
  # return dataframe. 
  return(df)
  
}

# Custom Diversity Index Function using Hill Numbers. 
# Exceptions for: 
# q = 1 - need to manually set to exp(shannon). 
# q = Inf - need to manually set to 1/(max_rf)

qD <- function(xs, q){
  xs <- xs[xs > 0]
  if(q == 1){
    x3 <- exp(diversity(xs, index = "shannon"))
  } else if (q == Inf){
    x3 <- 1/(max(xs))
  } else {
    x1 <- xs^q
    x2 <- sum(x1)
    x3 <- x2^(1/(1 - q))
  }
  return(x3)
}

# Calculate the qD dissimilarity, as in Jasinska
# et al. 2020 online methods. 
# = qD(mean) * qM(eff)
# Takes a dataframe where columns are samples, 
# and columns are relative frequencies. 

# qD(mean) = (1/M * sum(qD(p))^(1-q))^(1/1-q)
# where M = number of samples, and
# qD(p) is the qD using q=q for the p'th sample.

# qM(eff) = qD(pooled)/qD(mean)
#         = qD for all pooled lineages, given q/
#           qD mean across all p samples. 

# n.b. need exceptions again for q=1 and q=Inf.

qD_diss <- function(xs_df, q){
  
  if(q == 1){
    
    # Collect log(xi/mean(x)), where
    # mean(x) = mean rel_freq x across all i's. 
    lxi_mx <- lapply(1:ncol(xs_df), function(i){
      log(xs_df[, i]/apply(xs_df, 1, mean))    
    })
    # Now want sum(xi * log(xi/mean(x)))
    lxi_mxxi <- sapply(1:ncol(xs_df), function(i){
      vals <- xs_df[, i]*lxi_mx[[i]]
      # Convert NAs to 0s
      vals[is.na(vals)] <- 0
      return(vals)
    })
    # Take sum across all x lineages
    mx1 <- apply(lxi_mxxi, 2, sum)
    # And now mean across all i populations
    mx2 <- mean(mx1)
    # And finally e^ this value  
    qM <- exp(mx2)
    
  } else if(q == Inf){
    
    # Here we only look at the most abundant
    # lineages. 
    # We just want the number of 
    # i populations that the most abundant
    # lineages are found in. Therefore can
    # simply run: 
    max_xs <- xs_df[apply(xs_df, 2, which.max) ,]
    qM <- nrow(max_xs)
    
  } else {
    
    # Collect a qD | q for each sample. 
    qDs <- apply(xs_df, 2, qD, q=q)
    # Now all ^(1-q)
    qDs_pw <- qDs^(1 - q)
    # Now take the mean
    qDm <- mean(qDs_pw)
    # Now back-transform using ^1/(1-q)
    qDm_bt <- qDm^(1/(1 - q))
    # Calculate a qD for all pooled samples and
    # re-normalise so sum=1.0 for rel_freqs. 
    p_qDs <- rowSums(xs_df)/ncol(xs_df)
    # Now calculate pooled qD
    pqD <- qD(p_qDs, q=q)
    # Finally, want pooled/mean
    qM <- pqD/qDm_bt
    
  }
  
  return(qM)
}


# Save the following given a count dataframe where columns = 
# samples (including Passage information), and rows = barcode
# lineage, values = counts...

#   i) q-Div value (where q = Hill number). 
#  ii) The qD_dissimilarity calculated using qD_diss function. 

collect_summ_stats <- function(bc_counts){
  
  dfc <- bc_counts
  # remove any lineages not found in any of the samples. 
  dfc <- dfc[rowSums(dfc) > 0 ,]
  # create a relative frequency table
  dfr <- data.frame(apply(dfc, 2, prop.table))
  # First, if only one barcode, need to manually set dfr to all = 1.0.
  if(nrow(dfc) == 1){
    dfr <- dfc
    dfr[1 ,] <- 1.0
  }
    
  # ii) q-div 
  q0 <- data.frame(t(apply(dfr, 2, qD, q=0)))
  q0["q"] <- 0
  q1 <- data.frame(t(apply(dfr, 2, qD, q=1)))
  q1["q"] <- 1
  q2 <- data.frame(t(apply(dfr, 2, qD, q=2)))
  q2["q"] <- 2
  qI <- data.frame(t(apply(dfr, 2, qD, q=3)))
  qI["q"] <- Inf
  
  qdf <- bind_rows(q0, q1, q2, qI)
  qdf <- gather(qdf, key = "sample", value = "qD", grep("q", colnames(qdf), invert = T))
  
  dfr_stats <- qdf
  
  # Remove passage info and save as seperate column along with sample type. 
  
  dfr_stats["Passage"] <- sub("\\D{2}\\d{1}_(P\\d{1})", replacement = "\\1", dfr_stats$sample)
  
  dfr_stats$sample <- sub("(\\D{2}\\d{1})_P\\d{1}", replacement = "\\1", dfr_stats$sample)
  
  dfr_stats["type"] <- sub(".*(\\D{2})\\d{1}\\>", replacement = "\\1", dfr_stats$sample)
  
  if(length(unique(dfr_stats$Passage)) > 1){
    stop("dfr_stats should only have a single Passage variable.")
  }
  if(length(unique(dfr_stats$type)) > 1){
    stop("dfr_stats should only have a single type variable.")
  }
  
  # v) Calulate the qD_diss - Dissimilarity of Hill numbers - using the 4 chosenvalues of q. 
  qD_diss_q0 <- qD_diss(dfr, q = 0)
  qD_diss_q1 <- qD_diss(dfr, q = 1)
  qD_diss_q2 <- qD_diss(dfr, q = 2)
  qD_diss_qI <- qD_diss(dfr, q = Inf)
  # Join to eachother then stats. 
  qD_diss_df <- data.frame(q = c(0, 1, 2, Inf), qD_diss = c(qD_diss_q0, qD_diss_q1, qD_diss_q2, qD_diss_qI))
  
  dfr_stats <- join(dfr_stats, qD_diss_df, type = "full")
  
  return(dfr_stats)
  
}



# Calculate the normalised, transformed qD distance:
# - the qD distance (x-axis) is taken in log10 transformed space and has a 
#   maximum value of 10^6. 
# - the qD_diss distance (y-axis) is taken in non-log10 space.
# - both qD and qD_diss are normalised so that the qD distance can only fall
#   between 0 and 1 (therefore x- and y-axis limited to [0, ~0.7071])

# The function assumes that the dataframe has: 
#   i) a qD column 
#  ii) a qD_diss column 
# iii) qD and qD_disss columns with a given suffix (e.g. _HCTbc) with the 
#      sequencing statistics to be compared.
# n_sub_pop is the maximum number of sub-populations/replicates - used to 
# normalise the qD_diss values. 

qD_distance <- function(summ_df, samp_suffix, n_sub_pop){
  
  try(if(max(summ_df[, "qD"]) > 10^6) stop("qD should not exceed 10^6."))
  try(if(max(summ_df[, paste0("qD", samp_suffix)]) > 10^6) stop("qD should not exceed 10^6."))
  
  # Normalised simulations qDs
  x2s <- (log10(summ_df[, "qD"] + 1) / log10(10^6 + 1)) * 0.7071
  # Normalised sequencing qDs
  x1s <- (log10(summ_df[, paste0("qD", samp_suffix)] + 1) / log10(10^6 + 1)) * 0.7071
  # Normalised simulation qD_diss s
  y2s <- ((summ_df[, "qD_diss"] - 1)/(n_sub_pop - 1)) * 0.7071
  # Normalised sequencing qD_diss s
  y1s <- ((summ_df[, paste0("qD_diss", samp_suffix)] - 1)/(n_sub_pop - 1)) * 0.7071
  # Finished normalised distance: 
  norm_qD_dist <- sqrt((x2s - x1s)^2 + (y2s - y1s)^2)
  # Add to the dataframe
  summ_df[, paste0("qD", samp_suffix, "_dist")] <- norm_qD_dist
  
  return(summ_df)
    
}



#################################################################
# Other Functions.  
#################################################################

stable_R_pheno <- function(mu, sig, del, br, dr, sig_fig_p=F, n_dig=6){
  
  # Constructing Quadratic Formula
  result <- function(a,b,c){
    if(delta(a,b,c) > 0){ # first case D>0
      x_1 = (-b+sqrt(delta(a,b,c)))/(2*a)
      x_2 = (-b-sqrt(delta(a,b,c)))/(2*a)
      result = c(x_1,x_2)
    }
    else if(delta(a,b,c) == 0){ # second case D=0
      x = -b/(2*a)
    }
    else {"There are no real roots."} # third case D<0
  }
  
  # Constructing delta
  delta<-function(a,b,c){
    b^2-4*a*c
  }
  
  if(mu == 0.0 & sig == 0.0 & del == 0.0){
    p <- 0.0
  } else if(del == 0.0){
    p = mu/(mu + sig)
  } else {
    lam_r <- (br - dr)
    lam_R <- lam_r - (lam_r * del)
    # Assume splitting cost between birth and death rates.
    # i.e. R_real = "l"
    bR <- br * (1-del)
    
    # Convert into quadratic terms;
    a <- (-lam_R + lam_r)
    b <- (lam_R - (2*mu*br) - (2*sig*bR) - lam_r)
    c <- (2*mu*br)
    
    # Solve the quadratic;
    p <- result(a, b, c)[[2]]
  }

  if(sig_fig_p == T){
    return(signif(p, digits = n_dig))
  }else{
    return(p) 
  }
  
}


#################################################################
# Plotting. 
#################################################################


# Custom log transformation and axis breaks for range [0.0, 1.0].

mylog10_trans <- function (base = 10) {
  trans <- function(x) log(x + 1e-09, base)
  inv <- function(x) base^x
  trans_new(paste0("log-", format(base)), trans, inv, log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

mybreaks <- c(0, (10^(-10:0) - 1e-09))


# Custom log transformation and axis breaks for counts [0, 1, 2...].

mylog10_2_trans <- function (base = 10) {
  trans <- function(x) log(x + 1, base)
  inv <- function(x) base^x
  trans_new(paste0("log-", format(base)), trans, inv, log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

mybreaks_2 <- c(0, 10^(0:10) - 1)


# Custom theme for plotting within jupyter notebooks. 

theme_FW <- function(text_size) {
  library(grid)
  library(ggthemes)
  theme_minimal() + 
      theme(panel.background = element_rect(colour = "grey80"),
            text = element_text(size = text_size), element_line(size = 0.8),
            panel.border = element_rect(colour = "black", fill=NA, size=2)
            )
}


# Custom colour palettes. 

pal <- colorRampPalette(brewer.pal("Spectral", n = 9))(9)[c(1:3, 6:9)]

pal2 <- colorRampPalette(brewer.pal("Spectral", n = 10))(10)[c(1:3, 6:10)]

# Custom function to deal with plotting numerical factors as 
# continuous variables. 

anac <- function(x){return(as.numeric(as.character(x)))}

# Coerce variables to factors and then order for plotting. 

order_factors <- function(df){
  df$p <- factor(df$p, levels = unique(df$p)[order(unique(as.numeric(as.character(df$p))))])
  df$mu <- factor(df$mu, levels = unique(df$mu)[order(unique(as.numeric(as.character(df$mu))))])
  df$sig <- factor(df$sig, levels = unique(df$sig)[order(unique(as.numeric(as.character(df$sig))))])
  df$del <- factor(df$del, levels = unique(df$del)[order(unique(as.numeric(as.character(df$del))))])
  df$q <- factor(df$q, levels = unique(df$q)[order(unique(as.numeric(as.character(df$q))))])
  df$eq_R <- factor(df$eq_R, levels = unique(df$eq_R)[order(unique(as.numeric(as.character(df$eq_R))))])
  return(df)
}




# Function that, given a rel_freq_df (wide-fromat), two samples (x-axis, y-axis) 
# returns a ggplot object. Plot_n denotes which of the 4x4 (total 16) plots in 
# the grid this current function refers to. Can then remove axis elements to 
# ensure that the final plot isn't cluttered with repeated axis elements for 
# each individual plot. 

plot_comp_scatter <- function(x_samp, y_samp, plot_n, rf_df, ax_max, text_size,
                              abs_count=F){
  
  # Pull out chosen samples
  sub_df <- rf_df[, c("Center", x_samp, y_samp)]
  # Remove those barcodes found in neither
  sub_df <- sub_df[rowSums(sub_df == 0) < 2 ,]
  # Factorise by shared/unique
  sub_df["shared_unique"] <- 1 - (rowSums(sub_df == 0)) 
  sub_df$shared_unique[sub_df$shared_unique == 0] <- "unique"
  sub_df$shared_unique[sub_df$shared_unique == 1] <- "shared"
  sub_df$shared_unique <- as.factor(sub_df$shared_unique)
  
  # Plot params
  group.colours <- c(unique = "red", shared ="blue")
  
  if(abs_count == F){
    
    ax_breaks <- 10^(-10:0)
    # Plot object for relative frequencies
    p <- ggplot(data = sub_df, aes_string(x = x_samp, y = y_samp, fill = "shared_unique")) + 
      geom_point(alpha = 0.4, size = 3.0, shape = 21, colour = "black") + 
      scale_x_continuous(trans = "mylog10", breaks = c(0, (ax_breaks - 1e-09)), 
                         name = x_samp, limits = c(0, ax_max)) + 
      scale_y_continuous(trans = "mylog10", breaks = c(0, (ax_breaks - 1e-09)), 
                         name = y_samp, limits = c(0, ax_max)) + 
      scale_fill_manual(values = group.colours) +
      theme_FW(text_size) + 
      theme(legend.position = "none") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) + 
      annotation_logticks(base = 10, sides = "bl", outside = T)
  }

  
  if(abs_count == T){
    
    ax_breaks <- 10^(0:10)
    # Plot object for absolute counts
    p <- ggplot(data = sub_df, aes_string(x = x_samp, y = y_samp, fill = "shared_unique")) + 
      geom_point(alpha = 0.4, size = 3.0, shape = 21, colour = "black") + 
      scale_x_continuous(trans = "mylog10_2", breaks = c(0, (ax_breaks - 1)), 
                         name = paste0("\n", x_samp), limits = c(0, ax_max)) + 
      scale_y_continuous(trans = "mylog10_2", breaks = c(0, (ax_breaks - 1)), 
                         name = paste0(y_samp, "\n"), limits = c(0, ax_max)) + 
      scale_fill_manual(values = group.colours) +
      theme_FW(text_size) + 
      theme(legend.position = "none") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) + 
      annotation_logticks(base = 10, sides = "bl", outside = T)
  }
  
  # Remove xlabels if NOT in xlab vec
  xlab_vec <- c(13, 14, 15, 16)
  if(sum(xlab_vec == plot_n) == 0){
    p <- p + theme(axis.text.x=element_blank(),
                   axis.title.x=element_blank())
  }
  # Remove ylabels if NOT in ylab vec
  ylab_vec <- c(1, 5, 9, 13)
  if(sum(ylab_vec == plot_n) == 0){
    p <- p + theme(axis.text.y=element_blank(),
                   axis.title.y=element_blank())
  }
  
  return(p)
  
}






################################################################################

# Large colour palettes for muller plots etc. 

swatch <- function(x) {
  # x: a vector of colours (hex, numeric, or string)
  par(mai=c(0.2, max(strwidth(x, "inch") + 0.4, na.rm = TRUE), 0.2, 0.4))
  barplot(rep(1, length(x)), col=rev(x), space = 0.1, axes=FALSE,
          names.arg=rev(x), cex.names=0.8, horiz=T, las=1)
}

# Example:
# swatch(colours()[1:10])
# swatch(iwanthue(5))
# swatch(1:4)

iwanthue <- function(n, hmin=0, hmax=360, cmin=0, cmax=180, lmin=0, lmax=100,
                     plot=FALSE, random=FALSE) {
  # Presently doesn't allow hmax > hmin (H is circular)
  # n: number of colours
  # hmin: lower bound of hue (0-360)
  # hmax: upper bound of hue (0-360)
  # cmin: lower bound of chroma (0-180)
  # cmax: upper bound of chroma (0-180)
  # lmin: lower bound of luminance (0-100)
  # lmax: upper bound of luminance (0-100)
  # plot: plot a colour swatch?
  # random: should clustering be random? (if FALSE, seed will be set to 1,
  #         and the RNG state will be restored on exit.)
  require(colorspace)
  stopifnot(hmin >= 0, cmin >= 0, lmin >= 0,
            hmax <= 360, cmax <= 180, lmax <= 100,
            hmin <= hmax, cmin <= cmax, lmin <= lmax,
            n > 0)
  if(!random) {
    if (exists(".Random.seed", .GlobalEnv)) {
      old_seed <- .GlobalEnv$.Random.seed
      on.exit(.GlobalEnv$.Random.seed <- old_seed)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv))
    }
    set.seed(1)
  }
  lab <- LAB(as.matrix(expand.grid(seq(0, 100, 1),
                                   seq(-100, 100, 5),
                                   seq(-110, 100, 5))))
  if (any((hmin != 0 || cmin != 0 || lmin != 0 ||
           hmax != 360 || cmax != 180 || lmax != 100))) {
    hcl <- as(lab, 'polarLUV')
    hcl_coords <- coords(hcl)
    hcl <- hcl[which(hcl_coords[, 'H'] <= hmax & hcl_coords[, 'H'] >= hmin &
                       hcl_coords[, 'C'] <= cmax & hcl_coords[, 'C'] >= cmin &
                       hcl_coords[, 'L'] <= lmax & hcl_coords[, 'L'] >= lmin), ]
    #hcl <- hcl[-which(is.na(coords(hcl)[, 2]))]
    lab <- as(hcl, 'LAB')
  }
  lab <- lab[which(!is.na(hex(lab))), ]
  clus <- kmeans(coords(lab), n, iter.max=50)
  if (isTRUE(plot)) {
    swatch(hex(LAB(clus$centers)))
  }
  hex(LAB(clus$centers))
}

# From: https://gist.github.com/johnbaums/45b49da5e260a9fc1cd7
#################

###############################################################################



