
# Cancer Barcode Simulation
# Freddie Whiting - 2021

# Plotting the qD statistics extracted using the CBC_Exp_Output_Summ_Stats 
# script. The qD statistics capture the within-replicate diversity (qD) and 
# the between-replicate diversity dissimilarity (qD_diss). 

################################################################################

rm(list = ls())

# Load in command line arguments: 

#args <- commandArgs(trailingOnly = TRUE)
#try(if(length(args) != 3) stop("There should be 3 command-line arguments."))


# TEMPORARY ####################################################################
setwd("C:/Users/whitin01/Google Drive/Barcode_Simulations/Cancer_BarCode_Sim_CBC/Outputs/")
#setwd("C:/Users/Freddie/Google Drive/Barcode_Simulations/Cancer_BarCode_Sim_CBC/Outputs/")
args <- c()
args[1] <- "CBC_Exp_Output_Summary_Stats_HCTbc_sim_del_psi_0.15.csv"
args[2] <- "../../Cancer_BarCode_Sim_CBC/CBC_Analysis_Custom_Functions.R"
args[3] <- "HCTbc_sim_del_psi_0.15"
args[4] <- "del"
args[5] <- "2"
args[6] <- T
args[7] <- "../../../Barcode_Analysis/QR_Experiment/QR_Summary_Stats/HCTbc_Output_Summary_Stats.csv"
  #args[7] <- "../../../Barcode_Analysis/QR_Experiment/QR_Summary_Stats/SW6bc_Output_Summary_Stats.csv"
################################################################################

# First argument - output summ_stats .csv. 
out_csv_path <- args[1]
# Second argument - custom function script. 
cust_fun_path <- args[2]
# Third argument - simulation set prefix for plotting directory.
simset_pref <- args[3]
# Fourth argument - resistance phenotype evolution parameter (controls resistance
# evolution in the simulations - either p, sig or del). 
res_param <- args[4]
# Fifth argument - passage to use for plotting the qD panel plots. 
passage <- args[5]

# Sixth argument - should a qD-dist plot be produced by comparing the 
# simulations to sequenced qD statistics? (bool).
qD_dist_plot <- args[6]
# Seventh argument - if qD_dist_plot == T, the name of the sequenced output 
# file with the analogous qD statistics. 
if(args[6] == T){
  seq_csv_path <- args[7]  
}

################################################################################
# Data Preparation 
################################################################################

# Load in the custom CBC functions. 
source(cust_fun_path)

# Read in the simulation summary statistics. 
summ_stats <- read.csv(out_csv_path, stringsAsFactors = F)

# Add a column which is the equilibrium frequency of resistance given the
# simulation parameter values - if looking at the pre-existing/de-novo outputs, 
# set the eq_R value = p (this is the proportion of resistance assigned at t=0).

if(res_param == "p"){
  summ_stats["eq_R"] <- summ_stats$p
} else {
  summ_stats["eq_R"] <- sapply(1:nrow(summ_stats), function(x){stable_R_pheno(summ_stats$mu[[x]], summ_stats$sig[[x]], summ_stats$del[[x]], 
                                                                              summ_stats$b[[x]], summ_stats$d[[x]],
                                                                              sig_fig_p = T, n_dig=1)})
}

# Convert relevant variables to factors and then order for plotting. 

summ_stats <- order_factors(summ_stats)

# Just look at the chosen passage, and separate the control and drug-treatments. 

COsub_df <- subset(summ_stats, Passage == paste0("P", passage) & type == "CO")
DTsub_df <- subset(summ_stats, Passage == paste0("P", passage) & type == "DT")

# Also only now look at q=2.

COsub_df <- subset(COsub_df, q == 2)
DTsub_df <- subset(DTsub_df, q == 2)

# Collect the mean values of the control qD to add to all other plots as a reference point. 

CO_df_summ <- data.frame(qD = mean(COsub_df$qD), 
                         qD_diss = mean(COsub_df$qD_diss))

# Add this to the drug-treatment subset as columns for plotting. 

DTsub_df["CO_qD"] <- CO_df_summ$qD
DTsub_df["CO_qD_diss"] <- CO_df_summ$qD_diss

# Replace the name of the res param (e.g. p, sig, del) with the 'res_param' - 
# will make calling this variable column in the plotting functions easier. 

colnames(DTsub_df)[grep(paste0("\\<", res_param, "\\>"), 
                        colnames(DTsub_df))] <- "res_param"

# Write the named facets for plotting. 

named_ps <- paste0(" \U03C1 = \n", unique(as.character(summ_stats$p)))
names(named_ps) <- unique(as.character(summ_stats$p))
named_mus <- paste0(" \U03BC = \n", unique(as.character(summ_stats$mu)))
names(named_mus) <- unique(as.character(summ_stats$mu))
named_sigs <- paste0(" \U03C3 = \n", unique(as.character(summ_stats$sig)))
names(named_sigs) <- unique(as.character(summ_stats$sig))
named_dels <- paste0(" \U03B4 = \n", unique(as.character(summ_stats$del)))
names(named_dels) <- unique(as.character(summ_stats$del))
named_qs <- paste0(" q = \n", unique(as.character(summ_stats$q)))
names(named_qs) <- unique(as.character(summ_stats$q))
named_eq_Rs <- paste0(" Req = \n", unique(as.character(summ_stats$eq_R)))
names(named_eq_Rs) <- unique(as.character(summ_stats$eq_R))

# Create a dataframe of symbols used in plotting with corresponding character
# names. 

symbol_df <- data.frame(param = c("p", "mu", "sig", "del", "q", "eq_R"),
                        symbols = c("\U03C1", "\U03BC", "\U03C3", "\U03B4", "q", "Req"))

# Need to change labels according to the resistance parameter chosen: 
if(res_param == "p"){
  label_facets_nl <- labeller("res_param" = named_ps, "mu" = named_mus, "q" = named_qs, "eq_R" = named_eq_Rs)
}
if(res_param == "sig"){
  label_facets_nl <- labeller("mu" = named_mus, "res_param" = named_sigs, "q" = named_qs, "eq_R" = named_eq_Rs)
}
if(res_param == "del"){
  label_facets_nl <- labeller("mu" = named_mus, "res_param" = named_dels, "q" = named_qs, "eq_R" = named_eq_Rs)
}

# Get a summary of the number of completed simulation runs per parameter set

ncomp_sims <- plyr::ddply(DTsub_df, .(mu, res_param), summarise, n_sim = length(unique(sim_N)))

# Get the simulation parameter sets that completed < 18 simulation runs (90%)

ncomp_sims_sub <- subset(ncomp_sims, n_sim < 18)

# Remove these from the drug-treatment dataframe

ncomp_sims_sub <- ncomp_sims_sub[, c("mu", "res_param")]

to_remove <- merge(ncomp_sims_sub, DTsub_df)

DTsub_df <- anti_join(DTsub_df, to_remove)

################################################################################
# Plotting 
################################################################################

# Move to Plots Directory 

setwd("../Plots/")

# First qD plot: all parameters, colour points by resistant parameter.  

p1 <- ggplot(data = DTsub_df, aes(x = qD, y = qD_diss, fill = res_param)) + 
  geom_point(shape = 21, alpha = 0.6, colour = "black", size = 6, stroke = 1) + 
  scale_x_continuous(trans = "mylog10_2", breaks = mybreaks_2[c(seq(1, length(mybreaks_2), 2))], limits = c(0, 10^6)) + 
  scale_y_continuous(limits = c(1, 4)) + 
  facet_grid(mu~res_param, labeller = label_facets_nl) + 
  scale_fill_manual(values = pal, name = symbol_df[grep(res_param, symbol_df$param) ,]$symbols) +
  theme_FW(18) + 
  xlab(expression(""^{"q=2"}~D)) + 
  ylab(expression(""^{"q=2"}~D~({"\U03B2"}))) +
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5, hjust=0.5)) + 
  geom_point(data = DTsub_df, aes(x = CO_qD, y = CO_qD_diss), 
             shape = 8, fill = "black", colour = "black", size = 2, alpha = 0.4, stroke = 1.5)

ggsave(paste0(args[3], "_qD_All_Pairwise_Colour_Param.png"), p1, width = 11, height = 10)

# Second qD plot: all parameters, colour points by equilibrium resistance 
# frequency. 

p2 <- ggplot(data = DTsub_df, aes(x = qD, y = qD_diss, fill = anac(eq_R))) + 
  geom_point(shape = 21, alpha = 0.6, colour = "black", size = 6, stroke = 1) + 
  scale_x_continuous(trans = "mylog10_2", breaks = mybreaks_2[c(seq(1, length(mybreaks_2), 2))], limits = c(0, 10^6)) + 
  scale_y_continuous(limits = c(1, 4)) + 
  facet_grid(mu~res_param, labeller = label_facets_nl) + 
  scale_fill_viridis_c(name = expression(R["eq"])) + 
  theme_FW(18) + 
  xlab(expression(""^{"q=2"}~D)) + 
  ylab(expression(""^{"q=2"}~D~({"\U03B2"}))) +
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5, hjust=0.5)) + 
  geom_point(data = DTsub_df, aes(x = CO_qD, y = CO_qD_diss), 
             shape = 8, fill = "black", colour = "black", size = 2, alpha = 0.4, stroke = 1.5)

ggsave(paste0(args[3], "_qD_All_Pairwise_Colour_Req.png"), p2, width = 10.9, height = 10)

# Third qD plot: a subset of parameters with eq_R < ~0.1, and the panel 
# columns are now arranged by eq_R. 

p3 <- ggplot(data = subset(DTsub_df, anac(eq_R) <= 0.1), aes(x = qD, y = qD_diss, fill = res_param)) + 
  geom_point(shape = 21, alpha = 0.6, colour = "black", size = 6, stroke = 1) + 
  scale_x_continuous(trans = "mylog10_2", breaks = mybreaks_2[c(seq(1, length(mybreaks_2), 2))], limits = c(0, 10^6)) + 
  scale_y_continuous(limits = c(1, 4)) + 
  facet_grid(mu~eq_R, labeller = label_facets_nl) + 
  scale_fill_manual(values = pal, name = symbol_df[grep(res_param, symbol_df$param) ,]$symbols) +
  theme_FW(18) + 
  xlab(expression(""^{"q=2"}~D)) + 
  ylab(expression(""^{"q=2"}~D~({"\U03B2"}))) +
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5, hjust=0.5)) + 
  geom_point(data = subset(DTsub_df, anac(eq_R) <= 0.1), aes(x = CO_qD, y = CO_qD_diss), 
             shape = 8, fill = "black", colour = "black", size = 2, alpha = 0.4, stroke = 1.5)

ggsave(paste0(args[3], "_qD_Sub_Pairwise_Colour_Param.png"), p3, width = 10, height = 9)

# If the qD_dist_plot argument is set to TRUE, then proceed with a distance plot

if(qD_dist_plot == T){
  
  # Load in the sequenced qD summary statistic dataframe
  seq_summ_df <- read.csv(seq_csv_path, stringsAsFactors = F)
  # Subset the drug-treatment rows and only now look at q=2.
  seqsub_df <- subset(seq_summ_df, type == "DT" & q == 2)
  # Not including P5 - remove if present
  seqsub_df <- subset(seqsub_df, Passage != "P5")
  
  # Distinguish colnames with '_seq' suffix. 
  colnames(seqsub_df)[grep("qD|qD_diss", colnames(seqsub_df))] <- paste0(colnames(seqsub_df)[grep("qD|qD_diss", colnames(seqsub_df))], "_seq")

  # Re-subset the drug-treatment dataframe
  DTdist_df <- subset(summ_stats, type == "DT" & q == 2)
  
  # Join to the drug-treatment simulation dataframe
  full_dist_df <- join(DTdist_df, seqsub_df, type = "full")
  
  # Calculate the qD distances
  full_dist_df <- qD_distance(full_dist_df, "_seq", 4)
  
  # 
  colnames(full_dist_df)[grep(paste0("\\<", res_param, "\\>"), 
                          colnames(full_dist_df))] <- "res_param"
  
  # Get a summary of the number of completed simulation runs per parameter set
  ncomp_sims <- plyr::ddply(full_dist_df, .(mu, res_param), summarise, n_sim = length(unique(sim_N)))
  # Get the simulation parameter sets that completed < 18 simulation runs (90%)
  ncomp_sims_sub <- subset(ncomp_sims, n_sim < 18)
  # Remove these from the drug-treatment dataframe
  ncomp_sims_sub <- ncomp_sims_sub[, c("mu", "res_param")]
  to_remove <- merge(ncomp_sims_sub, full_dist_df)
  full_dist_df <- anti_join(full_dist_df, to_remove)
  
  # Make/change directory to a sub-dir given the sequenced values provided
  setwd(seq_summ_df$cell_line)
  
  # Fourth plot: heatmap of qD distances by each Passage
  
  p4 <- ggplot(full_dist_df, aes(x = res_param, y = mu, fill = qD_seq_dist)) + 
    geom_tile() + 
    scale_fill_viridis_c(direction = -1, limits = c(0.0, 0.8), name = "qD-Dist.",
                         option = "magma") + 
    ylab("\U03BC") +
    xlab(symbol_df[grep(res_param, symbol_df$param) ,]$symbols) +
    facet_wrap(~Passage, nrow = 1) +
    theme_minimal() + 
    theme(panel.background = element_rect(colour = "grey80"),
          text = element_text(size = 18), element_line(size = 0.0),
          panel.border = element_rect(colour = "black", fill=NA, size=3),
          panel.spacing = unit(-0.5, "lines"),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  
  ggsave(paste0(args[3], "_qD_Dist_All_Pass.jpeg"), p4, width = 13, height = 4)
  
  # Take a mean value of the qD distances across all Passages
  
  qD_dist_summ <- plyr::ddply(full_dist_df, 
                              .(mu, res_param),
                              summarise,
                              mean = mean(qD_seq_dist),
                              sd = sd(qD_seq_dist))
  

  # Fifth plot: heatmap of qD distances averaged across each Passage
  
  p5 <- ggplot(qD_dist_summ, aes(x = 1, y = 1, fill = mean)) + 
    geom_tile() + 
    scale_fill_viridis_c(direction = -1, limits = c(0.0, 0.8), name = "qD-Dist.",
                         option = "magma") + 
    ylab("\U03BC") +
    xlab(symbol_df[grep(res_param, symbol_df$param) ,]$symbols) +
    facet_grid(mu~res_param, labeller = label_facets_nl) +
    theme_minimal() + 
    theme(panel.background = element_rect(colour = "grey80"),
          text = element_text(size = 18), element_line(size = 0.0),
          panel.border = element_rect(colour = "black", fill=NA, size=3),
          panel.spacing = unit(-0.5, "lines"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) 
  
  ggsave(paste0(args[3], "_Avg_qD_Dist.jpeg"), p5, width = 7.5, height = 6)
  
  setwd("../")
  
}

setwd("../Outputs/")

################################################################################

