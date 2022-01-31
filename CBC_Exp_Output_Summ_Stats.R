
# Cancer Barcode Simulation
# Freddie Whiting - 2021

# Collecting the Summary Statistics from the CBC_Exp simulations.   

################################################################################

rm(list = ls())

# Load in command line arguments: 

args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 4) stop("There should be 4 command-line arguments."))

# First argument - output directory
out_dir_path <- args[1]
# Second argument - run directory n (Int, > 0)
run_dir_n <- as.numeric(args[2])
# Second argument - custom function script. 
cust_fun_path <- args[3]
# Third argument - output-specific identifier prefix (e.g. 'sigs_psi_0.0_')
out_id <- args[4]

################################################################################

# Load in the custom CBC functions. 
source(cust_fun_path)

# Set directory to results directory. 
setwd(out_dir_path)
  
# Initiate empty master_df. 
master_df <- data.frame()
ts_count <- 1

# Extract all simulation runs 
run_dirs <- grep("CBC_Exp_out_", list.files(), value = T)

run_dir <- run_dirs[run_dir_n]
  
setwd(run_dir)

sim_dirs <- grep("Sim_", list.files(), value = T)

for(sim_dir in sim_dirs){
  
  setwd(sim_dir)
  
  param_df <- dir_to_params(run_dir)
  
  # If there are no bc_count dataframes in the chosen sim_dir, need to create an 
  # empty summ_stat dataframe. 
  if(length(grep("bc_counts_tot_mean_R", list.files(), value = T)) == 0){
    
    bc_summ_stats <- data.frame(q = NA, sample = NA, qD = NA, Passage = NA, 
                                type = NA, qD_diss = NA)
    
  } else {
    
    # Read in count and total + mean-R output. 
    bc_counts_Rs <- sapply(grep("bc_counts_tot_mean_R", list.files() , value = T), read.csv, simplify = F, USE.NAMES = T, stringsAsFactors = F)
    # Replace names with passage names
    names(bc_counts_Rs) <- paste0("P", sub(".*R_P([0-9]).csv", replacement = "\\1", names(bc_counts_Rs)))
    # Split into control and drug-treatment groups and drop barcode id columns. 
    CO_dfs <- lapply(bc_counts_Rs, function(x){x[,grep("CO", colnames(x))]})
    DT_dfs <- lapply(bc_counts_Rs, function(x){x[,grep("DT", colnames(x))]})
    # Insert Passage information into colnames. 
    for(i in seq_along(CO_dfs)){
      colnames(CO_dfs[[i]]) <- sub('(?<=.{3})', paste0("_", names(CO_dfs)[[i]]), colnames(CO_dfs[[i]]), perl=TRUE)
    }
    for(i in seq_along(DT_dfs)){
      colnames(DT_dfs[[i]]) <- sub('(?<=.{3})', paste0("_", names(DT_dfs)[[i]]), colnames(DT_dfs[[i]]), perl=TRUE)
    }
    # Remove barcode counts not present in any of the chosen group.
    CO_dfs <- lapply(CO_dfs, function(x){x[rowSums(x == 0) < ncol(x) ,]})
    DT_dfs <- lapply(DT_dfs, function(x){x[rowSums(x == 0) < ncol(x) ,]})
    # Counts
    bc_CO_counts <- lapply(CO_dfs, function(x){x[, grep("_mR|_R", colnames(x), invert = T)]})
    bc_DT_counts <- lapply(DT_dfs, function(x){x[, grep("_mR|_R", colnames(x), invert = T)]})
    
    
    ###############################################################################################
    
    # Pass the counts dataframes through the summary stats collection function. 
    # (found in 'CBC_Analysis_Custom_Functions').
    
    bc_CO_summ_stats <- lapply(bc_CO_counts, collect_summ_stats)
    bc_CO_summ_stats <- bind_rows(bc_CO_summ_stats)
    bc_DT_summ_stats <- lapply(bc_DT_counts, collect_summ_stats)
    bc_DT_summ_stats <- bind_rows(bc_DT_summ_stats)
    bc_summ_stats <- bind_rows(bc_CO_summ_stats, bc_DT_summ_stats)
  
    ###############################################################################################
    
  }
  

  rep_params <- param_df[rep(seq_len(nrow(param_df)), each = nrow(bc_summ_stats)) ,]
  bc_summ_stats <- bind_cols(bc_summ_stats, rep_params)
  master_df <- bind_rows(master_df, bc_summ_stats)

  cat(ts_count)
  ts_count <- ts_count + 1
  
  setwd("../")
  
}
  
setwd("../")
  

# Save the master_df per output directory. 
csv_name <- paste0("CBC_", out_id, "_run-", run_dir_n, 
                      "_Exp_Output_Summary_Stats.csv")

# Save the master_df.
write.csv(master_df, csv_name, row.names = F)


###############################################################################################
