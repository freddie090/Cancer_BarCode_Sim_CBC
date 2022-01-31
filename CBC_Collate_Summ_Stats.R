
# Cancer Barcode Simulation
# Freddie Whiting - 2021

# Collating the Summary Statistics from the CBC_Exp simulations. 

################################################################################

rm(list = ls())

# Load in command line arguments: 

args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 3) stop("There should be 3 command-line arguments."))

# First argument - output directory
out_dir_path <- args[1]
# Second argument - custom function script. 
cust_fun_path <- args[2]
# Third argument - output-specific identifier prefix (e.g. 'sig_psi_0.0_')
out_id <- args[3]

################################################################################

# Load in the custom CBC functions. 
source(cust_fun_path)

# Set directory to results directory. 
setwd(out_dir_path)

master_df <- data.frame()

ts_count <- 1

# Save the output csv name 
csv_name <- paste0("CBC_Exp_Output_Summary_Stats_", out_id, ".csv")

try(if(file.exists(csv_name) == T) stop("The collated output csv already exists."))

# Need to go through each run directory and combine the summary statistic 
# dataframe from each. 

# Extract all simulation runs 

run_csvs <- grep("CBC_.*Summary_Stats.csv", list.files(), value = T)

for(run_csv in run_csvs){

  run_csv <- read.csv(run_csv, stringsAsFactors = F)
  
  master_df <- bind_rows(master_df, run_csv)

  cat(ts_count)
  
  ts_count <- ts_count + 1
  
}

# Write the finished dataframe. 

setwd("../")

write.csv(master_df, csv_name, row.names = F)



