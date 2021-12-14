
# Cancer Barcode Simulation
# Freddie Whiting - 2021

# Plotting the qD statistics extracted using the CBC_Exp_Output_Summ_Stats 
# script. The qD statistics capture the within-replicate diversity (qD) and 
# the between-replicate diversity dissimilarity (qD_diss). 

################################################################################

rm(list = ls())

# Load in command line arguments: 

args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 3) stop("There should be 3 command-line arguments."))

# First argument - output directory
out_dir_path <- args[1]
# Second argument - custom function script. 
cust_fun_path <- args[2]
# Third argument - output-specific identifier prefix (e.g. 'sigs_psi_0.0_')
out_id <- args[3]

################################################################################
