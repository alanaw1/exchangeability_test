###################################################################
########### Conversion of WVS Data to Analysable Format ###########
############### Created by Alan Aw on 5th June 2021 ###############
###################################################################

# commandArgs picks up the variables you pass from the command line
args <- commandArgs(trailingOnly = TRUE)

# Load libraries
library(tidyverse)

# Source script containing R functions
source("helper_functions.R")

# Load the RDS file downloaded from WVS online portal
cat(">>> Reading in WVS Data...\n")
#wvs <- readRDS("F00007762-WV6_Data_R_v20180912.Rds") # [!] Change to argument
wvs <- readRDS(args[1])

# Read Michael's file containing included variables for cultural distances
cat(">>> Reading in MM Codebook CSV file...\n")
#MM_metadata <- read.csv("MM_2020/included-variables.csv") # [!] Change to argument
MM_metadata <- read.csv(args[2])

# Perform Filtering and Recoding 
wv_num_filtered_ <- keepCulDistFeatures(wvs_obj = wvs, culdist_metadata = MM_metadata)

cat(">>> Saving Data as CSV file...\n")
write.csv(wv_num_filtered_, file = "all_MM_num_filtered_recode.csv")
cat(">>> Saving Data as RData file...\n")
save(wv_num_filtered_, file = "all_MM_num_filtered_recode.RData")