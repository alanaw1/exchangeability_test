################################################################################
########### Application of Exchangeability Test to 1000 Genomes Data ###########
############### Created by Alan Aw on 25th July 2021 ############################
################################################################################

# commandArgs picks up the variables you pass from the command line
# ARG1 = directory to plink input files
# ARG2 = directory to plink
# ARG3 = number of permutations
# Ex: Rscript test_exchangeability.R /Users/alanaw/Documents/research/pgs/030121/pop_files/GBR/1000G_GBR /Users/alanaw/Documents/research/pgs/101220/plink2
args <- commandArgs(trailingOnly = TRUE)
plink_input_dir <- args[1] # Ex: plink_input_dir <- "/Users/alanaw/Documents/research/pgs/030121/pop_files/GBR/1000G_GBR"
plink_dir <- args[2] # Ex: plink_dir <- "/Users/alanaw/Documents/research/pgs/101220/plink2"
resample_n <- args[3] # Ex: 1000

# Load libraries
library(tidyverse)
library(snpStats)
library(flintyR)
library(doParallel)
registerDoParallel() # use available multiple cores to speed up permutation task 

# For each chromosome
message(date(), ": Generating distance matrices for each chromosome...\n")
dist_list <- list()
for (b in 1:22) {
  message(paste("Computing pairwise distances for Chr", b))
  start_time <- Sys.time()
  # Use PLINK to filter genos to that chrom only
  system(paste0(plink_dir, " --bfile ", 
                plink_input_dir, " --chr ", b, 
                " --make-bed --out ", 
                plink_input_dir, "_", b))
  
  # Load plink file into R using snpStats::read.plink
  geno_matrix <- snpStats::read.plink(paste0(plink_input_dir, "_", b))
  
  # Compute pairwise distances
  pairwise_dists <- as.matrix(dist(geno_matrix$genotypes, "manhattan"))
  
  # Save to dist_list
  dist_list[[b]] <- pairwise_dists
  
  # Remove geno_matrix object loaded by read_plink to save memory
  remove(geno_matrix)
  
  # Force clear garbage 
  gc()
  
  # Print time elapsed
  print(Sys.time() - start_time)
}

# Removing all intermediate files generated (ex: chr-specific bed files and their logs)
message(date(), ": Removing all intermediate files generated...\n")
system(paste0("rm -f ", plink_input_dir, "_*"))

# Compute exchangeability test p-value using 1000 permutations
message(date(), ": Computing exchangeability test p-value...\n")
start_time <- Sys.time()
p_value <- distDataPValue(dist_list, nruns = resample_n)
print(paste0("Exchangeability p-value: ", p_value))
message(paste0("Exchangeability p-value: ", p_value))
print(Sys.time() - start_time)
