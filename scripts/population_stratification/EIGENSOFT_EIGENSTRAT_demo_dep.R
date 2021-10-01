## Function to perform normalization and regression on a genotype matrix,
## and to replace the original columns by regression residuals, epsilon = Y - X*beta.
## Here k is the number of columns to use (see eq 12 of PPR2006). When k = 0,
## this function can be used to normalize the binary SNP matrix.
regressMatrix <- function(X_data, k) {
  # Normalize according to eq 2 of PPR2006
  mu_matrix <- t(matrix(colMeans(X_data$data), length(colMeans(X_data$data)), nrow(X_data$data)))
  M_matrix <- (X_data$data - mu_matrix) / sqrt(mu_matrix * (1-mu_matrix))
  if (any(is.na(M_matrix))) {
    message("NAs found.")
  }
  if (k == 0) {
    return(M_matrix)
  } else {
    # Create matrix to return
    to_return <- matrix(ncol = ncol(X_data$data), nrow = nrow(X_data$data))
    
    # Fill in columns 
    to_return[,1:k] <- M_matrix[,1:k]
    
    for (j in (k+1):ncol(M_matrix)) {
      # Subset matrix
      X_local <- M_matrix[,(j-k):j]
      # Regress 
      fit <- lm(X_local[,k+1] ~ X_local[,-(k+1)]-1)
      # Update
      to_return[,j] <- fit$residuals
    }
    
    # Return
    return(to_return)
  }
}

## Function to compute significance of largest eigenvalue 
## in the correlation matrix. We closely follow p. 2078 of
## PPR2006.
require(RMTstat)
getTWPValue <- function(M) {
  # 1. Assume M is computed as in Equations 1, 2, and 3.
  # Use regressMatrix function to remove LD if necessary
  m <- nrow(M)
  
  # 2. Compute X = MM'. X is m x m
  X <- M %*% t(M)
  
  # 3. Order the eigenvalues of X so that ...
  eigen_X <- eigen(X, symmetric = TRUE, only.values = TRUE)$values
  eigen_X <- eigen_X[1:(m-1)]
  
  # 4. Estimate n' from Equation 10
  n_prime <- (m + 1) * sum(eigen_X)^2 / (((m-1) * sum(eigen_X^2)) - sum(eigen_X)^2)
  
  # 5. The largest eigenvalue of M is lambda_1. Set l
  lambda_1 <- max(eigen_X)
  l <- (m-1) * lambda_1 / sum(eigen_X)
  
  # 6. Normalize l with Equations 5-7, where the effective number of markers n' replaces n
  mu_mn <- (sqrt(n_prime - 1) + sqrt(m))^2 / n_prime
  sigma_mn <- (sqrt(n_prime - 1) + sqrt(m))/n_prime * (1/sqrt(n_prime - 1) + 1/sqrt(m))^(1/3)
  tw_stat <- (l - mu_mn) / sigma_mn
  
  # Compute and return p-value
  # Using two-tailed test for conservativeness
  p_val <- 2 * min(ptw(tw_stat, beta = 1, lower.tail = T), ptw(tw_stat, beta = 1, lower.tail = F))
  return(p_val)
}

## Set up a function calling msprime to generate 
## synthetic chromosomes
reticulate::conda_install(packages = "msprime")
reticulate::source_python('msprime_model.py')

## Function to generate N individual haplotypes 
## using a backend coalescent simulator. The second
## argument, num_blocks, decides the number of 
## independent chromosomes to simulate.
getExHaplotypes <- function(N, num_blocks) {
  # Build the matrix to be returned
  data_set <- matrix(nrow = N,
                     ncol = 0)

  # Build the column names to be added later
  col_names <- c()
  # Build the block delimiter vector to be returned
  block_bounds <- c(1)

  # Modify the matrix with simulated msprime data
  for (b in 1:num_blocks) {
    new_subarray <- getHaplotypes(N, 1e5)
    if (b < num_blocks) {
      block_bounds  <- c(block_bounds,
                         block_bounds[length(block_bounds)]+dim(new_subarray)[2])
    }
    data_set <- cbind(data_set, new_subarray)
    col_names <- c(col_names, paste0("rs", b, "_", 1:dim(new_subarray)[2]))
  }

  # Include column names for dataset
  colnames(data_set) <- col_names

  # Store output as a list
  to_return <- list(data = data_set,
                    bounds = block_bounds)
  return(to_return)
}

## Function to generate N individual haplotypes 
## using simple coin flips, with heads probability 
## determined by the allele frequency. Allele frequencies
## are arbitrarily drawn from a uniform distribution and fixed
## beforehand. The second argument, num_ind_markers, decides 
## the number of independent alleles.
# getExHaplotypes <- function(N, num_ind_markers) {
#   # Create population allele frequencies
#   allele_freqs <- runif(n = num_ind_markers, min = 0.2, max = 0.8)
#   
#   # Generate array
#   out_array <- do.call("cbind",
#                        lapply(allele_freqs, function(x) {rbinom(N,1,x)}))
#   colnames(out_array) <- paste0("rs", 1:num_ind_markers)
#   
#   # Return
#   to_return <- list(data = out_array,
#                     bounds = 1:num_ind_markers)
#   return(to_return)
# }

## Load some required packages
library(genio)
library(AssocTests)
library(doParallel)
registerDoParallel()
library(dplyr)
library(flintyR)
library(ggplot2)

## Set simulation parameters
N_SIMULATIONS <- 200 
N_CHROM <- 10
N_IND <- 100

## Create dataframe for storing simulated results
results_df <- data.frame(seed = numeric(),
                         EIGENSTRATp_a1 = numeric(),
                         EIGENSTRATp_a5 = numeric(),
                         EIGENSTRATmom = numeric(),
                         EIGENSOFT_a1 = numeric(),
                         EIGENSOFT_a5 = numeric(),
                         flinty = numeric())

## Perform experiment N_SIMULATIONS (= 200) times
sink("EIGENSOFT_dep_log.Rout", split = TRUE)
for (i in 1:N_SIMULATIONS) {
  # Set seed
  set.seed(i * 10)
  
  # Generate data
  X_data <- getExHaplotypes(N = N_IND, num_blocks = N_CHROM)
  #X_data <- getExHaplotypes(N = N_IND, num_ind_markers = N_CHROM)
  
  # Perform LD pruning 
  message(date(), ": Pruning genotype matrix on simulated seed ", i*10, "...")
  genio::write_plink("X_data", t(X_data$data))
  system("./plink --bfile X_data --indep-pairwise 50 10 0.1")
  system("./plink --bfile X_data --extract plink.prune.in --make-bed --out pruned_X_data")
  
  # Perform regression on genotype matrix
  X_norm_reg0 <- regressMatrix(X_data, k = 5)
  
  # Run EIGENSTRAT with LD pruning (suggested by PPR2006)
  message(date(), ": Running EIGENSTRAT with LD pruning on simulated seed ", i*10, "...")
  pruned_X_data <- read_plink("pruned_X_data")
  write.table(pruned_X_data$X, file = "eigenstratG.eg.txt", quote = FALSE,
              sep = "", row.names = FALSE, col.names = FALSE)
  eigenstrat_result <- eigenstrat(genoFile = "eigenstratG.eg.txt", outFile.Robj = "eigenstrat.result.list",
                                  outFile.txt = "eigenstrat.result.txt", rm.marker.index = NULL,
                                  rm.subject.index = NULL, miss.val = 9, num.splits = 10,
                                  topK = NULL, signt.eigen.level = 0.05, signal.outlier = FALSE,
                                  iter.outlier = 5, sigma.thresh = 6)
  EIGENSTRATp_K_a005 <- eigenstrat_result$topK
  eigenstrat_result <- eigenstrat(genoFile = "eigenstratG.eg.txt", outFile.Robj = "eigenstrat.result.list",
                                  outFile.txt = "eigenstrat.result.txt", rm.marker.index = NULL,
                                  rm.subject.index = NULL, miss.val = 9, num.splits = 10,
                                  topK = NULL, signt.eigen.level = 0.01, signal.outlier = FALSE,
                                  iter.outlier = 5, sigma.thresh = 6)
  EIGENSTRATp_K_a001 <- eigenstrat_result$topK
  
  # Run self-implemented conservative algorithm with regression (suggested by PPR2006)
  message(date(), ": Running EIGENSTRAT with regression (K=1,2,3,4,5) on simulated seed ", i*10, "...")
  EIGENSTRAT_p_value <- getTWPValue(X_norm_reg0)
  
  # Run Patterson et al.'s EIGENSOFT (smartpca submodule)
  message(date(), ": Running EIGENSOFT/smartpca on simulated seed ", i*10, "...")
  system("./EIG-7.2.1/bin/smartpca -p X_data.pca.par > EIGENSOFT_out") # run EIGENSOFT
  system("foo=$(grep -n 'Tracy-Widom' EIGENSOFT_out | grep : | cut -f 1 -d:); head -n $(($foo + 11)) EIGENSOFT_out | tail -11 > twstats.txt") # read twstats
  twstats_df <- data.table::fread("twstats.txt") %>% as.data.frame()
  EIGENSOFT_K_a005 <- sum(twstats_df$`p-value` < 0.05)
  EIGENSOFT_K_a001 <- sum(twstats_df$`p-value` < 0.01)
  
  # Run flintyR
  message(date(), ": Running flinty on simulated seed ", i*10, "...")
  flinty_p_value <- getPValue(X = X_data$data, block_boundaries = X_data$bounds)
  
  # Save results 
  message(date(), ": Saving results for simulated seed ", i*10, "...")
  results_df <- rbind(results_df,
                      data.frame(seed = i * 10,
                                 EIGENSTRATp_a1 = EIGENSTRATp_K_a001,
                                 EIGENSTRATp_a5 = EIGENSTRATp_K_a005,
                                 EIGENSTRATmom = EIGENSTRAT_p_value,
                                 EIGENSOFT_a1 = EIGENSOFT_K_a001,
                                 EIGENSOFT_a5 = EIGENSOFT_K_a005,
                                 flinty = flinty_p_value))
  print(data.frame(seed = i * 10,
                   EIGENSTRATp_a1 = EIGENSTRATp_K_a001,
                   EIGENSTRATp_a5 = EIGENSTRATp_K_a005,
                   EIGENSTRATmom = EIGENSTRAT_p_value,
                   EIGENSOFT_a1 = EIGENSOFT_K_a001,
                   EIGENSOFT_a5 = EIGENSOFT_K_a005,
                   flinty = flinty_p_value))

  # Clear garbage
  gc()
}

sink()

write.csv(results_df, 
          file = "EIGENSOFT_results_dep_df.csv")

## Compute empirical FPR
EIGENSTRATp_a1_fpp <- mean(sapply(results_df$EIGENSTRATp_a1, function(x) {ifelse(x > 1, 1, 0)}))
EIGENSTRATp_a5_fpp <- mean(sapply(results_df$EIGENSTRATp_a5, function(x) {ifelse(x > 1, 1, 0)}))
EIGENSTRATmom_fpp <- mean(sapply(results_df$EIGENSTRATmom, function(x) {ifelse(x < 0.05, 1, 0)}))
EIGENSOFT_a1_fpp <- mean(sapply(results_df$EIGENSOFT_a1, function(x) {ifelse(x > 1, 1, 0)}))
EIGENSOFT_a5_fpp <- mean(sapply(results_df$EIGENSOFT_a5, function(x) {ifelse(x > 1, 1, 0)}))
flinty_fpp <- mean(sapply(results_df$flinty, function(x) {ifelse(x < 0.05, 1, 0)}))

comparison <- c(EIGENSTRATp_a1_fpp, EIGENSTRATp_a5_fpp,
                EIGENSTRATmom_fpp, EIGENSOFT_a1_fpp,
                EIGENSOFT_a5_fpp, flinty_fpp)
names(comparison) <- c("EIGENSTRATp_a1", "EIGENSTRATp_a5", 
                       "EIGENSTRATmom", "EIGENSOFT_a1", 
                       "EIGENSOFT_a5", "flinty")

## Print FPRs for comparison
print(comparison)


################################
##### Analysis of Results ######
################################

#EIGENSOFT_df <-read.csv("vignettes/EIGENSOFT_results_df.csv")
#hist(EIGENSOFT_df$EIGENSTRATmom, main = "Histogram of p-values (Independent Markers)", xlab = "p-value")
#sum(EIGENSOFT_df$EIGENSOFT_a5> 0) / nrow(EIGENSOFT_df)
#sum(EIGENSOFT_df$EIGENSTRATp_a5 > 0) / nrow(EIGENSOFT_df)
#EIGENSOFT_df %>% subset(EIGENSOFT_a5 >0)

# # Summarize results for Alan's implementation of the MOM estimator
# sum(EIGENSOFT_df$EIGENSTRATmom < 0.05) / nrow(EIGENSOFT_df) # FPR estimate 
# hist(EIGENSOFT_df$EIGENSTRATmom, main = "Histogram of p-values (Independent Markers)", xlab = "p-value")
# qqplot(qunif(ppoints(200)), EIGENSOFT_df$EIGENSTRATmom, col = "blue", pch = 1, cex = 0.7, 
#        xlab = "Theoretical Quantile", ylab = "Empirical Quantile", main = "QQ Plot of p-values Under Null") 
# abline(a = 0, b = 1, lty = "dashed")
# 
# # Summarize results for running EIGENSOFT's smartpca
# visual_df <- data.frame(FPR = c(sum(EIGENSOFT_df$EIGENSOFT_a1> 0) / nrow(EIGENSOFT_df), sum(EIGENSOFT_df$EIGENSOFT_a5> 0) / nrow(EIGENSOFT_df)),
#                         THRESHOLD = c(0.01,0.05))
# visual_df <- visual_df %>% mutate(UPPER = sapply(FPR, function(x) {x + 1.96 * sqrt(x * (1-x) / 200)}),
#                                   LOWER = sapply(FPR, function(x) {max(0, x - 1.96 * sqrt(x * (1-x) / 200))}))
# ggplot(visual_df, aes(x = THRESHOLD, y = FPR)) +
#   geom_point(colour = "blue") +
#   geom_errorbar(aes(ymin = LOWER, ymax = UPPER), size = 0.5, colour = "blue") +
#   geom_abline(slope = 1, intercept = 0, lty = "dashed") +
#   xlim(c(-0.01,0.106)) + ylim(c(-0.01,0.106)) +
#   theme_bw() +
#   xlab(expression(paste("Significance Level, ", alpha))) + 
#   ylab("Type I Error Rate") +
#   ggtitle("smartpca FPR Estimates") +
#   theme(plot.title = element_text(face="bold", hjust = 0.5))
# 
# #ggsave(filename = "smartpca_FPR.jpg", 
# #       width = 4, height = 4.2, 
# #       dpi = 400)
# 
# ALPHA_VEC <- seq(from=0,to=1,length.out=1e3)
# FPR_VEC <- c()
# for (alpha in ALPHA_VEC) {
#   FPR_VEC <- c(FPR_VEC, sum(EIGENSOFT_df$EIGENSTRATmom < alpha) / length(EIGENSOFT_df$EIGENSTRATmom))
# }
# 
# fpr_df <- data.frame(ALPHA = ALPHA_VEC,
#                      FPR = FPR_VEC)
# 
# fpr_df <- fpr_df %>% mutate(UPPER = sapply(FPR, function(x) {x + 1.96 * sqrt(x * (1-x) / 200)}),
#                                   LOWER = sapply(FPR, function(x) {max(0, x - 1.96 * sqrt(x * (1-x) / 200))}))
# 
# ggplot(fpr_df, aes(x = ALPHA, y = FPR)) + 
#   geom_line(colour = "blue") +
#   geom_errorbar(aes(ymin = LOWER, ymax = UPPER), size = 0.5, colour = "#ece7f2", alpha = 0.3) +
#   geom_abline(slope = 1, intercept = 0, lty = "dashed") +
#   theme_bw() +
#   xlab(expression(paste("Significance Level, ", alpha))) +
#   ylab("Type I Error Rate") +
#   ggtitle("R Implementation FPR Estimates") +
#   theme(plot.title = element_text(face="bold", hjust = 0.5))
# 
# ggsave(filename = "R_implement_FPR.jpg", 
#               width = 4, height = 4.2, 
#               dpi = 400)
