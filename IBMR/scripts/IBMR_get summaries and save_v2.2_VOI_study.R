################################################################
### Delta smelt compile IBMR runs, get summaries, save results #
### William Smith (USFWS; BDFWO); 21 June 2022 #################
################################################################
## Modified by Brian Mahardja 3/25/2025
input_path <- here::here("data/data_raw/demo_inputs")
action_path <- here::here("data/data_processed/")
output_path <- here::here("output/model_outputs/outputs_VOI/")
# output_path <- here::here("output/model_outputs/outputs_AdjHist/")

### Load functions and data ###

# Read all files and filter to rds files
fp_abund <- dir(here(output_path), full.names = TRUE)
fp_abund2 <- data.frame(fp_abund) %>%
  filter(grepl("outz_", fp_abund))

fp_abund_char <- as.character(fp_abund2$fp_abund)
extracted_fp_abund_char <- gsub(".*outz_(.*?)\\.rds.*", "\\1", fp_abund_char)

# List of rds output
ls_abund <- map(fp_abund_char, readRDS)

# Name model run
names(ls_abund) <- extracted_fp_abund_char

# Automate lamAB calculation
calculate_statistics <- function(outz) {
  # Initialize a list to store results
  results <- vector("list", length(outz))
  
  # Iterate over each matrix in the list
  for (i in seq_along(outz)) {
    lamAB <- matrix(NA, 19, 11)
    
    for (t in 1:19) { 
      ratio <- outz[[i]][t + 1, 1, ] / outz[[i]][t, 1, ]
      
      lamAB[t, 1] <- mean(ratio, na.rm = TRUE)
      lamAB[t, 2] <- min(ratio, na.rm = TRUE)
      lamAB[t, 3] <- max(ratio, na.rm = TRUE)
      lamAB[t, 4] <- quantile(ratio, 0.025, na.rm = TRUE)
      lamAB[t, 5] <- quantile(ratio, 0.05, na.rm = TRUE)
      lamAB[t, 6] <- quantile(ratio, 0.25, na.rm = TRUE)
      lamAB[t, 7] <- quantile(ratio, 0.5, na.rm = TRUE)
      lamAB[t, 8] <- quantile(ratio, 0.75, na.rm = TRUE)
      lamAB[t, 9] <- quantile(ratio, 0.90, na.rm = TRUE)
      lamAB[t, 10] <- quantile(ratio, 0.975, na.rm = TRUE)
      lamAB[t, 11] <- (sd(ratio, na.rm = TRUE) / mean(ratio, na.rm = TRUE))
    }
    
    colnames(lamAB) <- c("mean", "min", "max", "2.5%", "5%", "25%", "50%", "75%", "90%", "97.5%", "CV")
    results[[i]] <- lamAB
  }
  
  return(results)
}

# Run the function
list_lamAB <- calculate_statistics(ls_abund)
names(list_lamAB) <- extracted_fp_abund_char

# Add water year
wtr.yr<-c(1,1,1,1,1,1,4,4,2,3,3,1,4,5,4,3,1,3,4,5,5) # Sacto WY type wet = 1, critical = 5

# Calculate lam.mn
calculate_geometric_means <- function(results) {
  lam.mn.list <- vector("list", length(results))  # Initialize a list to store results for each matrix
  
  # Iterate over each matrix in the results
  for (i in seq_along(results)) {
    lamAB <- results[[i]]  # Get the current matrix
    
    lam.mn <- numeric(9)  # Initialize a numeric vector of length 9
    
    lam.mn[1] <- exp(mean(log(lamAB[, 7]), na.rm = TRUE))  # Geometric mean 1995-2014
    lam.mn[2] <- exp(mean(log(lamAB[12:19, 7]), na.rm = TRUE))  # Geometric mean 2007-2014
    lam.mn[3] <- exp(mean(log(lamAB[10:19, 7]), na.rm = TRUE))  # Geometric mean 2005-2014
    lam.mn[4] <- exp(mean(log(lamAB[1:11, 7]), na.rm = TRUE))  # Geometric mean 1995-2006
    lam.mn[5] <- exp(mean(log(lamAB[which(wtr.yr[1:19] <= 2), 7]), na.rm = TRUE))  # Geometric mean AN and wet years
    lam.mn[6] <- exp(mean(log(lamAB[which(wtr.yr[1:19] >= 4), 7]), na.rm = TRUE))  # Geometric mean dry and critical years
    lam.mn[7] <- exp(mean(log(lamAB[3:19, 7]), na.rm = TRUE))  # Geometric mean 1997-2014
    lam.mn[8] <- exp(quantile(log(lamAB[3:19, 7]), 0.025, na.rm = TRUE))  # 95% CI lower
    lam.mn[9] <- exp(quantile(log(lamAB[3:19, 7]), 0.975, na.rm = TRUE))  # 95% CI upper
    
    lam.mn.list[[i]] <- lam.mn  # Store the results in the list
  }
  
  return(lam.mn.list)
}

# Run the geom mean function
list_lam.mn <- calculate_geometric_means(list_lamAB)
names(list_lam.mn) <- extracted_fp_abund_char

# Write the tables
write_lamAB_to_txt <- function(listobj, file_prefix = "lamAB") {
  # Iterate over each item in the listobj
  for (i in names(listobj)) {
    # Create a filename based on the list index
    filename <- paste0(file_prefix, "_", i, ".txt")
    
    # Write the vector to a .txt file
    write.table(listobj[[i]], file.path(output_path,filename), row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

write_lamAB_to_txt(list_lamAB)

write_lam_mn_to_txt <- function(lam.mn.list, file_prefix = "lamABmn") {
  # Iterate over each item in the lam.mn.list
  for (i in names(list_lam.mn)) {
    # Create a filename based on the list index
    filename <- paste0(file_prefix, "_", i, ".txt")
    
    # Write the vector to a .txt file
    write.table(lam.mn.list[[i]], file.path(output_path,filename), row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

write_lam_mn_to_txt(list_lam.mn)














