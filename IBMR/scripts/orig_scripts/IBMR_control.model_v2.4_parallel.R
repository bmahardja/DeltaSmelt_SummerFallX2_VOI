############################################################
### Delta smelt IBMR control model #########################
### William Smith (USFWS; BDFWO); 26 June 2022 #############
############################################################
# v2: updated V and K values, based on latest DSIBM documentation from Rose and Kimmerer
# -added M.act.mult and Dist.act.array to account for actions that scale natural mortality (e.g. contaminants) or fish distributions (e.g. Franks Tract)
# -moved parameters to new file 'IBMR_parameters.R'
# -moved crash.max from pop.model to code, reduced to 2.3e5 to catch more explodey sims before R has problems
# -reduced super.ad to 200 to avoid simulating too many larvae
# v2.2: added FranksTract and Contaminant action effects
# v2.3: added PD multipler by prey type
# v2.4: BoRVOI inputs

### Set working directory to source file location
rm(list=ls(all=TRUE))
# setwd('D:/FWS/R code/Rose BEM/Peterson-Smith-Rose BEM, 2020 revisions/IBMR')
setwd(here::here())
input_path <- here::here("data/data_raw/demo_inputs")
action_path <- here::here("data/data_processed/demo_action_inputs")

### Model parameters ###
# 1. Calibration
# 2. Model dimensions
# 3. Load functions and data
# 4. Define management actions
# 5. Run and compile population model

### 1. Calibration ###
#Mmult<- 1.262 #  monthly natural mortality adjustments
Mmult<- 1.2753 # calibrated to LCME PL2 latent states

### 2. Model dimensions ###
# spatial strata names: Yolo Sac SDelta EDelta LowSac LowSJ  Conf SSuisunE NSuisunE Marsh SSuisunW NSuisunW
# zoop names: limno othcaljuv pdiapjuv othcalad acartela othclad allcopnaup daphnia othcyc other eurytem pdiapfor
# data dimensions: c(year, month, spatial strata, prey)
yearz <- c(1980:2014) # Peterson et al. modeled years 1980:1999
first <- 16 # begin simulation in year 1995
super.ad <- 200 # number of pre-spawn adult super-individual
crash.max <- 2.3e+5 # maximum number of larval production before stopping ('exploding' simulated pop.)
sims <- 330 # number of simulations to run

### 3. Load functions and data ###
source("scripts/IBMR_parameters_v2.R")
source("scripts/ds_functions_v9_2.R") # This loads R libraries and functions to summarize, OMR, smelt distribution, food, and egg survival, v3=observed dist v4=occ predicted dist
source("scripts/IBMR_pop.model_v2.4.R")
move.matrix<-read.table(file.path(input_path,'move.matrix.12strata_v3.txt')) # Load movement rules
DOY<-c(15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349) # day of year midpoints for each month
days.in.month<-read.table(file.path(input_path, 'days.in.month.txt'),header=T)
wtr.yr<-c(1,1,1,1,1,1,4,4,2,3,3,1,4,5,4,3,1,3,4,5,5) # Sacto WY type wet = 1, critical = 5

### 4. Define management actions ###
# Set environmental variables to change, given management action #
OMR.action <- F
Temp.action <- F
Food.spp.action <- F # specified as taxa-specific multiplier
Move.manual.action <- F # specified by directly changing distributions 'by hand'
Secchi.action <- F

action <- 'AltF80'

move.act.mn <- matrix(NA,(n.years+1),12)
Dist.act.array <- M.act.mult <- PD.mult <- secchi.act.array <- Temp.act.array <- X2.act.sal.array <- array(NA,dim=c((n.years+1),12,n.strata))

# Load files defining covariates associated with actions
# OMR
OMR.act <- read.csv(file.path(action_path,'IBMR_OMR_SummerFallX2VOI_input.csv'),header=T)
OMR.act <- OMR.act[which(OMR.act==paste0(action)),]
rownames(OMR.act) <- OMR.act[,2]
OMR.act <- OMR.act[,3:14]
#OMR.act[,c(1:8,11:12)] <- NA

# Temp
Temp.act <- read.table(file.path(action_path,'Wetlands_2_hi_temp_temp.txt'),header=F)

for (y in 1:(n.years+1)) {
for (m in 1:12) {
 Temp.act.array[y,m,] <- as.matrix(Temp.act[,(m+(y-1)*12)])
 }}

# Secchi
secchi.act <- read.table(file.path(action_path,'AWC_Yolo_secchi.txt'),header=F)
for (y in 1:(n.years+1)) {
for (m in 1:12) {
 secchi.act.array[y,m,] <- as.matrix(secchi.act[,(m+(y-1)*12)])
 }}

# X2
X2.dist.act <- read.csv(file.path(action_path,'IBMR_X2_SummerFallX2VOI_input.csv'),header=T)
X2.dist.act <- X2.dist.act[which(X2.dist.act==paste0(action)),]
rownames(X2.dist.act) <- X2.dist.act[,2]
X2.target <- X2.dist.act[,3:14]

# Prey
PD.act <- read.csv(file.path(action_path,'zoop_scalar_output_2024-08-01.csv'),header=T)

PD.act.col <- which(colnames(PD.act)==paste0("sal_",action,"_median"))

PD.mult_spp <- array(NA,dim=c((n.years+1),12,n.strata,n.prey))
PD.mult_table <- with(PD.act, tapply(PD.act[,PD.act.col], list(PD.act$year,PD.act$month,PD.act$region,PD.act$IBMR), FUN=mean, na.rm='T'))
PD.mult_spp[(1:n.years),,(7:12),] <- PD.mult_table[,,c("Confluence","SE Suisun","NE Suisun","Suisun Marsh","SW Suisun", "NW Suisun"),
 c("limno","othcaljuv","pdiapjuv","othcalad","acartela","othclad","allcopnaup","daphnia","othcyc","other","eurytem","pdiapfor")]

# Population distribution
#DS.dist.act <- read.table('Action effects/FirstFlush_dist.txt',header=F)

#for (y in 1:(n.years+1)) {
#for (m in 1:12) {
# Dist.act.array[y,m,] <- as.matrix(DS.dist.act[,(m+(y-1)*12)])
# }}

# X2 move action
for (y in 1:n.years) {
for (m in 1:12) {
 Dist.act.array[y,m,] <- (smlt.dist(m,(15+max(which(abs(X2[(1:20),(m+1)]-X2.target[y,m])==min(abs(X2[(1:20),(m+1)]-X2.target[y,m]))))))
							+smlt.dist(m,(15+max(which(abs(X2[(1:20),(m+1)]-X2.target[y,m])==sort(abs(X2[(1:20),(m+1)]-X2.target[y,m]))[2])))))/2
 }}

#summX2 <- T
#fallX2 <- F
#X2.target <- c(70,75,80,80) # Jul,Aug,Sep.Oct

#for (y in 1:(n.years+1)) {
#if (summX2 == T) {
#if(wtr.yr[y] < 3 & (X2[y,"Jul"]-X2.target[1])>2) {
# Dist.act.array[y,7,] <- (smlt.dist(7,(15+which(abs(X2[,"Jul"]-X2.target[1])==min(abs(X2[,"Jul"]-X2.target[1])))))
#							+smlt.dist(7,(15+which(abs(X2[,"Jul"]-X2.target[1])==sort(abs(X2[,"Jul"]-X2.target[1]))[2]))))/2
# }
#if(wtr.yr[y] < 3 & (X2[y,"Aug"]-X2.target[2])>2) {
# Dist.act.array[y,8,] <- (smlt.dist(7,(15+which(abs(X2[,"Aug"]-X2.target[2])==min(abs(X2[,"Aug"]-X2.target[2])))))
#							+smlt.dist(7,(15+which(abs(X2[,"Aug"]-X2.target[2])==sort(abs(X2[,"Aug"]-X2.target[2]))[2]))))/2
# }}
#if (fallX2 == T) {
#if(wtr.yr[y] < 3 & (X2[y,"Sep"]-X2.target[3])>2) {
# Dist.act.array[y,9,] <- (smlt.dist(7,(15+which(abs(X2[,"Sep"]-X2.target[3])==min(abs(X2[,"Sep"]-X2.target[3])))))
#							+smlt.dist(7,(15+which(abs(X2[,"Sep"]-X2.target[3])==sort(abs(X2[,"Sep"]-X2.target[3]))[2]))))/2
# }
#if(wtr.yr[y] < 3 & (X2[y,"Oct"]-X2.target[4])>2) {
# Dist.act.array[y,10,] <- (smlt.dist(7,(15+which(abs(X2[,"Oct"]-X2.target[4])==min(abs(X2[,"Oct"]-X2.target[4])))))
#							+smlt.dist(7,(15+which(abs(X2[,"Oct"]-X2.target[4])==sort(abs(X2[,"Oct"]-X2.target[4]))[2]))))/2
# }}}

# Doing this first to grab the seeds provided to each worker to replicate in the sequential step
### 5. Parallel Processing Setup ###
library(future)
library(future.apply)
library(progressr)

parallel_start <- Sys.time()
# Does take time to set up the overhead, i.e., more time savings the more sims you run

# Configure future package
# Reserve 2 cores for system (or crank to 0 if you are going to leave the computer)
plan(multisession, workers = availableCores() - 1)

# Create wrapper function for single simulation
run_single_simulation <- function(sim_number) {

  # Saving seed state before starting. However, this is totally dependent on
  # how delta_smelt_run() is formulated. The function is super complex so I am
  # not diving into it unless it's needed. Regardless, this is how you should be
  # able to save seeds for reproducibility purposes.
  rng_state <- .Random.seed

  tryCatch({
    result <- delta_smelt_run()
    return(list(
      success = TRUE,
      data = result,
      sim_number = sim_number,
      seed_state = rng_state,
      error = NULL
    ))
  }, error = function(e) {
    return(list(
      success = FALSE,
      data = NULL,
      sim_number = sim_number,
      seed_state = rng_state,
      error = e$message
    ))
  })
}

# Initialize storage for failed simulations
failed_sims <- list()

# Run simulations with progress tracking
# Keep in mind that the progress bar only updates once the simulation is done
# This means that you won't see movement until close to the end since all
# cores are running simulataneously (and assuming at the same speed)
with_progress({
  p <- progressor(steps = sims)

  # Run all simulations in parallel
  results <- future_lapply(1:sims, function(x) {
    result <- run_single_simulation(x)
    p()  # Update progress once done with model
    return(result)
  }, future.seed = TRUE)
})

# Process results and handle errors
successful_runs <- list()
failed_runs <- list()

for (i in seq_along(results)) {
  if (results[[i]]$success) {
    successful_runs[[length(successful_runs) + 1]] <- list(
      data = results[[i]]$data,
      seed.state = results[[i]]$seed_state
    )
  } else {
    failed_runs[[length(failed_runs) + 1]] <- list(
      sim.number = results[[i]]$sim.number,
      error = results[[i]]$error
    )
  }
}

# Report on failed simulations?
if (length(failed_runs) > 0) {
  message("\nFailed simulations:")
  for (fail in failed_runs) {
    message(sprintf("Simulation %d failed with error: %s",
                    fail$sim.number, fail$error))
  }
}

# Combine Results
if (length(successful_runs) > 0) {
  # Combine successful results
  outz1p <- abind(lapply(successful_runs, function(x) x$data[[1]]), along=3)
  outz2p <- abind(lapply(successful_runs, function(x) x$data[[2]]), along=4)

  # Discard burn-in iterations
  outzp <- outz1p[(5:24),,]

  # # Save results if run time is limiting
  # save(outz, outz1, outz2, file = paste0("results_", format(Sys.time(), "%Y%m%d_%H%M"), ".RData"))

  # Report success
  message("\nSuccessfully completed ", length(successful_runs), " out of ", sims, " simulations")
} else {
  stop("No successful simulations completed")
}

parallel_runtime <- difftime(Sys.time(), parallel_start, units = "mins")

# Clean up
# Reset future plan to sequential (default R, no parallelization)
plan(sequential)
gc()

# Default RNGkind() is Mersenne-Twister. Will have to change this to "L'Ecuyer-CMRG"
# This is what's used by futures

RNGkind("L'Ecuyer-CMRG")
### 5. Run and compile population model ###
sequential_start <- Sys.time()
# Loop through sims without parallel processing
# Loading the seed state from first run in parallel
.Random.seed <- results[[1]]$seed_state
outz<-delta_smelt_run()
outz1<-outz[[1]] # nAB, febFL, ET risk
outz2<-outz[[2]] # spatial distributions
for (s in 2:sims){
  .Random.seed <- results[[s]]$seed_state
  new.outz<-delta_smelt_run()
  new.outz1<-new.outz[[1]]
  new.outz2<-new.outz[[2]]
  outz1<-abind(outz1,new.outz1,along=3)
  outz2<-abind(outz2,new.outz2,along=4)
}
# discard burn-in iterations
outz<-outz1[(5:24),,]
sequential_runtime <- difftime(Sys.time(), sequential_start, unit = "mins")

# Reset rng to default
RNGkind("Mersenne-Twister")
