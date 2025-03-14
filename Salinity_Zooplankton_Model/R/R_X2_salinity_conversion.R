#Script to convert X2 to salinity data for IBMR regions
root <- "~/GitHub/DeltaSmelt_SummerFallX2_VOI"
setwd(root)

path_data <- file.path(root,"Salinity_Zooplankton_Model","Data")
path_r <- file.path(root,"Salinity_Zooplankton_Model","R")
path_output <- file.path(root,"Salinity_Zooplankton_Model","Output")
path_hydro <- file.path(root,"CalSim_output")

library(reshape)
library(tidyverse)
library(stringr)
library(lubridate)
library(conflicted)
library(wql)
conflict_prefer("rename", "dplyr")

# Load X2 data from CalSim3
AltF80_data <- read.csv(file.path(path_hydro,"AltF80_CalSim3_data.csv")) 
AltF74_data <- read.csv(file.path(path_hydro,"AltF74_CalSim3_data.csv")) 
AltS74_data <- read.csv(file.path(path_hydro,"AltS74_data_CalSim3_data.csv")) 
AltS74F80_data <- read.csv(file.path(path_hydro,"AltS74F80_data_CalSim3_data.csv")) 
AltNoX2_data <- read.csv(file.path(path_hydro,"AltNoX2_CalSim3_data.csv")) 

# Combine X2 data
x2_data <- AltF80_data %>% select(Date,X2_current) %>% mutate(Scenario="AltF80") %>%
  bind_rows((AltF74_data %>% select(Date,X2_current) %>% mutate(Scenario="AltF74"))) %>%
  bind_rows((AltS74_data %>% select(Date,X2_current) %>% mutate(Scenario="AltS74"))) %>%
  bind_rows((AltS74F80_data %>% select(Date,X2_current) %>% mutate(Scenario="AltS74F80"))) %>%
  bind_rows((AltNoX2_data %>% select(Date,X2_current) %>% mutate(Scenario="AltNoX2"))) 

x2_data <- na.omit(x2_data) %>% rename(X2 = X2_current) %>% mutate(Month=month(Date))

# Create X2 data frame for all relevant regions
x2_data_expanded <- crossing(x2_data, Region=c("NW Suisun","SW Suisun","NE Suisun","SE Suisun","Confluence", "Suisun Marsh"))

# Load final salinity-X2 model
salX2mod <- readRDS(file.path(path_r,"model_sal_X2.Rdata")) 

# Add the proper predictors (change month and region to factors)
x2_data_expanded <- x2_data_expanded %>% mutate(month_f = as.factor(Month),region_f= as.factor(Region))

# Use the CSAMP X2-Salinity model to convert CalSim3 X2 values to salinity
x2_data_expanded$salinity<-predict(salX2mod,x2_data_expanded,type="response")

# Ensure that there will be no negative salinity values and use the minimum value in Sam's conversion table
summary(salX2mod)
x2_data_expanded$salinity <- ifelse(x2_data_expanded$salinity<0.1,0.1,x2_data_expanded$salinity)

# Finalize data format
x2_data_reformat <- x2_data_expanded %>%
  mutate(year=year(Date)) %>% rename(region=Region, month=Month) %>% select(-Date,-month_f,-region_f,-X2) %>%
  spread(Scenario,salinity) 

#####
# Load Montezuma slough salinity data from CalSim
BD_data <- AltF80_data %>% select(Date,BD_EC_current) %>% mutate(Scenario="AltF80") %>%
  bind_rows((AltF74_data %>% select(Date,BD_EC_current) %>% mutate(Scenario="AltF74"))) %>%
  bind_rows((AltS74_data %>% select(Date,BD_EC_current) %>% mutate(Scenario="AltS74"))) %>%
  bind_rows((AltS74F80_data %>% select(Date,BD_EC_current) %>% mutate(Scenario="AltS74F80"))) %>%
  bind_rows((AltNoX2_data %>% select(Date,BD_EC_current) %>% mutate(Scenario="AltNoX2"))) %>%
  # Convert data to salinity units per discretewq package (https://github.com/InteragencyEcologicalProgram/discretewq)
  mutate(BD_salinity=wql::ec2pss(.data$BD_EC_current / 1000, t = 25))

# Create data for Suisun Marsh
SM_data <- BD_data %>% mutate(year=year(Date),month=month(Date),region="Suisun Marsh") %>% rename(scenario=Scenario) %>%
  select(month,region,year,scenario,BD_salinity) %>% spread(scenario,BD_salinity)
#####

# Remove original Suisun Marsh salinity from X2-salinity model
x2_data_reformat_SM_edit <- x2_data_reformat %>% dplyr::filter(region!="Suisun Marsh") %>% bind_rows(SM_data) %>%
  dplyr::filter(!is.na(AltF74))
  
# Rename column names to sal_
colnames(x2_data_reformat_SM_edit)[4:ncol(x2_data_reformat_SM_edit)] <- paste("sal", colnames(x2_data_reformat_SM_edit)[4:ncol(x2_data_reformat_SM_edit)] , sep = "_")


#Export output file for  model input
write.csv(x2_data_reformat_SM_edit,file.path(path_data,"converted_salinity_data.csv"),row.names=F)
