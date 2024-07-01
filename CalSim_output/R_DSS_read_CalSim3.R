# Script to convert necessary CalSim3 dss information into csv
# Identify your working directory for saving outputs of interest
root <- "~/GitHub/DeltaSmelt_SummerFallX2_VOI"
setwd(root)

path_hydro <- file.path(root,"CalSim_output")

# The following libraries need to be installed and loaded
# NOTE: You also need to have HEC-DSSVue installed on your computer
# See: https://www.hec.usace.army.mil/software/hec-dssvue/downloads.aspx

library(tidyverse)
library(stringr)
library(lubridate)
library(rJava)


#############
#Read DSS file

# The following function for is used for turning CalSim time stamps into R dates. 

from_time_stamp <- function(x) {
  day_ref <- as.Date("1899-12-30")
  return(day_ref+x/1440)
}



# Run this workaround if your R session crashes when running .jinit() - below
# This issue occurs in base R versions 4.2 and later
# In lieu of this workaround, you can also install a patched R version 
# E.g., https://cran.r-project.org/bin/windows/base/rpatched.html

replacement <- function(category = "LC_ALL") {
  
  if (identical(category, "LC_MESSAGES"))
    return("")
  
  category <- match(category, .LC.categories)
  if (is.na(category)) 
    stop("invalid 'category' argument")
  .Internal(Sys.getlocale(category))
  
}
base <- asNamespace("base")
environment(replacement) <- base
unlockBinding("Sys.getlocale", base)
assign("Sys.getlocale", replacement, envir = base)
lockBinding("Sys.getlocale", base)


# This code establishes your java connection with HEC-DSSVue

# Specify your own location for 'HEC-DSSVue'
dss_location <- "C:\\Program Files\\HEC\\HEC-DSSVue\\" 

# Specify your own location for the 'jar' sub-folder
# This identifies all possible java executables that meet be needed
jars <- c(list.files("C:\\Program Files\\HEC\\HEC-DSSVue\\jar")) 

jars <- paste0(dss_location, "jar/", jars)

# Specify your own location for the 'lib' sub-folder
libs <- "-Djava.library.path=C:\\Program Files\\HEC\\HEC-DSSVue\\lib\\"

.jinit(classpath = jars, parameters = libs)


##########
# Function to assemble the dataset

# Identify the DSS file you want to access with dss_input

dss_data_pull<-function(dss_input="D:\\2023-01-06 - CalSim3 example file for ReROC\\CalSim3_2040MED_120722_DRAFT_wDWRv705update_wCCdraftBC\\DSS\\output\\CS3_L2020_DV_2021_ext_2040MED"){
  # Open the DSS file through rJava
  dssFile <- .jcall("hec/heclib/dss/HecDss", "Lhec/heclib/dss/HecDss;",   method="open", dss_input)
  #Delta Outflow (OUT in Dayflow)
  java.NDOI_MIN <- dssFile$get("/CALSIM/NDOI_MIN/FLOW//1MON/L2020A/") 
  java.NDOI_ADD <- dssFile$get("/CALSIM/NDOI_ADD/FLOW//1MON/L2020A/") 
  OUTFLOW=data.frame(Date=java.NDOI_MIN$times %>% from_time_stamp,OUTFLOW=java.NDOI_MIN$values+java.NDOI_ADD$values)
  #Old and Middle River flow (OMR)
  java.OMR <- dssFile$get("/CALSIM/C_OMR014/CHANNEL//1MON/L2020A/") 
  OMR=data.frame(Date=java.OMR$times %>% from_time_stamp,OMR=java.OMR$values)
  #X2 (previous month)
  java.X2 <- dssFile$get("/CALSIM/X2_PRV/X2-POSITION-PREV//1MON/L2020A/") 
  X2=data.frame(Date=java.X2$times %>% from_time_stamp,X2_prev=java.X2$values)
  
  final_data_frame= OUTFLOW %>% left_join(OMR) %>% left_join(X2) %>%
    mutate(X2_current=lead(X2_prev,n=1))
  return(final_data_frame)
}


#Use the function to create data frame
AltF80_data <- dss_data_pull(dss_input="D:\\2024-07-01 - Delta Smelt SummerFallX2 VOI CalSim Runs\\CS3_DS_VOI_Alt_F80\\DSS\\output\\CS3_LTO_NAA_2022MED_03102024_L2020A_DV_dp")
AltF74_data <- dss_data_pull(dss_input="D:\\2024-07-01 - Delta Smelt SummerFallX2 VOI CalSim Runs\\CS3_DS_VOI_Alt_F74\\DSS\\output\\CS3_DS_VOI_F74_DV_dp_Aug80")
AltS74_data <- dss_data_pull(dss_input="D:\\2024-07-01 - Delta Smelt SummerFallX2 VOI CalSim Runs\\CS3_DS_VOI_Alt_S74\\DSS\\output\\CS3_DS_VOI_S74_DV_dp")
AltS74F80_data <- dss_data_pull(dss_input="D:\\2024-07-01 - Delta Smelt SummerFallX2 VOI CalSim Runs\\CS3_DS_VOI_Alt_S74F80\\DSS\\output\\CS3_DS_VOI_S74F80_DV_dp")
AltNoX2_data <- dss_data_pull(dss_input="D:\\2024-07-01 - Delta Smelt SummerFallX2 VOI CalSim Runs\\CS3_DS_VOI_Alt_NoX2\\DSS\\output\\CS3_DS_VOI_NoX2_DV_dp")


#Export DSS output files for Delta Smelt LCM model input
write.csv(AltF80_data,file.path(path_hydro,"AltF80_CalSim3_data.csv"),row.names=F)
write.csv(AltF74_data,file.path(path_hydro,"AltF74_CalSim3_data.csv"),row.names=F)
write.csv(AltS74_data,file.path(path_hydro,"AltS74_data_CalSim3_data.csv"),row.names=F)
write.csv(AltS74F80_data,file.path(path_hydro,"AltS74F80_data_CalSim3_data.csv"),row.names=F)
write.csv(AltNoX2_data,file.path(path_hydro,"AltNoX2_CalSim3_data.csv"),row.names=F)

##### Calculate flow input for the Delta Smelt IBMR

AltF80_data <- AltF80_data %>% mutate(scenario="AltF80")
AltF74_data <- AltF74_data %>% mutate(scenario="AltF74")
AltS74_data <- AltS74_data %>% mutate(scenario="AltS74")
AltS74F80_data <- AltS74F80_data %>% mutate(scenario="AltS74F80")
AltNoX2_data <- AltNoX2_data %>% mutate(scenario="AltNoX2")

combined_data <- bind_rows(AltF80_data,AltF74_data,AltS74_data,AltS74F80_data,
                           AltNoX2_data)
data_flow_IMBR <- combined_data %>% rename(X2=X2_current) %>% select(-X2_prev) %>% mutate(year=year(Date), month=month(Date))

#Export data
write.csv(data_flow_IMBR,file.path(path_hydro,"IBMR_FlowData_DeltaSmelt_SummerFallX2_VOI_CalendarYear.csv"),row.names=F)

##### Calculate flow input for the Delta Smelt LCME

#Per Will Smith's excel documentation
#-time is indexed by cohort year, with the first month of the year beginning in April

#Summer Delta Outflow covariate
data_DeltaOutflow <- combined_data %>% mutate(Cohort_Year=ifelse(month(Date)<4,year(Date)-1,year(Date)),Month=month(Date)) %>% 
  #Filter June-August per Smith et al. 2021
  dplyr::filter(Month %in% c(6:8)) %>% 
  #Filter year 1995-2016 for Delta Smelt LCM
  dplyr::filter(Cohort_Year %in% c(1995:2016)) %>%
  #sum flow by month
  mutate(SumOutflow_Jun_Aug=case_when(Month==6 ~ OUTFLOW*30,
                                            Month==7 ~OUTFLOW*31,
                                            Month==8 ~OUTFLOW*31)) %>%
  #summarize by year
  group_by(Cohort_Year,scenario) %>% summarise(Outflow_Jun0Aug0=sum(SumOutflow_Jun_Aug))

#Apr-May OMR covariate
data_OMR_Apr_May <- combined_data %>% mutate(Cohort_Year=ifelse(month(Date)<4,year(Date)-1,year(Date)),Month=month(Date)) %>% 
  #Filter Apr-May per Smith et al. 2021
  dplyr::filter(Month %in% c(4:5)) %>% 
  #Filter year 1995-2016 for Delta Smelt LCM
  dplyr::filter(Cohort_Year %in% c(1995:2016)) %>%
  #summarize by year
  group_by(Cohort_Year,scenario) %>% summarise(OMR_AprMar=mean(OMR))

#June OMR covariate
data_OMR_June <- combined_data %>% mutate(Cohort_Year=ifelse(month(Date)<4,year(Date)-1,year(Date)),Month=month(Date)) %>% 
  #Filter June per Smith et al. 2021
  dplyr::filter(Month %in% c(6)) %>% 
  #Filter year 1995-2016 for Delta Smelt LCM
  dplyr::filter(Cohort_Year %in% c(1995:2016)) %>%
  rename(OMR_Jun=OMR) %>% select(Cohort_Year,scenario,OMR_Jun)

#Dec-Jan OMR covariate
data_OMR_Dec_Jan <- combined_data %>% 
  #Use water year instead for this one
  mutate(Cohort_Year=ifelse(month(Date)<4,year(Date)-1,year(Date)),Month=month(Date)) %>% 
  #Filter Dec-Jan per Smith et al. 2021
  dplyr::filter(Month %in% c(12,1)) %>% 
  #Filter year 1995-2016 for Delta Smelt LCM
  dplyr::filter(Cohort_Year %in% c(1995:2016)) %>%
  group_by(Cohort_Year,scenario) %>% summarise(OMR_DecJan=mean(OMR))

#February OMR covariate
data_OMR_Feb <- combined_data %>% mutate(Cohort_Year=ifelse(month(Date)<4,year(Date)-1,year(Date)),Month=month(Date)) %>% 
  #Filter February per Smith et al. 2021
  dplyr::filter(Month %in% c(2)) %>% 
  #Filter year 1995-2016 for Delta Smelt LCM
  dplyr::filter(Cohort_Year %in% c(1995:2016)) %>%
  rename(OMR_Feb=OMR) %>% select(Cohort_Year,scenario,OMR_Feb)

#March OMR covariate
data_OMR_Mar <- combined_data %>% mutate(Cohort_Year=ifelse(month(Date)<4,year(Date)-1,year(Date)),Month=month(Date)) %>% 
  #Filter March per Smith et al. 2021
  dplyr::filter(Month %in% c(3)) %>% 
  #Filter year 1995-2016 for Delta Smelt LCM
  dplyr::filter(Cohort_Year %in% c(1995:2016)) %>%
  rename(OMR_Mar=OMR) %>%  select(Cohort_Year,scenario,OMR_Mar)

##Merge datasets
data_flowinput <- data_DeltaOutflow %>% left_join(data_OMR_Apr_May) %>% left_join(data_OMR_June) %>% left_join(data_OMR_Dec_Jan) %>%
  left_join(data_OMR_Feb) %>%
  left_join(data_OMR_Mar)

#Export data
write.csv(data_flowinput,file.path(path_hydro,"LCME_FlowData_DeltaSmelt_SummerFallX2_VOI_CohortYear.csv"),row.names=F)
