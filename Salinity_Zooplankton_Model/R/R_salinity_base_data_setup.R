#Create a new salinity base data to allow inclusion of 2015-2022 data not present in the CSAMP dataset
root <- "~/GitHub/DeltaSmelt_SummerFallX2_VOI"
setwd(root)

path_data <- file.path(root,"Salinity_Zooplankton_Model","Data")
path_r <- file.path(root,"Salinity_Zooplankton_Model","R")
path_output <- file.path(root,"Salinity_Zooplankton_Model","Output")

library(reshape)
library(tidyverse)
library(stringr)
library(lubridate)
library(conflicted)
library(deltamapr)
library(sf)
library(AICcmodavg)
#devtools::install_github("InteragencyEcologicalProgram/discretewq")
require(discretewq)
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("summarise", "dplyr")

#Use same datasets as CSAMP Delta Smelt SDM
wq_data <-  wq(Sources = c("EMP", "STN", "FMWT", "EDSM", "DJFMP",
                           "SKT", "20mm", "Suisun", 
                           "Baystudy", "USBR", "USGS_SFBS")) %>%
  filter(!is.na(Latitude)&!is.na(Longitude)) %>%
  st_as_sf(coords=c("Longitude", "Latitude"),crs=st_crs("WGS84"))
#Sam's package goes until CY 2022

#Delta Smelt IBM layer
IBM_layer <- R_DSIBM

#Join IBM layer with wq data
wq_data_IBM <- st_transform(wq_data, crs = st_crs(IBM_layer))
wq_data_IBM<- st_join(wq_data_IBM,IBM_layer)
st_geometry(wq_data_IBM) <- NULL # remove geometry, coerce to data.frame

#Summarize salinity by region, month, and year
wq_data_sum <- wq_data_IBM %>% group_by(Year,Month,SUBREGION) %>% summarise(Salinity=mean(Salinity,na.rm = T)) %>% rename(year=Year,month=Month,region=SUBREGION) %>%
  filter(year %in% c(1994:2022)&!is.na(region)) %>% rename(sal_base=Salinity) %>% ungroup()


# Pull Dayflow data
dayflow_data_post97<-read.csv(file=file.path(path_data,"dayflow-results-1997-2023.csv")) %>%
  mutate(Date=as.Date(Date,"%m/%d/%Y"))

dayflow_data_pre97<-read.csv(file=file.path(path_data,"dayflowCalculations_pre97.csv")) %>%
  mutate(Date=as.Date(Date,"%m/%d/%Y"))

# Create X2 data by year and month
x2_data <- bind_rows(dayflow_data_pre97,dayflow_data_post97) %>% mutate(year=year(Date),month=month(Date)) %>% group_by(year, month) %>%
  summarise(X2=mean(X2)) %>%
  # Remove calendar years prior to 1994
  filter(year>=1994)

# Add X2 data to wq data
wq_data_sum <- wq_data_sum %>% left_join(x2_data) %>%
  # Remove data points with no X2
  filter(!is.na(X2)) %>%
  # Add water year as a variable
  mutate(WY=ifelse(month>9,year+1,year))

#### Construct models
# Convert month, region, and water year to factors
wq_data_sum <- wq_data_sum %>% mutate(month_f = as.factor(month),region_f = as.factor(region),wy_f = as.factor(WY))

# Z-standardize X2
X2_mean = mean(wq_data_sum$X2)
X2_sd = sd(wq_data_sum$X2)
wq_data_sum <- wq_data_sum %>% mutate(X2_z= (X2-X2_mean)/X2_sd) %>%
  rename(sal=sal_base) %>% mutate(X2_z_quad = X2_z^2)

# Model selection, using raw X2 seems to work just fine, so we will stick with that
fit0 <- glm(sal ~ 1,data=wq_data_sum,family=Gamma(link="log"))
fit1 <- glm(sal ~ X2,data=wq_data_sum,family=Gamma(link="log"))
fit2 <- glm(sal ~ X2 + region_f,data=wq_data_sum,family=Gamma(link="log"))
fit3 <- glm(sal ~ X2 + month_f + region_f,data=wq_data_sum,family=Gamma(link="log"))
fit4 <- glm(sal ~ X2 + month_f + region_f + X2*region_f,data=wq_data_sum,family=Gamma(link="log"))
fit5 <- glm(sal ~ X2 + month_f + region_f + X2*region_f + month_f*region_f,data=wq_data_sum,family=Gamma(link="log"))
fit6 <- glm(sal ~ X2 + I(X2^2) + month_f + region_f + region_f:(X2 + I(X2^2)) + month_f*region_f, data=wq_data_sum ,family=Gamma(link="log"))

lrModels <- list(fit0, fit1, fit2, fit3, fit4, fit5, fit6)
lrNames <- c("Null","fit1","fit2","fit3","fit4","fit5","fit6")
aicWt <- aictab(cand.set=lrModels, modnames=lrNames, sort=TRUE, c.hat=1, second.ord=TRUE)
aicWt <- aicWt %>% mutate(model_description = case_when(Modnames == "Null" ~ "1",
                                                        Modnames == "fit1" ~ "X2",
                                                        Modnames == "fit2" ~ "X2 + region_f",
                                                        Modnames == "fit3" ~ "X2 + month_f + region_f",
                                                        Modnames == "fit4" ~ "X2 + month_f + region_f + X2*region_f",
                                                        Modnames == "fit5" ~ "X2 + month_f + region_f + X2*region_f + month_f*region_f",
                                                        Modnames == "fit6" ~ "X2 + I(X2^2) + month_f + region_f + region_f:(X2 + I(X2^2)) + month_f*region_f"))

aicWt

# Export AICc data for table
write.csv(aicWt,file.path(path_output,"X2_Salinity_model_AICc.csv"),row.names=T)

# Model 6 is the best 
summary(fit6)
# McFadden's R squared
1-(fit6$deviance/fit6$null.deviance)
with(summary(fit6), 1 - deviance/null.deviance)

# Create prediction data
pred_data <- wq_data_sum 
pred_data$pred_m6 = predict(fit6,newdata=pred_data,type="response")

# Create plot
p.sal6 <- ggplot(pred_data, aes(x=X2,y=sal)) + 
  geom_point(size=0.4, alpha=0.5) +
  geom_line(data=pred_data, aes(x = X2, y = pred_m6, group=month_f), size=0.7, color="red", alpha=0.5) +
  labs(x= "X2",y="Salinity (PSU)") + 
  scale_x_continuous(breaks=c(seq(44,94,10))) + facet_wrap(~region_f) +
  theme_classic() + theme(panel.background = element_rect(color="black"),
                          axis.text.y=element_text(size=11),axis.text.x=element_text(size=11),axis.title=element_text(size=12,face="bold"),
                          legend.position=c(0.1,0.9),legend.title=element_blank(), plot.margin=unit(c(0.3,0.6,0.4,0.4),"cm"))
p.sal6

# Export plot
tiff(filename=file.path(path_output,"Figure_X2_salinity_model.tiff"),
     type="cairo",
     units="in", 
     width=8, #10*1, 
     height=6, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
p.sal6
dev.off()

# Export data for potential use later on
write.csv(wq_data_sum,file.path(path_data,"base_salinity_data_for_X2salmodel.csv"),row.names=T)

# Save best model for use
saveRDS(fit6,file.path(path_r,"model_sal_X2.Rdata"))
