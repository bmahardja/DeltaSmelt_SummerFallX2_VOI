############################################################
### Delta smelt IBMR data summary ##########################
### revised by William Smith (USFWS; BDFWO) 17 Feb 2022 ####
############################################################
### each section reads in a particular type of data, summarizes into year-month-strata means, and 
### ends with a function to return processed data for IBMR
# 1. make.OMR(first), make.WQ(first) - Create tables of physical conditions OMR, Temp, and Secch
# 2. smlt.dist(t,yr) - Calculate DS distribution based on observed distributions
# 3. make.food(t,yr) - Get food variables to make prey density estimates
# 4. make.temp(t,yr) - Create hybrid temp dataset from DSM2 temps and fish survey data
# 5. make.sal(t,yr) - Create hybrid salinity dataset from DSM2 salinities and fish survey data
# 6. egg2larv(TempCv) - Calculate egg to larvae survival based on Bennett 2005
#################################################################
# t = month index(1:12)
# yr = year index(1=1980; 16=1995)
# first = first year index for simulation (usually, first=16)
# TempCv = n.strata-length vector of monthly mean temperature
# spatial strata names: Yolo Sac SDelta EDelta LowSac LowSJ  Conf SSuisunE NSuisunE Marsh SSuisunW NSuisunW
#################################################################
# Load libraries
require(reshape)
require(abind)
require(LaplacesDemon)
require(chron)
require(ks)
library(here)

#################################################################
# 1. Create tables of physical conditions X2, OMR, Temp, and Secchi
# make.OMR function returns an n.year x n.month matrix of mean OMR
# make.WQ function returns n.year x n.month x n.strata matrices of mean observed Temp and Secchi 
### X2 ###
input_path <- here::here("data/data_raw/demo_inputs")

X2 <- read.table(file.path(input_path, 'X2.txt'),header=T)

make.X2 <- function(t,yr) {
 if(t == 1) { X2.out <- X2[(yr-14),(t+1)] } else { # January is last month of each cohort, so use calendar year+1 in Jan
  X2.out <- X2[(yr-15),(t+1)] }
 return(X2.out)
 }

### OMR ###
# Load daily OMR file
OMR_daily<-read.table(file.path(input_path, 'OMR_daily.txt'),header=T)
OMR_table <- with(OMR_daily, tapply(OMR_daily$OMR.final, list(OMR_daily$Year,OMR_daily$Month), FUN=mean, na.rm='T'))

make.OMR <- function(first) {
 OMR <- matrix(NA,length(yearz),12)
 for (y in first:length(yearz)) {
  OMR[y,1] <- OMR_table[paste0(yearz[y]+1),1] # January is last month of each cohort, so use calendar year+1 in Jan
  OMR[y,(2:12)] <- OMR_table[paste0(yearz[y]),2:12]
  }
 return(OMR)
 }

### Secchi depth and Temp, from survey data ###
# load Secchi and Temp file
# table temp, Secchi, and salinity to get year-month-strata means
WQ_dat=read.csv(file=file.path(input_path, 'delta_water_quality_data_with_strata_v2.csv'),header=T,stringsAsFactors = F)
WQ_dat$Year <- sapply(1:nrow(WQ_dat), function(i) ifelse(WQ_dat$Month[i]==12, WQ_dat$Year[i]-1, WQ_dat$Year[i]))  # Corrects December year assignment
Secchi_table <- with(WQ_dat, tapply(Secchi_pos, list(Year,Month,Strata_code), FUN=mean, na.rm=T))
Temp_table <- with(WQ_dat, tapply(Temperature, list(Year,Month,Strata_code), FUN=mean, na.rm=T))
Temp_table.sd <- with(WQ_dat, tapply(Temperature, list(Year,Month,Strata_code), FUN=sd, na.rm=T))
Sal_table <- with(WQ_dat, tapply(Salinity, list(Year,Month,Strata_code), FUN=mean, na.rm=T))

# fit glm to predict missing Secchi values for strata 1,4, and 9, as function of remaining data in other strata
Secchi.mod.dat <- cbind(as.numeric(Secchi_table[,,1]),as.numeric(Secchi_table[,,2]),as.numeric(Secchi_table[,,3]),as.numeric(Secchi_table[,,4]),
 as.numeric(Secchi_table[,,5]),as.numeric(Secchi_table[,,6]),as.numeric(Secchi_table[,,7]),as.numeric(Secchi_table[,,8]),
 as.numeric(Secchi_table[,,9]),as.numeric(Secchi_table[,,10]),as.numeric(Secchi_table[,,11]),as.numeric(Secchi_table[,,12]))

Secchi.mod1.dat <- Secchi.mod.dat
Secchi.mod1.dat[568,] <- NA # remove outliers and points with unusually high leverage
Secchi.mod1 <- glm(Secchi.mod1.dat[,1]~Secchi.mod1.dat[,2]+Secchi.mod1.dat[,3]+Secchi.mod1.dat[,5]+Secchi.mod1.dat[,6])

Secchi.mod4.dat <- Secchi.mod.dat
Secchi.mod4.dat[670,] <- NA # remove outliers and points with unusually high leverage
Secchi.mod4 <- glm(Secchi.mod4.dat[,4]~Secchi.mod4.dat[,2]+Secchi.mod4.dat[,3]+Secchi.mod4.dat[,6]+Secchi.mod4.dat[,11])

Secchi.mod9.dat <- Secchi.mod.dat
Secchi.mod9.dat[604,] <- NA # remove outliers and points with unusually high leverage 
Secchi.mod9 <- glm(Secchi.mod9.dat[,9]~Secchi.mod9.dat[,2]+Secchi.mod9.dat[,7]+Secchi.mod9.dat[,8]+Secchi.mod9.dat[,12])
 
# replace missing Secchi values with glm predictions
for (i in 1:dim(Secchi_table)[1]) {
 for (j in 1:dim(Secchi_table)[2]) {
  if(is.na(Secchi_table[i,j,1])) { Secchi_table[i,j,1]  <- (Secchi.mod1$coefficients[1]+Secchi.mod1$coefficients[2]*Secchi_table[i,j,2]+
   Secchi.mod1$coefficients[3]*Secchi_table[i,j,3]+Secchi.mod1$coefficients[4]*Secchi_table[i,j,5]+Secchi.mod1$coefficients[5]*Secchi_table[i,j,6])
   }
  if(is.na(Secchi_table[i,j,4])) { Secchi_table[i,j,4]  <- (Secchi.mod4$coefficients[1]+Secchi.mod4$coefficients[2]*Secchi_table[i,j,2]+
   Secchi.mod4$coefficients[3]*Secchi_table[i,j,3]+Secchi.mod4$coefficients[4]*Secchi_table[i,j,6]+Secchi.mod4$coefficients[5]*Secchi_table[i,j,11])
   }
  if(is.na(Secchi_table[i,j,9])) { Secchi_table[i,j,9]  <- (Secchi.mod9$coefficients[1]+Secchi.mod9$coefficients[2]*Secchi_table[i,j,2]+
   Secchi.mod9$coefficients[3]*Secchi_table[i,j,7]+Secchi.mod9$coefficients[4]*Secchi_table[i,j,8]+Secchi.mod9$coefficients[5]*Secchi_table[i,j,12])
   }}}

# fit glm to predict missing Sal values for strata 1,4, and 9, as function of remaining data in other strata
Sal.mod.dat <- cbind(as.numeric(Sal_table[,,1]),as.numeric(Sal_table[,,2]),as.numeric(Sal_table[,,3]),as.numeric(Sal_table[,,4]),
 as.numeric(Sal_table[,,5]),as.numeric(Sal_table[,,6]),as.numeric(Sal_table[,,7]),as.numeric(Sal_table[,,8]),
 as.numeric(Sal_table[,,9]),as.numeric(Sal_table[,,10]),as.numeric(Sal_table[,,11]),as.numeric(Sal_table[,,12]))

Sal.mod1.dat <- Sal.mod.dat
Sal.mod1.dat[597,] <- NA # remove outliers and points with unusually high leverage
Sal.mod1 <- glm(Sal.mod1.dat[,1]~Sal.mod1.dat[,c(3,7,10,11)],family = Gamma(link = "log"))

Sal.mod4.dat <- Sal.mod.dat
Sal.mod4.dat[c(419,652),] <- NA # remove outliers and points with unusually high leverage
Sal.mod4 <- glm(Sal.mod1.dat[,4]~Sal.mod4.dat[,c(2,3,5,6,7)],family = Gamma(link = "log"))

Sal.mod9.dat <- Sal.mod.dat
Sal.mod9.dat[c(32:33,412),] <- NA # remove outliers and points with unusually high leverage 
Sal.mod9 <- glm(Sal.mod9.dat[,9]~Sal.mod9.dat[,c(3,6,10,11,12)],family = Gamma(link = "log"))
 
# replace missing Sal values with glm predictions
for (i in 1:dim(Sal_table)[1]) {
 for (j in 1:dim(Sal_table)[2]) {
  if(is.na(Sal_table[i,j,1])) { Sal_table[i,j,1]  <- exp(Sal.mod1$coefficients[1]+Sal.mod1$coefficients[2]*Sal_table[i,j,3]+
   Sal.mod1$coefficients[3]*Sal_table[i,j,7]+Sal.mod1$coefficients[4]*Sal_table[i,j,10]+Sal.mod1$coefficients[5]*Sal_table[i,j,11])
   }
  if(is.na(Sal_table[i,j,4])) { Sal_table[i,j,4]  <- exp(Sal.mod4$coefficients[1]+Sal.mod4$coefficients[2]*Sal_table[i,j,2]+
   Sal.mod4$coefficients[3]*Sal_table[i,j,3]+Sal.mod4$coefficients[4]*Sal_table[i,j,5]+Sal.mod4$coefficients[5]*Sal_table[i,j,6]+Sal.mod4$coefficients[6]*Sal_table[i,j,7])
   }
  if(is.na(Sal_table[i,j,9])) { Sal_table[i,j,9]  <- exp(Sal.mod9$coefficients[1]+Sal.mod9$coefficients[2]*Sal_table[i,j,3]+
   Sal.mod9$coefficients[3]*Sal_table[i,j,6]+Sal.mod9$coefficients[4]*Sal_table[i,j,10]+Sal.mod9$coefficients[5]*Sal_table[i,j,11]+Sal.mod9$coefficients[6]*Sal_table[i,j,12])
   }}}

make.WQ <- function(first) {
 WQ <- array(NA,dim=c(length(yearz),12,n.strata,3))
 for (y in first:length(yearz)) {
  WQ[y,1,(1:n.strata),1] <- Secchi_table[paste0(yearz[y]+1),1,(1:n.strata)] # January is last month of each cohort, so use calendar year+1 in Jan
  WQ[y,1,(1:n.strata),2] <- Temp_table[paste0(yearz[y]+1),1,(1:n.strata)]
  WQ[y,1,(1:n.strata),3] <- Sal_table[paste0(yearz[y]+1),1,(1:n.strata)]
  WQ[y,(2:12),(1:n.strata),1] <- Secchi_table[paste0(yearz[y]),2:12,(1:n.strata)]
  WQ[y,(2:12),(1:n.strata),2] <- Temp_table[paste0(yearz[y]),2:12,(1:n.strata)]
  WQ[y,(2:12),(1:n.strata),3] <- Sal_table[paste0(yearz[y]),2:12,(1:n.strata)]
  }
 return(WQ)
 }

#################################################################
# 2. Calculate DS distribution based on observed distributions
# smlt.dist function returns a vector with proportion smelt in each spatial strata

# load the files containing the observed DS distributions
DS.dist=read.table(file.path(input_path, 'DS_distribution_JanDec_1995_2015.txt'),header=T)
PEL=read.table(file.path(input_path, 'PEL_DSLCM.txt'),header=T)

# convert observed densities to occupancy probabilities
vol<-c(100510941,144940444,225381539,65449147,89444941,259500691,163419100,153797952,76796487,121672916,107178813,184941122) # mean DSM2 volumes
cov.mns<-c(16.233575,2.36091,54.575090,4.127917,74.23726)
cov.sds<-c(4.591553,3.706014,29.654643,10.041023,11.00871)
Year<-DS.dist[,1]
Month<-DS.dist[,2]
PEL.Year<-PEL[,1]
temp1.distr<-temp2.distr<-array(NA,dim=c(12,12,length(yearz)))
distr<-array(NA,dim=c(12,12,length(yearz)))
for (i in first:length(yearz)) {
for (j in c(1:12)) {
	ifelse(Month[i]==1,temp1.distr[(1:12),j,i] <- as.matrix(DS.dist[which(Year==1979+i+1 & Month==j),(3:14)]), # January is the 12th month of the life cycle, so it gets the next year's value
		temp1.distr[(1:12),j,i] <- as.matrix(DS.dist[which(Year==1979+i & Month==j),(3:14)]))
		}}
for (i in first:length(yearz)) {
for (j in c(1:12)) {
	for (k in c(1:12)) {
		temp1.distr[k,j,i] <- vol[k]*temp1.distr[k,j,i] # volumetric expansion to abundance estimate
		}}
for (j in c(1:12)) {
	for (k in c(1:12)) {
		temp2.distr[k,j,i] <- max(1e-6,temp1.distr[k,j,i]/sum(temp1.distr[(1:n.strata),j,i])) # estimate occupancy probability (non-zero) 
		}}
	distr[3,1,i] <- PEL[which(PEL.Year==1979+i),'Jan']
	distr[3,2,i] <- PEL[which(PEL.Year==1979+i-1),'Feb']
	distr[3,3,i] <- PEL[which(PEL.Year==1979+i-1),'Mar1']
	distr[3,4,i] <- PEL[which(PEL.Year==1979+i),'Apr']
	distr[3,5,i] <- PEL[which(PEL.Year==1979+i),'May']
	distr[3,6,i] <- PEL[which(PEL.Year==1979+i),'Jun']
	distr[3,12,i] <- PEL[which(PEL.Year==1979+i),'Dec']
for (j in c(1:6,12)) {
	for (k in c(1:2,4:12)) {
		distr[k,j,i] <- (1-distr[3,j,i])*(temp2.distr[k,j,i]/sum(temp2.distr[c(1:2,4:12),j,i]))
		}}
for (j in 7:11) {
	distr[,j,i] <- temp2.distr[,j,i]
	}
	}

smlt.dist<-function(t,yr){

prbs<-distr[(1:n.strata),t,yr]

names(prbs)<-c("Yolo", "Sac", "SDelta", "EDelta", "LowSac", "LowSJ",  "Conf", "SSuisunE", 
               "NSuisunE", "Marsh", "SSuisunW", "NSuisunW")

return(prbs)
}
 
#######################################################################
# 3. Get food variables to make prey density estimates
# make food function returns an array with prey density, 
# sd(log(prey density), and mean(proportion 0s) in each spatial strata

# Load the files containing the observed zoop distributions and densities
zoopdat1=read.table(file.path(input_path, 'Zoop.limno_____.txt'),header=T)
zoopdat2=read.table(file.path(input_path, 'Zoop.othcaljuv_.txt'),header=T)
zoopdat3=read.table(file.path(input_path, 'Zoop.pdiapjuv__.txt'),header=T)
zoopdat4=read.table(file.path(input_path, 'Zoop.othcalad__.txt'),header=T)
zoopdat5=read.table(file.path(input_path, 'Zoop.acartela__.txt'),header=T)
zoopdat6=read.table(file.path(input_path, 'Zoop.othclad___.txt'),header=T)
zoopdat7=read.table(file.path(input_path, 'Zoop.allcopnaup.txt'),header=T)
zoopdat8=read.table(file.path(input_path, 'Zoop.daphnia___.txt'),header=T)
zoopdat9=read.table(file.path(input_path, 'Zoop.othcyc____.txt'),header=T)
zoopdat10=read.table(file.path(input_path, 'Zoop.other_____.txt'),header=T)
zoopdat11=read.table(file.path(input_path, 'Zoop.eurytem___.txt'),header=T)
zoopdat12=read.table(file.path(input_path, 'Zoop.pdiapfor__.txt'),header=T)
zoopdatall<-abind(zoopdat1[,3:14],zoopdat2[,3:14],zoopdat3[,3:14],zoopdat4[,3:14],zoopdat5[,3:14],zoopdat6[,3:14],zoopdat7[,3:14],zoopdat8[,3:14],zoopdat9[,3:14],zoopdat10[,3:14],zoopdat10[,3:14],zoopdat12[,3:14],along=3)
zoopdat<-array(NA,dim=c(nrow(zoopdatall),ncol(zoopdatall),n.prey))  # Switch last column (Yolo) to first column
zoopdat[,1,]<-zoopdatall[,12,]
zoopdat[,2:12,]<-zoopdatall[,1:11,]

zoopP01=read.table(file.path(input_path, 'Zoop.P0.limno_____.txt'),header=T)
zoopP02=read.table(file.path(input_path, 'Zoop.P0.othcaljuv_.txt'),header=T)
zoopP03=read.table(file.path(input_path, 'Zoop.P0.pdiapjuv__.txt'),header=T)
zoopP04=read.table(file.path(input_path, 'Zoop.P0.othcalad__.txt'),header=T)
zoopP05=read.table(file.path(input_path, 'Zoop.P0.acartela__.txt'),header=T)
zoopP06=read.table(file.path(input_path, 'Zoop.P0.othclad___.txt'),header=T)
zoopP07=read.table(file.path(input_path, 'Zoop.P0.allcopnaup.txt'),header=T)
zoopP08=read.table(file.path(input_path, 'Zoop.P0.daphnia___.txt'),header=T)
zoopP09=read.table(file.path(input_path, 'Zoop.P0.othcyc____.txt'),header=T)
zoopP010=read.table(file.path(input_path, 'Zoop.P0.other_____.txt'),header=T)
zoopP011=read.table(file.path(input_path, 'Zoop.P0.eurytem___.txt'),header=T)
zoopP012=read.table(file.path(input_path, 'Zoop.P0.pdiapfor__.txt'),header=T)
zoopP0all<-abind(zoopP01[,3:14],zoopP02[,3:14],zoopP03[,3:14],zoopP04[,3:14],zoopP05[,3:14],zoopP06[,3:14],zoopP07[,3:14],zoopP08[,3:14],zoopP09[,3:14],zoopP010[,3:14],zoopP010[,3:14],zoopP012[,3:14],along=3)
zoopP0<-array(NA,dim=c(nrow(zoopP0all),ncol(zoopP0all),n.prey))  # Switch last column (Yolo) to first column
zoopP0[,1,]<-zoopP0all[,12,]
zoopP0[,2:12,]<-zoopP0all[,1:11,]

Year<-zoopdat1[,1]
Day<-zoopdat1[,2]
hldr<-array(NA,dim=c(12,n.prey,3)) # strata x prey type x mean, sd, P0

make.food<-function(t,yr) {
dayz.temp <- vector()
for (mn in 1:12) { # get number of days in each month, varies in leap years
 dayz.temp[mn] <- days.in.month[yr,(mn+1)]
 }
for (h in 1:12) {
for (i in 1:n.prey) {
 if (t==1) { hldr[h,i,1] <- mean(exp(as.matrix(zoopdat[which(Year==1979+yr+1 & Day < (dayz.temp[t]+1)),h,i]))) # January observations # January is the 12th month of the life cycle, so it gets the next year's value
  hldr[h,i,2] <- sd(as.matrix(zoopdat[which(Year==1979+yr+1 & Day < (dayz.temp[t]+1)),h,i]))
  hldr[h,i,3] <- mean(as.matrix(zoopP0[which(Year==1979+yr+1 & Day < (dayz.temp[t]+1)),h,i])) }
 if (t>1) { hldr[h,i,1] <- mean(exp(as.matrix(zoopdat[which(Year==1979+yr & Day > sum(dayz.temp[1:(t-1)]) & Day < (sum(dayz.temp[1:t])+1)),h,i])))
  hldr[h,i,2] <- sd(as.matrix(zoopdat[which(Year==1979+yr & Day > sum(dayz.temp[1:(t-1)]) & Day < (sum(dayz.temp[1:t])+1)),h,i])) 
  hldr[h,i,3] <- mean(as.matrix(zoopP0[which(Year==1979+yr & Day > sum(dayz.temp[1:(t-1)]) & Day < (sum(dayz.temp[1:t])+1)),h,i])) }
  }}
  
 return(hldr)
 }

#######################################################################
# 4. Create hybrid temp dataset from DSM2 temps and fish survey data
# make temp function reads in the DSM2 temperature database,
# fits a model of DSM2.temp=f(fish survey temps), and
# returns mean monthly temps in each spatial strata

# Load DSM2 daily temps
tempDSM=read.table(file.path(input_path,'temp.txt'),header=T)

# DSM2 temps for 1995-2010
temp.string.names <- c('Year','jday','Temp','Strata_code')
tempDSM <- cbind(tempDSM[,1:2],tempDSM[,14],tempDSM[,3:13]) # Switch last column (Yolo) to first column
tempDSM.string <- cbind(tempDSM[,1:2],tempDSM[,3],rep(1,nrow(tempDSM))) # format DSM2 daily temps for tapply() function
colnames(tempDSM.string) <- temp.string.names
for (i in 2:n.strata) {
	a1 <- cbind(tempDSM[,1:2],tempDSM[,(i+2)],rep(i,nrow(tempDSM)))
	colnames(a1) <- temp.string.names
	tempDSM.string <- rbind(tempDSM.string,a1)
	}
tempDSM.string <- cbind(tempDSM.string,month.day.year(tempDSM.string[,2],origin=c(month=1,day=1,year=tempDSM.string[,1])))
Temp_table.DSM2.mean <- with(tempDSM.string, tapply(Temp, list(Year,month,Strata_code), FUN=mean, na.rm=T))
Temp_table.DSM2.sd <- with(tempDSM.string, tapply(Temp, list(Year,month,Strata_code), FUN=sd, na.rm=T))

# make Temp_table.DSM2.mean into a matrix for glm
first.yr <- c(2,32,52) # row index for [first full year of DSM2 temp data in DSM2 data, in fish survey data, first year to predict missing DSM2]
last.yr <- c(22,52,60)# row index for [last full year of DSM2 temp data in DSM2 data, in fish survey data, lastst year to predict missing DSM2]
Temp_table.mat.names <- c('Year','Month','Strata_code','Temp')
a2 <- rep(1,nrow(Temp_table.DSM2.mean[(first.yr[1]:last.yr[1]),,]))
a3 <- rep(1,nrow(Temp_table.DSM2.mean[(first.yr[1]:last.yr[1]),,]))
Temp_table.DSM2.mat <- cbind(rownames(Temp_table.DSM2.mean[(first.yr[1]:last.yr[1]),,]),a2,a3,Temp_table.DSM2.mean[(first.yr[1]:last.yr[1]),1,1])
colnames(Temp_table.DSM2.mat) <- Temp_table.mat.names

# i=months, months are grouped based on glm model selection
i1 <- c(1:3,10:12) # winter and fall months
i2 <- c(4:6,rep(NA,3)) # spring months
i3 <- c(7:9,rep(NA,3)) # summer months
season.index<-cbind(i1,i2,i3)
i.elements <- c(6,3,3)
# j=strata, strata are grouped based on glm model selection
j1 <- c(1:2,4:6,8:10,12) # most strata are grouped
j2 <- c(3,7,11,rep(NA,6)) # strata 3, 7, and 11 get unique regression coefficient
strata.index<-cbind(j1,j2)
j.elements <- c(9,3)

for (g in 1:3) {
for (i in season.index[(1:i.elements[g]),g]) {
for (j in strata.index[(1:j.elements[1]),1]) {
	a2 <- rep(g,nrow(Temp_table.DSM2.mean[(first.yr[1]:last.yr[1]),,])) # seasonal index
	a3 <- rep(1,nrow(Temp_table.DSM2.mean[(first.yr[1]:last.yr[1]),,])) # strata grouping index
	a4 <- cbind(rownames(Temp_table.DSM2.mean[(first.yr[1]:last.yr[1]),,]),a2,a3,Temp_table.DSM2.mean[(first.yr[1]:last.yr[1]),i,j]) # append seasonal and strata grouping indices
	colnames(a4) <- Temp_table.mat.names
	Temp_table.DSM2.mat <- rbind(Temp_table.DSM2.mat,a4)
	}
for (j in strata.index[(1:j.elements[2]),2]) { # strata 3, 7, and 11 get unique regression coefficient
	a2 <- rep(g,nrow(Temp_table.DSM2.mean[(first.yr[1]:last.yr[1]),,]))
	a3 <- rep(j,nrow(Temp_table.DSM2.mean[(first.yr[1]:last.yr[1]),,]))
	a4 <- cbind(rownames(Temp_table.DSM2.mean[(first.yr[1]:last.yr[1]),,]),a2,a3,Temp_table.DSM2.mean[(first.yr[1]:last.yr[1]),i,j])
	colnames(a4) <- Temp_table.mat.names
	Temp_table.DSM2.mat <- rbind(Temp_table.DSM2.mat,a4)
	}}}
Temp_table.DSM2.mat<- as.matrix(Temp_table.DSM2.mat)

# make Temp_table into a matrix for glm
a2 <- rep(1,nrow(Temp_table[(first.yr[2]:last.yr[2]),,]))
a3 <- rep(1,nrow(Temp_table[(first.yr[2]:last.yr[2]),,]))
Temp_table.mat <- cbind(rownames(Temp_table[(first.yr[2]:last.yr[2]),,]),a2,a3,Temp_table[(first.yr[2]:last.yr[2]),1,1])
colnames(Temp_table.mat) <- Temp_table.mat.names
for (g in 1:3) {
for (i in season.index[(1:i.elements[g]),g]) {
for (j in strata.index[(1:j.elements[1]),1]) {
	a2 <- rep(g,nrow(Temp_table[(first.yr[2]:last.yr[2]),,]))
	a3 <- rep(j,nrow(Temp_table[(first.yr[2]:last.yr[2]),,]))
	a4 <- cbind(rownames(Temp_table[(first.yr[2]:last.yr[2]),,]),a2,a3,Temp_table[(first.yr[2]:last.yr[2]),i,j])
	colnames(a4) <- Temp_table.mat.names
	Temp_table.mat <- rbind(Temp_table.mat,a4)
	}
for (j in strata.index[(1:j.elements[2]),2]) {
	a2 <- rep(g,nrow(Temp_table[(first.yr[2]:last.yr[2]),,]))
	a3 <- rep(j,nrow(Temp_table[(first.yr[2]:last.yr[2]),,]))
	a4 <- cbind(rownames(Temp_table[(first.yr[2]:last.yr[2]),,]),a2,a3,Temp_table[(first.yr[2]:last.yr[2]),i,j])
	colnames(a4) <- Temp_table.mat.names
	Temp_table.mat <- rbind(Temp_table.mat,a4)
	}}}
Temp_table.mat <- as.matrix(Temp_table.mat)

# model DSM2 temps as a function of fish monitoring temps
mod <- glm(as.numeric(Temp_table.DSM2.mat[,4]) ~ as.factor(Temp_table.DSM2.mat[,2]) + as.factor(Temp_table.DSM2.mat[,3]) + as.numeric(Temp_table.mat[,4]))

# predict DSM2 monthly temps, given fish monitoring temps in years 2011-2014
for (j in 1:1) { # strata index; 1:2,4:6,8:10,12 share intercept, 3,7,11 get unique factor coefficient (=unique intercept)
 Surv2DSM1 <- mod$coefficients[1] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),1:3,j]
 Surv2DSM2 <- mod$coefficients[1] + mod$coefficients[2] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),4:6,j]
 Surv2DSM3 <- mod$coefficients[1] + mod$coefficients[3] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),7:9,j]
 Surv2DSM4 <- mod$coefficients[1] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),10:12,j]
 
 Surv2DSM <- cbind(Surv2DSM1,Surv2DSM2,Surv2DSM3,Surv2DSM4)
 }
for (j in 2:2) {
 Surv2DSM1 <- mod$coefficients[1] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),1:3,j]
 Surv2DSM2 <- mod$coefficients[1] + mod$coefficients[2] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),4:6,j]
 Surv2DSM3 <- mod$coefficients[1] + mod$coefficients[3] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),7:9,j]
 Surv2DSM4 <- mod$coefficients[1] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),10:12,j]
 
 Surv2DSM.temp <- cbind(Surv2DSM1,Surv2DSM2,Surv2DSM3,Surv2DSM4)
 Surv2DSM <- abind(Surv2DSM,Surv2DSM.temp,along=3)
 }
for (j in 3:3) {
 Surv2DSM1 <- mod$coefficients[1] + mod$coefficients[5] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),1:3,j]
 Surv2DSM2 <- mod$coefficients[1] + mod$coefficients[2] + mod$coefficients[5] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),4:6,j]
 Surv2DSM3 <- mod$coefficients[1] + mod$coefficients[3] + mod$coefficients[5] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),7:9,j]
 Surv2DSM4 <- mod$coefficients[1] + mod$coefficients[5] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),10:12,j]
 
 Surv2DSM.temp <- cbind(Surv2DSM1,Surv2DSM2,Surv2DSM3,Surv2DSM4)
 Surv2DSM <- abind(Surv2DSM,Surv2DSM.temp,along=3)
 }
for (j in 4:6) {
 Surv2DSM1 <- mod$coefficients[1] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),1:3,j]
 Surv2DSM2 <- mod$coefficients[1] + mod$coefficients[2] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),4:6,j]
 Surv2DSM3 <- mod$coefficients[1] + mod$coefficients[3] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),7:9,j]
 Surv2DSM4 <- mod$coefficients[1] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),10:12,j]
 
 Surv2DSM.temp <- cbind(Surv2DSM1,Surv2DSM2,Surv2DSM3,Surv2DSM4)
 Surv2DSM <- abind(Surv2DSM,Surv2DSM.temp,along=3)
 }
for (j in 7:7) {
 Surv2DSM1 <- mod$coefficients[1] + mod$coefficients[6] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),1:3,j]
 Surv2DSM2 <- mod$coefficients[1] + mod$coefficients[2] + mod$coefficients[6] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),4:6,j]
 Surv2DSM3 <- mod$coefficients[1] + mod$coefficients[3] + mod$coefficients[6] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),7:9,j]
 Surv2DSM4 <- mod$coefficients[1] + mod$coefficients[6] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),10:12,j]
 
 Surv2DSM.temp <- cbind(Surv2DSM1,Surv2DSM2,Surv2DSM3,Surv2DSM4)
 Surv2DSM <- abind(Surv2DSM,Surv2DSM.temp,along=3)
 }
for (j in 8:10) {
 Surv2DSM1 <- mod$coefficients[1] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),1:3,j]
 Surv2DSM2 <- mod$coefficients[1] + mod$coefficients[2] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),4:6,j]
 Surv2DSM3 <- mod$coefficients[1] + mod$coefficients[3] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),7:9,j]
 Surv2DSM4 <- mod$coefficients[1] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),10:12,j]
 
 Surv2DSM.temp <- cbind(Surv2DSM1,Surv2DSM2,Surv2DSM3,Surv2DSM4)
 Surv2DSM <- abind(Surv2DSM,Surv2DSM.temp,along=3)
 }
for (j in 11:11) {
 Surv2DSM1 <- mod$coefficients[1] + mod$coefficients[4] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),1:3,j]
 Surv2DSM2 <- mod$coefficients[1] + mod$coefficients[2] + mod$coefficients[4] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),4:6,j]
 Surv2DSM3 <- mod$coefficients[1] + mod$coefficients[3] + mod$coefficients[4] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),7:9,j]
 Surv2DSM4 <- mod$coefficients[1] + mod$coefficients[4] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),10:12,j]
 
 Surv2DSM.temp <- cbind(Surv2DSM1,Surv2DSM2,Surv2DSM3,Surv2DSM4)
 Surv2DSM <- abind(Surv2DSM,Surv2DSM.temp,along=3)
 }
for (j in 12:12) {
 Surv2DSM1 <- mod$coefficients[1] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),1:3,j]
 Surv2DSM2 <- mod$coefficients[1] + mod$coefficients[2] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),4:6,j]
 Surv2DSM3 <- mod$coefficients[1] + mod$coefficients[3] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),7:9,j]
 Surv2DSM4 <- mod$coefficients[1] + mod$coefficients[7]*Temp_table[(first.yr[3]:last.yr[3]),10:12,j]
 
 Surv2DSM.temp <- cbind(Surv2DSM1,Surv2DSM2,Surv2DSM3,Surv2DSM4)
 Surv2DSM <- abind(Surv2DSM,Surv2DSM.temp,along=3)
 }

# function to retrieve temp summaries
make.temp<-function(t,yr) {
 temp.hldr <- matrix(NA,n.strata,2)
if (yr<first.yr[2]) {
 for (h in 1:12) {
 if (t == 1) { temp.hldr[h,1] <- Temp_table.DSM2.mean[(yr-9),t,h] # January observations # January is the 12th month of the life cycle, so it gets the next year's value
  temp.hldr[h,2] <- Temp_table.DSM2.sd[(yr-9),t,h] } else { 
  temp.hldr[h,1] <- Temp_table.DSM2.mean[(yr-10),t,h] 
   temp.hldr[h,2] <- Temp_table.DSM2.sd[(yr-10),t,h] }
   }} else {
  for (h in 1:12) {
  if (t == 1) { temp.hldr[h,1] <- Surv2DSM[(yr-30),t,h] # January observations # January is the 12th month of the life cycle, so it gets the next year's value
   temp.hldr[h,2] <- mean(Temp_table.DSM2.sd[,t,h],na.rm=T) } else { # use month-strata mean sd for 2011-2014
  temp.hldr[h,1] <- Surv2DSM[(yr-31),t,h] 
   temp.hldr[h,2] <- mean(Temp_table.DSM2.sd[,t,h],na.rm=T) }
   }}
   
 if (yr==first.yr[2]) {
 if (t < 11) {
 if (t >1) {
 for (h in 1:12) {
 if (t == 1) { temp.hldr[h,1] <- Temp_table.DSM2.mean[(yr-9),t,h] # January observations # January is the 12th month of the life cycle, so it gets the next year's value
  temp.hldr[h,2] <- Temp_table.DSM2.sd[(yr-9),t,h] } else { temp.hldr[h,1] <- Temp_table.DSM2.mean[(yr-10),t,h] 
   temp.hldr[h,2] <- Temp_table.DSM2.sd[(yr-10),t,h] }
   }}}}
	
return(temp.hldr)
	}
	
#################################################################
# 5. Salinity data
# function to retrieve sal summaries
make.sal<-function(t,yr) {
 sal.hldr <- matrix(NA,n.strata,1)
  if (t==1) {
  for (h in 1:12) {
   sal.hldr[h,1] <- Sal_table[(yr+21),t,h]
   }} else {
  for (h in 1:12) {
   sal.hldr[h,1] <- Sal_table[(yr+20),t,h]
   }}
  #for (h in 1:12) {
  # sal.hldr[h,1] <- Sal_table[(yr+20),t,h]
  # }
 return(sal.hldr)
 }

#################################################################
# 6. Calculate egg to larvae survival based on Bennett 2005
# egg2larv function returns an n.strata-length vector of survival probabilities
# TempCv = mean monthly temperature C 

egg2larv<-function(TempCv){
EggDay=40.1-1.5*TempCv # days from spawn to feeding
LrvDays=30-EggDay # days as yolk-sac larvae
Se=invlogit(-2.35+0.45*TempCv-0.016*TempCv^2) # proportion hatching
#Se=invlogit(-0.0731*TempCv^2+2.0536*TempCv-13.02)

Se*exp(-0.035*LrvDays)
}
