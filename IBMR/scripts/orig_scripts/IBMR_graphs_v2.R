############################################################
### Delta smelt IBMR graphs ################################
### revised by William Smith (USFWS; BDFWO) 29 Jun 2022 ####
############################################################

# v2: added July and January observed abundance vs. IBMR sim
# -switched FWS.abundances from observed to LCME latent states

### 1. Individual growth (length) plots ###
### 2. Abundance and population growth plots ###
### 3. Enrainment risk plots ###
### 4. Spatial distribution plots ###

yr.seq <- seq((min(yearz)+first-1),max(yearz),by=1)

### 1. Individual growth (length) plots ###
### Plot mean length after 12 months of growth ###
# mean observed February length in MWT 1992-2001, SKT 2002-2018
obs.Febln<-c(63.55,64.60,62.55,65.31,64.24,67.72,61.56,64.37,72.06,62.87,68.44,66.57,62.29,65.94,65.13,66.67,67.64,63.29,69.30,68.22,66.53,68.67,68.88)
obs.Febln.salvage<-c(66.71,68.5,66,74,69.02,72.05,68.26,71.08,70.56,68.32,68.5,71,65.5,68.7,66,71.5,78,69.87,72.42,NA,64)
obs.Marln<-obs.Febln+(76.1-obs.Febln)*(1-exp(-2.98*28/365))

par(family='serif')
plot(yr.seq,apply(outz[(1:nrow(outz)),5,],1,median,na.rm=T),col="red",ylim=c(55,80),type='l',ylab='',xlab='')
lines(yr.seq,apply(outz[(1:nrow(outz)),5,],1,quantile, probs=0.025,na.rm=T),col="red",lty="dotted")
lines(yr.seq,apply(outz[(1:nrow(outz)),5,],1,quantile, probs=0.975,na.rm=T),col="red",lty="dotted")
lines(yr.seq[1:23],obs.Febln[1:23])
lines(yr.seq[1:19],obs.Febln.salvage[1:19],lty="dotted")
points(yr.seq[21],obs.Febln.salvage[21])
legend("bottomleft",legend=c("Observed, surveys","Observed, salvage","Median simulated","95% interval of simulations"),lty=c("solid","dotted","solid","dashed"),col=c("black","black","red","red"))
mtext('Fork length (mm)', side=2,line=2.3,cex=1.13)
mtext('Year', side=1,line=2.3,cex=1.13)
#title('Mean length of adults in February',line=-1.5)

### 2. Abundance and population growth plots ###
### Histogram of change in abundance among simulations ###
hist(outz[20,3,]/outz[5,3,],xlab="Proportional decline",main=" ",ylim=c(0,20))
mtext('Frequency', side=2,line=2.3,cex=1.13)
mtext('Growth', side=1,line=2,cex=1)
title('1999-2014 population growth',line=-1.5)

### Abundance plots ###
outz.sim<-outz

#FWS.abundance<-read.table('FWS.abundance.txt',header=F) # FWS abundance for 1995-2015
FWS.abundance<-read.table(file.path(input_path, 'FWS.abundance_LCME.txt'),header=T)
FWS.abundance<-cbind(FWS.abundance[,2],FWS.abundance[,3],FWS.abundance[,4],FWS.abundance[,5])
TNS.bias <- 1 #0.19
MWT.bias <- 1 #0.42
AB.start <- 1 #first-15
AB.end <- 20 #length(yearz)-15
OMR.start <- c((first-10),(first-9))
OMR.end <- c((length(yearz)-10),(length(yearz)-9))

par(mfrow=c(2,2))
par(family='serif',mar = c(2,3,2,3), oma = c(3,3, 0.5, 1),cex.main=1.5,cex.axis=1.5)

# June
plot(yr.seq[1:20],FWS.abundance[AB.start:AB.end,1],type="l",ylim=c(0.15*min(FWS.abundance[AB.start:AB.end,1]),1.6*max(FWS.abundance[AB.start:AB.end,1])),xlim=c(1994,2019),ylab='ln(TMM)',xlab='Year')
#for (i in 1:dim(outz)[3]) {
# lines(yr.seq,outz.sim[(1:nrow(outz)),1,i]/super[1],col="gray",lty="dotted")
# }
lines(yr.seq,apply(outz[(1:nrow(outz)),1,],1,median,na.rm=T)/super[1],col="red")
lines(yr.seq,apply(outz[(1:nrow(outz)),1,],1,quantile, probs=0.025,na.rm=T)/super[1],col="gray",lty="dotted")
lines(yr.seq,apply(outz[(1:nrow(outz)),1,],1,quantile, probs=0.975,na.rm=T)/super[1],col="gray",lty="dotted")
legend("right",legend=c("LCME estimate","Median simulated","95% interval"),lty=c("solid","solid","dashed"),col=c("black","red","gray"))
mtext("Abundance", side=2,line=1,cex=1.75,outer=T)
title('June',line=-1.5)

# August
plot(yr.seq[1:20],(FWS.abundance[AB.start:AB.end,2]/MWT.bias),type="l",ylim=c(0.85*min(FWS.abundance[AB.start:AB.end,2]/MWT.bias),1.6*max(FWS.abundance[AB.start:AB.end,2]/MWT.bias)),xlim=c(1994,2019),ylab='ln(TNS)',xlab='Year')
#for (i in 1:dim(outz)[3]) {
# lines(yr.seq,outz.sim[(1:nrow(outz)),2,i]/super[1],col="gray",lty="dotted")
# }
lines(yr.seq,apply(outz[(1:nrow(outz)),2,],1,median,na.rm=T)/super[1],col="red")
lines(yr.seq,apply(outz[(1:nrow(outz)),2,],1,quantile, probs=0.025,na.rm=T)/super[1],col="gray",lty="dotted")
lines(yr.seq,apply(outz[(1:nrow(outz)),2,],1,quantile, probs=0.975,na.rm=T)/super[1],col="gray",lty="dotted")
title('August',line=-1.5)

#November
plot(yr.seq[1:20],(FWS.abundance[AB.start:AB.end,3]/TNS.bias),type="l",ylim=c(0.25*min(FWS.abundance[AB.start:AB.end,3]/TNS.bias),1.25*max(FWS.abundance[AB.start:AB.end,3]/TNS.bias)),xlim=c(1994,2019),ylab='ln(FMWT)',xlab='Year')
#for (i in 1:dim(outz)[3]) {
# lines(yr.seq,outz.sim[(1:nrow(outz)),3,i]/super[1],col="gray",lty="dotted")
# }
lines(yr.seq,apply(outz[(1:nrow(outz)),3,],1,median,na.rm=T)/super[1],col="red")
lines(yr.seq,apply(outz[(1:nrow(outz)),3,],1,quantile, probs=0.025,na.rm=T)/super[1],col="gray",lty="dotted")
lines(yr.seq,apply(outz[(1:nrow(outz)),3,],1,quantile, probs=0.975,na.rm=T)/super[1],col="gray",lty="dotted")
title('November',line=-1.5)

#January
plot(yr.seq[1:20],c((FWS.abundance[AB.start:(AB.start+5),4]/MWT.bias),(FWS.abundance[(AB.start+6):AB.end,4])),type="l",ylim=c(0.25*min(c((FWS.abundance[AB.start:(AB.start+5),4]/MWT.bias),(FWS.abundance[(AB.start+6):AB.end,4]))),1.75*max(c((FWS.abundance[AB.start:(AB.start+5),4]/MWT.bias),(FWS.abundance[(AB.start+6):AB.end,4])))),xlim=c(1994,2019),ylab='ln(FMWT)',xlab='Year')
#for (i in 1:dim(outz)[3]) {
# lines(yr.seq,outz.sim[(1:nrow(outz)),4,i]/super[1],col="gray",lty="dotted")
# }
lines(yr.seq,apply(outz[(1:nrow(outz)),4,],1,median,na.rm=T)/super[1],col="red")
lines(yr.seq,apply(outz[(1:nrow(outz)),4,],1,quantile, probs=0.025,na.rm=T)/super[1],col="gray",lty="dotted")
lines(yr.seq,apply(outz[(1:nrow(outz)),4,],1,quantile, probs=0.975,na.rm=T)/super[1],col="gray",lty="dotted")
title('January',line=-1.5)

mtext('Year', outer=T, side=1,line=1,cex=1.75)

### Pop. growth plots ###
#outz1<-readRDS('D:/FWS/R code/Rose BEM/Peterson-Smith-Rose BEM, 2020 revisions/outz_1_19_2021.rds')
FWS.abundance<-read.table(file.path(input_path,'FWS.abundance_LCME.txt'),header=F) # FWS abundance for 1995-2015
FWS.abundance[,2]<-FWS.abundance[,2]/TNS.bias
FWS.abundance[,3]<-FWS.abundance[,3]/MWT.bias
FWS.abundance[(1:6),3]<-FWS.abundance[(1:6),4]/MWT.bias
TNS.bias <- 0.19
MWT.bias <- 0.42

lamAB<-vector()
lamLCM<-vector()
for (t in 1:19) {
 lamAB[t] <- median(outz[t+1,3,]/outz[t,3,])
 lamLCM[t] <- FWS.abundance[t+1,4]/FWS.abundance[t,4]
 }

par(mfrow=c(2,2))
plot(seq(1996,2014,by=1),lamAB,typ='l',col='red',ylim=c(0.18,3.05))
lines(seq(1996,2014,by=1),lamLCM)
plot(lamLCM,lamAB,ylab='Simulated',xlab='Observed')
hist(lamAB,xlim=c(0,2))

# Lambda vs. Sacto WY type
par(mfrow=c(2,1))
par(family='serif',mar = c(2,2,1,1), oma = c(3,3, 0,0),cex.main=1.25,cex.axis=1.25)

hist(lamAB[which(wtr.yr<=2)],xlim=c(0,2.5),ylim=c(0,7),xlab='',main='')
title("Wet and above normal",line=-2)
hist(lamAB[which(wtr.yr>=4)],xlim=c(0,2.5),ylim=c(0,7),xlab='',main='')
title("Dry and critical",line=-2)

mtext('Population growth rate', outer=T, side=1,line=1,cex=1.5)

### 3. Enrainment risk plots ###
OMR <- make.OMR(first)
par(mfrow=c(2,4))
par(family='serif',mar = c(2,3,2,3), oma = c(3,3, 0.5, 3),cex.main=2,cex.axis=1.5)

# April
plot(yr.seq,outz[(1:nrow(outz)),10,1],ylim=c(1*min(outz[(1:nrow(outz)),6:11,],na.rm=T),1*max(outz[(1:nrow(outz)),6:11,],na.rm=T)),col="gray",lty="dotted",type='l',ylab=' ',xlab=' ')
for (i in 2:dim(outz)[3]) {
 lines(yr.seq,outz[(1:nrow(outz)),10,i],col="gray",lty="dotted")
 }
lines(yr.seq,apply(outz[(1:nrow(outz)),10,],1,median,na.rm=T),col="red")
legend(2004,0.35,legend=c("Old and Middle River flow","Median entrainment risk","Simulated entrainment risk"),lty=c("solid","solid","dashed"),col=c("black","red","gray"),cex=1.2)
par(new = T)
plot(yr.seq,OMR[(first:length(yearz)),4],yaxt='n',xaxt='n',typ="l",ylab=' ',xlab=' ',ylim=c(-9500,25000))
abline(-3500,0,lty='dashed')
axis(side = 4)
title("April",line=-2)
# May
plot(yr.seq,outz[(1:nrow(outz)),11,1],ylim=c(1*min(outz[(1:nrow(outz)),6:11,],na.rm=T),1*max(outz[(1:nrow(outz)),6:11,],na.rm=T)),col="gray",lty="dotted",type='l',ylab=' ',xlab=' ')
for (i in 2:dim(outz)[3]) {
 lines(yr.seq,outz[(1:nrow(outz)),11,i],col="gray",lty="dotted")
 }
lines(yr.seq,apply(outz[(1:nrow(outz)),11,],1,median,na.rm=T),col="red")
par(new = T)
plot(yr.seq,OMR[(first:length(yearz)),5],yaxt='n',xaxt='n',typ="l",ylab=' ',xlab=' ',ylim=c(-9500,25000))
abline(-3500,0,lty='dashed')
axis(side = 4)
title("May",line=-2)
# June
plot(yr.seq,outz[(1:nrow(outz)),12,1],ylim=c(1*min(outz[(1:nrow(outz)),6:11,],na.rm=T),1*max(outz[(1:nrow(outz)),6:11,],na.rm=T)),col="gray",lty="dotted",type='l',ylab=' ',xlab=' ')
for (i in 2:dim(outz)[3]) {
 lines(yr.seq,outz[(1:nrow(outz)),12,i],col="gray",lty="dotted")
 }
lines(yr.seq,apply(outz[(1:nrow(outz)),12,],1,median,na.rm=T),col="red")
par(new = T)
plot(yr.seq,OMR[(first:length(yearz)),6],yaxt='n',xaxt='n',typ="l",ylab=' ',xlab=' ',ylim=c(-9500,25000))
abline(-3500,0,lty='dashed')
axis(side = 4)
title("June",line=-2)

frame()

plot(yr.seq,outz[(1:nrow(outz)),6,1],ylim=c(1*min(outz[(1:nrow(outz)),6:11,],na.rm=T),1*max(outz[(1:nrow(outz)),6:11,],na.rm=T)),col="gray",lty="dotted",type='l',ylab=' ',xlab=' ')
for (i in 2:dim(outz)[3]) {
 lines(yr.seq,outz[(1:nrow(outz)),6,i],col="gray",lty="dotted")
 }
lines(yr.seq,apply(outz[(1:nrow(outz)),6,],1,median,na.rm=T),col="red")
par(new = T)
plot(yr.seq,OMR[(first:length(yearz)),12],yaxt='n',xaxt='n',typ="l",ylab=' ',xlab=' ',ylim=c(-9500,25000))
abline(-5000,0,lty='dashed')
axis(side = 4)
title("December",line=-2)

plot(yr.seq,outz[(1:nrow(outz)),7,1],ylim=c(1*min(outz[(1:nrow(outz)),6:11,],na.rm=T),1*max(outz[(1:nrow(outz)),6:11,],na.rm=T)),col="gray",lty="dotted",type='l',ylab=' ',xlab=' ')
for (i in 2:dim(outz)[3]) {
 lines(yr.seq,outz[(1:nrow(outz)),7,i],col="gray",lty="dotted")
 }
lines(yr.seq,apply(outz[(1:nrow(outz)),7,],1,median,na.rm=T),col="red")
par(new = T)
plot(yr.seq,OMR[(first:length(yearz)),1],yaxt='n',xaxt='n',typ="l",ylab=' ',xlab=' ',ylim=c(-9500,25000))
abline(-5000,0,lty='dashed')
axis(side = 4)
title("January",line=-2)

plot(yr.seq,outz[(1:nrow(outz)),8,1],ylim=c(1*min(outz[(1:nrow(outz)),6:11,],na.rm=T),1*max(outz[(1:nrow(outz)),6:11,],na.rm=T)),col="gray",lty="dotted",type='l',ylab=' ',xlab=' ')
for (i in 2:dim(outz)[3]) {
 lines(yr.seq,outz[(1:nrow(outz)),8,i],col="gray",lty="dotted")
 }
lines(yr.seq,apply(outz[(1:nrow(outz)),8,],1,median,na.rm=T),col="red")
par(new = T)
plot(yr.seq,OMR[(first:length(yearz)),2],yaxt='n',xaxt='n',typ="l",ylab=' ',xlab=' ',ylim=c(-9500,25000))
abline(-5000,0,lty='dashed')
axis(side = 4)
title("February",line=-2)

plot(seq((min(yearz)+first-1),max(yearz),,by=1),outz[(1:nrow(outz)),9,1],ylim=c(1*min(outz[(1:nrow(outz)),6:11,],na.rm=T),1*max(outz[(1:nrow(outz)),6:11,],na.rm=T)),col="gray",lty="dotted",type='l',ylab=' ',xlab=' ')
for (i in 2:dim(outz)[3]) {
 lines(seq((min(yearz)+first-1),max(yearz),,by=1),outz[(1:nrow(outz)),9,i],col="gray",lty="dotted")
 }
lines(yr.seq,apply(outz[(1:nrow(outz)),9,],1,median,na.rm=T),col="red")
par(new = T)
plot(yr.seq,OMR[(first:length(yearz)),3],yaxt='n',xaxt='n',typ="l",ylab=' ',xlab=' ',ylim=c(-9500,25000))
abline(-5000,0,lty='dashed')
axis(side = 4)
title("March",line=-2)

mtext("Mean entrainment risk", outer=T,side=2,line=1,cex=1.3)
mtext("Old and Middle River flow (cfs)", outer=T,side=4,line=1,cex=1.3)
mtext('Cohort year', outer=T, side=1,line=1,cex=1.3)

### 4. Spatial distribution plots ###
### Use sp.dist to graph spatial distributions ###
ex.yrs<-c(1995:1997) # choose 3
ex.mns<-c(7:12) # choose 4
n.yrs <- 1+(length(yearz)-first)

str.names <- c('Yolo','Sac','SDel','EDel','LSac','LSJR','Conf','SESu','NESu','Marsh','SWSu','NWSu')
# DS.dist=read.table('Data/DS_distribution_JanDec_1995_2015.txt',header=T)
PEL=read.table('PEL_DSLCM.txt',header=T)

# convert observed densities to occupancy probabilities
vol<-c(100510941,144940444,225381539,65449147,89444941,259500691,163419100,153797952,76796487,121672916,107178813,184941122) # mean DSM2 volumes
Year<-DS.dist[,1]
Month<-DS.dist[,2]
PEL.Year<-PEL[,1]
temp.distr<-array(NA,dim=c(12,12,length(yearz)))
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

smlt.dist2<-function(t,yr){
if(t==1) { prbs<-distr[,1,yr] } # January observations
if(t==2) { prbs<-distr[,2,yr] } # February observations
if(t==3) { prbs<-as.vector(distr[,3,yr]) } # March observations
if(t==4) { prbs<-distr[,4,yr] }  # April observations
if(t==5) { prbs<-distr[,5,yr] }  # May observations
if(t==6) { prbs<-distr[,6,yr] }  # June observations
if(t==7) { prbs<-distr[,7,yr] } # July observations
if(t==8) { prbs<-distr[,8,yr] }  # August observations
if(t==9) { prbs<-distr[,9,yr] } # September observations
if(t==10) { prbs<-distr[,10,yr] } # October observations
if(t==11) { prbs<-distr[,11,yr] } # November observations
if(t==12) { prbs<-distr[,12,yr] }  # December observations

names(prbs)<-c("Yolo", "Sac", "SDelta", "EDelta", "LowSac", "LowSJ",  "Conf", "SSuisunE",
               "NSuisunE", "Marsh", "SSuisunW", "NSuisunW")

return(prbs)
}

# Make distribution data
Dist<-array(NA,dim=c(n.strata,12,n.yrs))

for (yr in first:length(yearz)) { # year loop
for (t in 1:12) {
 Dist[,t,(yr-15)]<-smlt.dist2(t,yr)
 }}


c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

par(mfrow=c(3,6))
par(family='serif',mar = c(1,1,0,1), oma = c(4, 4, 0.5, 2),cex=1.1,cex.axis=1)
ax <- seq(0,12,by=1)
barplot(apply(outz2[(ex.yrs[1]-1994),ex.mns[1],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[1]-1994),ex.mns[1],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[1]-1994),ex.mns[1],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[1],(ex.yrs[1]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=F,las=2)
title(paste0("Month ",ex.mns[1]),line=-2)

barplot(apply(outz2[(ex.yrs[1]-1994),ex.mns[2],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[1]-1994),ex.mns[2],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[1]-1994),ex.mns[2],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[2],(ex.yrs[1]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=F,las=2)
title(paste0("Month ",ex.mns[2]),line=-2)

barplot(apply(outz2[(ex.yrs[1]-1994),ex.mns[3],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[1]-1994),ex.mns[3],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[1]-1994),ex.mns[3],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[3],(ex.yrs[1]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=F,las=2)
title(paste0("Month ",ex.mns[3]),line=-2)

barplot(apply(outz2[(ex.yrs[1]-1994),ex.mns[4],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[1]-1994),ex.mns[4],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[1]-1994),ex.mns[4],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[4],(ex.yrs[1]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=F,las=2)
title(paste0("Month ",ex.mns[4]),line=-2)

barplot(apply(outz2[(ex.yrs[1]-1994),ex.mns[5],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[1]-1994),ex.mns[5],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[1]-1994),ex.mns[5],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[5],(ex.yrs[1]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=F,las=2)
title(paste0("Month ",ex.mns[5]),line=-2)

barplot(apply(outz2[(ex.yrs[1]-1994),ex.mns[6],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[1]-1994),ex.mns[6],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[1]-1994),ex.mns[6],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[6],(ex.yrs[1]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=F,las=2)
mtext(paste0(ex.yrs[1]), side=4,line=0.5,cex=1.6)
title(paste0("Month ",ex.mns[6]),line=-2)


barplot(apply(outz2[(ex.yrs[2]-1994),ex.mns[1],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[2]-1994),ex.mns[1],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[2]-1994),ex.mns[1],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[1],(ex.yrs[2]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=F,las=2)

barplot(apply(outz2[(ex.yrs[2]-1994),ex.mns[2],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[2]-1994),ex.mns[2],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[2]-1994),ex.mns[2],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[2],(ex.yrs[2]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=F,las=2)

barplot(apply(outz2[(ex.yrs[2]-1994),ex.mns[3],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[2]-1994),ex.mns[3],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[2]-1994),ex.mns[3],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[3],(ex.yrs[2]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=F,las=2)

barplot(apply(outz2[(ex.yrs[2]-1994),ex.mns[4],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[2]-1994),ex.mns[4],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[2]-1994),ex.mns[4],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[4],(ex.yrs[2]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=F,las=2)

barplot(apply(outz2[(ex.yrs[2]-1994),ex.mns[5],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[2]-1994),ex.mns[5],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[2]-1994),ex.mns[5],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[5],(ex.yrs[2]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=F,las=2)

barplot(apply(outz2[(ex.yrs[2]-1994),ex.mns[6],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[2]-1994),ex.mns[6],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[2]-1994),ex.mns[6],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[6],(ex.yrs[2]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=F,las=2)
mtext(paste0(ex.yrs[2]), side=4,line=0.5,cex=1.6)


barplot(apply(outz2[(ex.yrs[3]-1994),ex.mns[1],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[3]-1994),ex.mns[1],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[3]-1994),ex.mns[1],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[1],(ex.yrs[3]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=str.names,las=2)

barplot(apply(outz2[(ex.yrs[3]-1994),ex.mns[2],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[3]-1994),ex.mns[2],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[3]-1994),ex.mns[2],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[2],(ex.yrs[3]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=str.names,las=2)

barplot(apply(outz2[(ex.yrs[3]-1994),ex.mns[3],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[3]-1994),ex.mns[3],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[3]-1994),ex.mns[3],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[3],(ex.yrs[3]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=str.names,las=2)

barplot(apply(outz2[(ex.yrs[3]-1994),ex.mns[4],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[3]-1994),ex.mns[4],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[3]-1994),ex.mns[4],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[4],(ex.yrs[3]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=str.names,las=2)

barplot(apply(outz2[(ex.yrs[3]-1994),ex.mns[5],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[3]-1994),ex.mns[5],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[3]-1994),ex.mns[5],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[5],(ex.yrs[3]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=str.names,las=2)

barplot(apply(outz2[(ex.yrs[3]-1994),ex.mns[6],,],1,median,na.rm=T),col=c1,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[3]-1994),ex.mns[6],,],1,quantile,probs=c(0.025),na.rm=T),pch=2)
points(seq(0.8,14,by=1.18),apply(outz2[(ex.yrs[3]-1994),ex.mns[6],,],1,quantile,probs=c(0.975),na.rm=T),pch=3)
barplot(Dist[,ex.mns[6],(ex.yrs[3]-1994)],col=c2,add=T,ylim=c(0,0.55),main='',xlab='',ylab='',xaxt='n')
axis(side = 1,at=seq(0.8,14,by=1.18),labels=str.names,las=2)
mtext(paste0(ex.yrs[3]), side=4,line=0.5,cex=1.6)

mtext('Proportion of population', outer=T, side=2,line=1.5,cex=1.6)
mtext('Spatial strata', outer=T, side=1,line=2.5,cex=1.6)
