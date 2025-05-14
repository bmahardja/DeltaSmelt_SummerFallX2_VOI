##################################################################
### IBMR population dynamics model ###############################
### after Rose et al. 2013a ######################################
### model of growth, reproduction, mortality and movement ########
### William Smith (USFWS; BDFWO); 21 June 2022 ###################

# v2: added place to apply natural mortality/contaminant actions M.act.mult
# -added place to substitute new years of DS.dist observations DS.dist.act
# -fixed stop function to break loop when pop. explodes, before R gets stuck
# v2.2: fixed Move.manual.action function to replace observed data
# v2.3: added PD multipler by prey type

# 1. Initialize
# 2. Get data
# 3. Spawning model
# 4. Movement model
# 5. Consumption and growth model
# 6. Get summaries

######### Bioenergetics Model starts here ########################
delta_smelt_run<-function(){

### 1. Initialize ###
 Wt<-ad.Wt <- rlnorm(super.ad,0.3,0.2)
 L<-ad.L <- (Wt/ln.a[2])^(1/ln.b[2])
 spawn<-ad.spawn <-rbinom(length(L),1,1/(1+exp(-0.25*(L-60))))
 females <- ad.females <- rbinom(length(L),1,0.48)
 YrlyTotls<-TMM_super<-TNS_super<-FMWT_super<-SKT_super<-meanL<-ET1<-ET2<-ET3<-ET4<-ET5<-ET6<-ET7<-mn.L<-lambdaAB4<-NULL
 sp.dist<-array(NA,dim=c((length(yearz)-first+1),12,12)) # year x month x strata
 storlen<-array(NA,dim=c((length(yearz)-first+1),12)) # year x month
 sumz1<-sumz2<-sumz3<-sumz4<-NULL

### 2. Get data ###
 # Get OMR data
 OMR <- make.OMR(first)

 # Get Secchi and Temp
 WQ_out <- make.WQ(first)

 strata <- ad.strata <- rcat(length(ad.L),(smlt.dist(2,first)))

 options(warn=1) # error checker need to set to zero which is default for regular run
 stop=F
 for(yr in c(rep(c(17,18),2),first:length(yearz))) { # year loop
  for(t in c(2:12,1)){ # month loop
   dayz=days.in.month[yr,(t+1)] # number of days in a month
   
   # Get physical data
   temp.dat <- make.temp(t,yr)
   secchi.dat <- WQ_out[yr,t,,1]
   X2.dat <- make.X2(t,yr)
   
   # Get management effects by replacing with a new value
   if(Temp.action == T) {
   if(yr>first & t==1) {
     temp.dat[,1] <- ifelse(is.na(Temp.act.array[(yr-(first-1)+1),t,]),temp.dat[,1],Temp.act.array[(yr-(first-1)+1),t,])
     } else {
     temp.dat[,1] <- ifelse(is.na(Temp.act.array[(yr-(first-1)),t,]),temp.dat[,1],Temp.act.array[(yr-(first-1)),t,])
     }}

   if(Secchi.action == T) {
   if(yr>first & t==1) {
     secchi.dat <- ifelse(is.na(secchi.act.array[(yr-(first-1)+1),t,]),secchi.dat,secchi.act.array[(yr-(first-1)+1),t,])
     } else {
     secchi.dat <- ifelse(is.na(secchi.act.array[(yr-(first-1)),t,]),secchi.dat,secchi.act.array[(yr-(first-1)),t,])
     }}

  if(OMR.action == T) { 
  if(t==1) { # January is last month of each cohort, so use calendar year+1 in Jan
    if(is.na(OMR.act[(yr-(first-1)+1),t])) {
     OMR[yr,t] <- OMR[yr,t]
      } else {
	   OMR[yr,t] <- OMR.act[(yr-(first-1)+1),t]
        }} 
    else {
	 if(is.na(OMR.act[yr-(first-1),t])) {
	  OMR[yr,t] <- OMR[yr,t]
     } else {
       OMR[yr,t] <- OMR.act[yr-(first-1),t]
	   }}}
	
### 3. Spawning model ###
   if(t==2 | t==3 | t==4){
    adults = length(ad.L)
    #fecund <-175.4*exp(ad.L/28.3) # from Rose (Bennett 2005?)
	if(t == 4){
	 fecund <- ifelse(rnorm(adults,temp.dat[ad.strata,1],temp.dat[ad.strata,2])<sp.temp.hi,1,0)*egg.surv*0.0183*((ad.Wt/ln.a[2])^(1/ln.b[2]))^2.7123 # from Damon et al 2017
	 } else {
      fecund <- egg.surv*0.0183*((ad.Wt/ln.a[2])^(1/ln.b[2]))^2.7123 # from Damon et al 2017
	  }
    eggcase<-rep(0,12)
    Eggs<-round(tapply(ad.females*fecund*ad.spawn,ad.strata,sum))
    eggcase[as.numeric(names(Eggs))]<-Eggs
	if(sum(eggcase)> crash.max) {stop=T; crash.yr=(length(SKT_super)+1); break} # stops processing when popn explodes or totally crashes
    if(t == 2){ stage=rep(5,length(L)) }
    }
		
   # assign life stage based on length
   stage <- ifelse(stage ==5,5,ifelse(L < stage.limit[1], 2, ifelse(L < stage.limit[2],3, 4)))
      
### 4. Movement model ###
   if(is.na(mean(L))) {stop=T; crash.yr=(length(SKT_super)+1); break} # stops processing when popn explodes or totally crashes
   last.strata<-strata
   exp.dist <- as.matrix(smlt.dist(t,yr))
   
   if (Move.manual.action == T) {
    if (t==1) {
	 for (j in 1:n.strata) {
      exp.dist[j,1] <- ifelse(is.na(Dist.act.array[(yr-(first-1)+1),t,j]),exp.dist[j,1],Dist.act.array[(yr-(first-1)+1),t,j])
	  }} else {
	  for (j in 1:n.strata) {
      exp.dist[j,1] <- ifelse(is.na(Dist.act.array[(yr-(first-1)),t,j]),exp.dist[j,1],Dist.act.array[(yr-(first-1)),t,j])
	  }}}
    
   new.dist <- matrix(NA,length(L),12)
   exp.dist[,1] <- pmax(1e-3,exp.dist[,1])
   new.dist <- t(exp.dist[,1]*t(move.matrix[last.strata,]))
   new.dist <- rdirichlet(length(L),10*new.dist)
   strata <- rcat(length(L),new.dist)
	  
   ############## Mortality model ############## 
   # Entrainment
   ad.E.risk <- ifelse(OMR[yr,t]<=E.ad.x[1],max.E[1],a.E[1]+b.E[1]*OMR[yr,t])
   ad.E.risk <- ifelse(OMR[yr,t]>=E.ad.x[2],0,ad.E.risk)
   pl.E.risk <- ifelse(OMR[yr,t]<=E.pl.x[1],max.E[2],a.E[2]+b.E[2]*OMR[yr,t])
   pl.E.risk <- ifelse(OMR[yr,t]>=E.pl.x[2],0,pl.E.risk)

   E <- ifelse(L<entrain.L,pl.E.risk,ad.E.risk)
   
   # Predation
   max.turb.fx <- ifelse(L<20,min1[1],ifelse(L>45,min1[2],turb.L.fx.mod$coefficients[1]+turb.L.fx.mod$coefficients[2]*L))
   pred <- (pcr-pcr*max.turb.fx)/(1+exp(-(b.turb[2]*(secchi.dat[strata]-a.turb))))+(1-(pcr-pcr*max.turb.fx))
   pr.surv <- exp(-dayz*(E*entrain.strata[strata]*entrain.month[t]+pred*Mmult*(-0.034+0.165*L^-0.322))) # length-based M from Rose et al. (2013b)
   pr.surv<-ifelse(pr.surv>0.999,0.999,pr.surv) # prevent impossible values of survival
   pr.surv<-ifelse(t>4 & stage==5,0,pr.surv) # kill adults after April
   
### 5. Consumption and growth model ###
   # get food values
   zoop<- make.food(t,yr)
   PD <- zoop[strata,,1]
   for (j in 1:n.prey) {
    PD[,j] <- PD[,j]*rbern(1,(1-zoop[strata,j,3]))
	}
	
   if(Food.spp.action == T) { #multiply food by PD.mult 
   for (j in 1:n.prey) {
    PD[,j] <- PD[,j]*ifelse(is.na(PD.mult_spp[(yr-(first-1)),t,strata,j]),1,PD.mult_spp[(yr-(first-1)),t,strata,j])
	}}
	
   # Bioenergetics model
   rnd.temp.dat <- vector() # make random summer temps
   if(t==7 | t==8 | t==9) {
     rnd.temp.dat <- rnorm(length(L),temp.dat[strata,1],temp.dat[strata,2])
	 } else {
	  rnd.temp.dat <- temp.dat[strata,1]
	  }
   R <- a.r[stage]*(Wt^b.r[stage])*exp(R.Q[stage]*rnd.temp.dat) # metabolism
   # Maximum consumption
   L.1<- exp((1/(T.0[stage]-CQ[stage]))*log(0.98*(1-CK.1[stage])/(CK.1[stage]*0.02))*(rnd.temp.dat-CQ[stage]))
   L.2 <- exp((1/(T.L[stage]-T.M[stage]))*log(0.98*(1-CK.4[stage])/(CK.4[stage]*0.02))*(T.L[stage]-rnd.temp.dat))
   K.A <- CK.1[stage]*L.1/(1+CK.1[stage]*(L.1-1))
   K.B <- CK.4[stage]*L.2/(1+CK.4[stage]*(L.2-1))
   temp.fx <- K.A*K.B
   turb.fx <- (1-max.turb.fx)/(1+exp(-(b.turb[1]*(secchi.dat[strata]-a.turb))))+max.turb.fx
   Cmax <- a.c[stage]*(Wt^b.c[stage])*temp.fx*turb.fx

   Food <- PD*t(V[,stage]/K[,stage])
   C.prey <- (Cmax*Food) /(1+rowSums(Food)) # total consumption of each prey type. function of max consumption, food
   energy.prey <- e.d*C.prey # total energy consumed, weighted by fraction of Limniothona 
   Limno <- energy.prey[,1]/rowSums(energy.prey) # fraction energy from Limno prey
   energy <- e.d[1]*Limno + e.d[2]*(1-Limno) # weighted average energy density
   C <- rowSums(C.prey) # realized consumption rate
   Feg <- F.a[stage]*C # egestion
   U <- U.a[stage]*(C-Feg) # excretion
   SDA <- S.d[stage]*(C-Feg) # specific dynamic action
   SL <- (Sp*spawn*Wt)/dayz # spawning weight loss
   Wt.gain <- Wt*(energy/e.s)*(C-R-Feg-U-SDA)-SL
   Wt.gain <- Wt.gain*dayz # scale weight gain by number days in month
   Wt.last <- Wt
   L.last <- L
   Wt<- Wt+Wt.gain # add the growth
   predL<-(Wt/ln.a[2])^(1/ln.b[2])
   L <- ifelse(L<predL,predL,L)
   
   skinny <- length(L[Wt <= stv*ln.a[2]*L^ln.b[2]])/length(L)
   z <- rbern(length(L),pr.surv) # stochastic mortality; 1=live, 0=die; all skinny fish die
   z <- ifelse(Wt <= stv*ln.a[2]*L^ln.b[2],0,z)

### 6. Get summaries ###
   # alive at beginning of month
   alive <- length(L)
   alive.strata <- strata

   # get mean entrainment risk
   if(t==12) { mean.E1 <- mean(dayz*(E*entrain.strata[strata]*entrain.month[t])*(1-exp(-dayz*(E*entrain.strata[strata]*entrain.month[t]+pred*Mmult*(-0.034+0.165*L^-0.322))))/(dayz*(E*entrain.strata[strata]*entrain.month[t]+pred*Mmult*(-0.034+0.165*L^-0.322)))) }
   if(t==1) { mean.E2 <- mean(dayz*(E*entrain.strata[strata]*entrain.month[t])*(1-exp(-dayz*(E*entrain.strata[strata]*entrain.month[t]+pred*Mmult*(-0.034+0.165*L^-0.322))))/(dayz*(E*entrain.strata[strata]*entrain.month[t]+pred*Mmult*(-0.034+0.165*L^-0.322)))) }
   if(t==2) { mean.E3 <- mean(dayz*(E*entrain.strata[strata]*entrain.month[t])*(1-exp(-dayz*(E*entrain.strata[strata]*entrain.month[t]+pred*Mmult*(-0.034+0.165*L^-0.322))))/(dayz*(E*entrain.strata[strata]*entrain.month[t]+pred*Mmult*(-0.034+0.165*L^-0.322)))) }
   if(t==3) { mean.E4 <- mean(dayz*(E*entrain.strata[strata]*entrain.month[t])*(1-exp(-dayz*(E*entrain.strata[strata]*entrain.month[t]+pred*Mmult*(-0.034+0.165*L^-0.322))))/(dayz*(E*entrain.strata[strata]*entrain.month[t]+pred*Mmult*(-0.034+0.165*L^-0.322)))) }
   if(t==4) { mean.E5 <- mean(dayz*(E*entrain.strata[strata]*entrain.month[t])*(1-exp(-dayz*(E*entrain.strata[strata]*entrain.month[t]+pred*Mmult*(-0.034+0.165*L^-0.322))))/(dayz*(E*entrain.strata[strata]*entrain.month[t]+pred*Mmult*(-0.034+0.165*L^-0.322)))) }
   if(t==5) { mean.E6 <- mean(dayz*(E*entrain.strata[strata]*entrain.month[t])*(1-exp(-dayz*(E*entrain.strata[strata]*entrain.month[t]+pred*Mmult*(-0.034+0.165*L^-0.322))))/(dayz*(E*entrain.strata[strata]*entrain.month[t]+pred*Mmult*(-0.034+0.165*L^-0.322)))) }
   if(t==6) { mean.E7 <- mean(dayz*(E*entrain.strata[strata]*entrain.month[t])*(1-exp(-dayz*(E*entrain.strata[strata]*entrain.month[t]+pred*Mmult*(-0.034+0.165*L^-0.322))))/(dayz*(E*entrain.strata[strata]*entrain.month[t]+pred*Mmult*(-0.034+0.165*L^-0.322)))) }
   
   # keep only the survivors
   Wt<-Wt[z==1]; L<-L[z==1]; stage<-stage[z==1]; strata<-strata[z==1]; females<-females[z==1]
   
   # egg to larvae transition, add to population
   if(t==2 | t==3 | t==4){
    larvae<-rbinom(12,eggcase,egg2larv(TempCv=temp.dat[,1]))
   
    Lrv<-c(rep((5.92-0.05*temp.dat[1,1]),larvae[1]),rep((5.92-0.05*temp.dat[2,1]),larvae[2]),rep((5.92-0.05*temp.dat[3,1]),larvae[3]),
	 rep((5.92-0.05*temp.dat[4,1]),larvae[4]),rep((5.92-0.05*temp.dat[5,1]),larvae[5]),rep((5.92-0.05*temp.dat[6,1]),larvae[6]),
     rep((5.92-0.05*temp.dat[7,1]),larvae[7]),rep((5.92-0.05*temp.dat[8,1]),larvae[8]),rep((5.92-0.05*temp.dat[9,1]),larvae[9]),
	 rep((5.92-0.05*temp.dat[10,1]),larvae[10]),rep((5.92-0.05*temp.dat[11,1]),larvae[11]),rep((5.92-0.05*temp.dat[12,1]),larvae[12]))
		
    WtLrv<- ln.a[1]*Lrv^ln.b[1]

    strataLrv<-c(rep(1,larvae[1]),rep(2,larvae[2]),rep(3,larvae[3]),
	 rep(4,larvae[4]),rep(5,larvae[5]),rep(6,larvae[6]),
     rep(7,larvae[7]),rep(8,larvae[8]),rep(9,larvae[9]),
	 rep(10,larvae[10]),rep(11,larvae[11]),rep(12,larvae[12]))
    femalesLrv <- rbinom(length(Lrv),1,0.48) # sex ratio 48% female, 52% male from SKT 2002-2018 using stages 4 and 5 only
   
    L<-c(L,Lrv) # add larvae to existing population
    Wt<-c(Wt,WtLrv)
    strata<-c(strata,strataLrv)
    stage<-c(stage,rep(1,length(Lrv)))
    spawn<-c(spawn,rep(0,length(Lrv)))
    females<- c(females,femalesLrv)
    }

   # calculate summary statistics
   hist.sim<-hist(strata,plot=F,breaks=seq(0,12,by=1))
   sp.prop<-hist.sim$counts/sum(hist.sim$counts)
   if(t==6) { sumz1=length(L) } # June AB
   if(t==8) { sumz2=length(L) } # Aug AB
   if(t==11){ sumz3=length(L) } # Nov AB
   if(length(L) < 2 | length(L)> crash.max) {stop=T; crash.yr=(length(SKT_super)+1); break} # stops processing when popn explodes or totally crashes

   # lets see what is going on each month
   # print(c(yr,t,length(L),mean(pr.surv),skinny,Limno),digits=3)
   
   # reduce pre-spawn adults to super-individuals using kernel density model
    if(t==1) {
	YrlyTotls<-rbind(YrlyTotls,c(yearz[yr],sum(Wt),length(L)))
	yr.L<-mean(L)
	sumz4=length(L)
    dat <- cbind(Wt,L)
	f1 <- kde(dat,gridsize=n.grid,xmin=apply(dat,2,min),xmax=apply(dat,2,max))
    joint.samp <- rkde(super.ad,f1)
    for (h in 1:clean.cyc) { # reject and replace out-of-range samples
     for (i in 1:nrow(joint.samp)) {
	  if(joint.samp[i,1] < min(Wt)) { joint.samp[i,] <- rkde(1,f1) }
	  if(joint.samp[i,1] > max(Wt)) { joint.samp[i,] <- rkde(1,f1) }
	  if(joint.samp[i,2] < min(L)) { joint.samp[i,] <- rkde(1,f1) }
	  if(joint.samp[i,2] > max(L)) { joint.samp[i,] <- rkde(1,f1) }
	  }}
	
    Wt <- joint.samp[,1]
	L <- joint.samp[,2]
	stage <- rep(4,length(L)) 
	females <- rbern(super.ad,0.48)
	strata <- rcat(length(L),sp.prop)
	}

   # sexual maturation at t+1
   if(t>=4) { spawn<-rep(0,length(L)) } # prevent spawning penalty outside of spawning season
   if(t<4) { spawn<-rbinom(length(L),1,1/(1+exp(-0.25*(L-60))))
    ad.L<-L[L>45]
    ad.Wt<-Wt[L>45]
    ad.strata<-strata[L>45]
    ad.spawn<-spawn[L>45]
    ad.females<-females[L>45] }
	else ad.L<-ad.Wt<-ad.strata<-ad.spawn<-ad.females<-NULL
   
   sp.dist[(yr-15),t,]<-sp.prop # store spatial distributions
   storlen[(yr-15),t]<-mean(L) # store mean monthly lengths for annual mean growth trajectories
    }
   lambdaAB4<-c(lambdaAB4,sumz4/super.ad)
   TMM_super<-c(TMM_super,sumz1)
   TNS_super<-c(TNS_super,sumz2)
   FMWT_super<-c(FMWT_super,sumz3)
   SKT_super<-c(SKT_super,sumz4)
   meanL<-c(meanL,yr.L)
   ET1<-c(ET1,mean.E1[1])
   ET2<-c(ET2,mean.E2[1])
   ET3<-c(ET3,mean.E3[1])
   ET4<-c(ET4,mean.E4[1])
   ET5<-c(ET5,mean.E5[1])
   ET6<-c(ET6,mean.E6[1])
   ET7<-c(ET7,mean.E7[1])
   
   if (stop==T) {TMM_super<-c(TMM_super[1:(crash.yr-1)],rep(NA,length(crash.yr:((length(yearz)-first+1)+4))))
    TNS_super<-c(TNS_super[1:(crash.yr-1)],rep(NA,length(crash.yr:((length(yearz)-first+1)+4))))
    FMWT_super<-c(FMWT_super[1:(crash.yr-1)],rep(NA,length(crash.yr:((length(yearz)-first+1)+4))))
	SKT_super<-c(SKT_super[1:(crash.yr-1)],rep(NA,length(crash.yr:((length(yearz)-first+1)+4))))
    meanL<-c(meanL[1:(crash.yr-1)],rep(NA,length(crash.yr:((length(yearz)-first+1)+4))))
    ET1<-c(ET1[1:(crash.yr-1)],rep(NA,length(crash.yr:((length(yearz)-first+1)+4))))
    ET2<-c(ET2[1:(crash.yr-1)],rep(NA,length(crash.yr:((length(yearz)-first+1)+4))))
    ET3<-c(ET3[1:(crash.yr-1)],rep(NA,length(crash.yr:((length(yearz)-first+1)+4))))
    ET4<-c(ET4[1:(crash.yr-1)],rep(NA,length(crash.yr:((length(yearz)-first+1)+4))))
    ET5<-c(ET5[1:(crash.yr-1)],rep(NA,length(crash.yr:((length(yearz)-first+1)+4))))
	ET6<-c(ET6[1:(crash.yr-1)],rep(NA,length(crash.yr:((length(yearz)-first+1)+4))))
	ET7<-c(ET7[1:(crash.yr-1)],rep(NA,length(crash.yr:((length(yearz)-first+1)+4))))
    break}
    }
	
   # Estimate abundance indices
   lambdaAB1<-lambdaAB2<-lambdaAB3<-est_TMM<-est_TNS<-est_FMWT<-est_SKT<-vector(length=24)
   lambdaAB1[1]<-lambdaAB2[1]<-lambdaAB3[1]<-NA 
   est_TMM[1]<-TMM_super[1]
   est_TNS[1]<-TNS_super[1]
   est_FMWT[1]<-FMWT_super[1]
   est_SKT[1]<-SKT_super[1]
   for (y in 2:24) {
   lambdaAB1[y]<-(SKT_super[y-1]/TMM_super[y-1])*(TMM_super[y]/super.ad)
   lambdaAB2[y]<-(SKT_super[y-1]/TNS_super[y-1])*(TNS_super[y]/super.ad)
   lambdaAB3[y]<-(SKT_super[y-1]/FMWT_super[y-1])*(FMWT_super[y]/super.ad)
   est_TMM[y]<-est_TMM[y-1]*lambdaAB1[y]
   est_TNS[y]<-est_TNS[y-1]*lambdaAB2[y]
   est_FMWT[y]<-est_FMWT[y-1]*lambdaAB3[y]
   est_SKT[y]<-est_SKT[y-1]*lambdaAB4[y]
   }
   
  est_all<-cbind(est_TMM,est_TNS,est_FMWT,est_SKT,meanL,ET1,ET2,ET3,ET4,ET5,ET6,ET7)
  return(list(est_all,sp.dist))
 } 
