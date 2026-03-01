rm(list=ls())
library(popbio)
library(matrixcalc)
library(readxl)
library(forecast)
library(fitdistrplus)
library(statmod)
library(viridis)
library(gtools)
options(max.print=5000)
setwd('C:\\Users\\fjohn\\OneDrive\\Documents\\PROJECTS\\Iceland\\Puffin\\Modeling\\2025 work')
source('C:\\Users\\fjohn\\OneDrive\\Documents\\MyFunctions.R')

#########################################################################
#                        PRELIMINARIES
#########################################################################
# sea temperatures
#########################################################################
sst <- read.csv('Sep-Aug SST.csv',header = F)
colnames(sst)<-c('year','sst')
sst <- sst[-nrow(sst),]
acf(sst$sst,lag.max=80)
data <- ts(sst$sst)
full_model <- auto.arima(data,stepwise=FALSE)
summary(full_model) # non-stationary
fit <- Arima(sst$sst,order=c(3,1,0),include.drift = TRUE)
summary(fit)
start <- 7
# a <- arima.sim(model=list(order=c(3,1,0),ar=fit$coef[1:3]),mean = fit$coef[4],n=80) + start
# plot(a,type='l')

cast <- forecast(fit,h=40)
series <- c(cast$model$x,cast$mean)
res <- 300
#dev.new(height=4*res,width=6*res,noRStudioGD = TRUE)
par(mar=c(5,6,2,2),cex.lab=1.5)
plot(cast$model$x,las=1,ylab='SST',xlab='Year',cex.lab=2,xaxt='n',main='',lwd=2,type='l',xlim=c(1,186),ylim=c(5,10))
lines(147:(147+39),cast$mean,col='blue',lwd=3)
lines(147:(147+39),cast$lower[,1],col='blue',lwd=2,lty=3)
lines(147:(147+39),cast$upper[,1],col='blue',lwd=2,lty=3)
abline(h=seq(5.5,10.5,0.5),v=seq(3,200,20),lty=2,col='lightgray')
lbl <- as.character(seq(1880,(1880+145+40),5))
axis(1,at=seq(3,(3+145+40),5),labels=lbl)
abline(h=c(6.2,7.7),lty=2,lwd=2)

sst.sub <- subset(sst,year>=1920) # eliminate early very cold period
data <- ts(sst.sub$sst)
best_model <- auto.arima(data,stepwise=FALSE) # now is stationary
summary(best_model)
coefs <- best_model$coef
sd.best_model <- sqrt(best_model$sigma2)
mean.sst <- mean(sst.sub$sst)
summary(sst.sub$sst)
p <- 7
hunk <- split(sst.sub$sst[-(1:6)], rep(1:(length(sst.sub$sst[-(1:6)])/round(p)), each = round(p))) # get SST in blocks of G years
CV<-NULL
for (i in 1:length(hunk)) CV[i]<- sd(hunk[[i]])/mean(hunk[[i]])
(cv <- mean(CV))

res <- 300
#dev.new(height=4*res,width=6*res,noRStudioGD = TRUE)
par(mar=c(5,6,2,2))
hist(sst.sub$sst,breaks=seq(6,9,0.25),main='',xlab='Sea surface temperature',las=1,cex.lab=1.5); box()
abline(v=c(6.1,7.8),lty=2,lwd=2) # temps of positive population growth
abline(v=6.966655,lty=3,lwd=2)

#dev.new(height=4*res,width=6*res,noRStudioGD = TRUE)
par(mar=c(5,6,2,2))
loess_model <- loess(sst ~ year, data = sst,span=0.25)
years<-1878+(0:145)
plot(years,sst$sst, type='l',las=1,ylim=c(5,9),
     xlab='Year',ylab='SST',xaxp=c(1875,2025,15),lwd=1,cex.lab=1.5)
abline(h=5:9,v=seq(1875,2025,10),lty=2,col='lightgray')
lines(sst$year,predict(loess_model),lwd=3,col='black')
lines(1920:2023,best_model$fitted,type='l',col='blue',lwd=3)


#########################################################################
# Relative production function based on Hansen et al 2021 in the Westmans
#########################################################################
pb3 <- pb4 <- 0.067; pb5 <- 0.75; pb6 <- 1 # probability of breeding from Erpur

Praw <- read_xlsx('gcb15665-sup-0002-datas1.xlsx',na='NA',sheet='Data',col_names=TRUE) 
chkndx <- subset(Praw, Year>=1878 & Year<=2005)
ppi <- 300
#dev.new(height=6*ppi,width=9*ppi,noRStudioGD = TRUE)
 par(mar=c(5,6,2,2))
 with(chkndx,plot(SST,(Prod/6.6),xlim=c(6,8.5),ylim=c(0,1.5),ylab='Relative production',
                 xlab='SST',las=1,cex.lab=1.5)) # normalize by max modeled Prod = 6.6 from Erpur 2021

#... Fit normal kernel 
  prod.norm <- chkndx$Prod/6.6
  SST.prod <- chkndx$SST
  fit <- nls( prod.norm ~ exp(-0.5 * ((SST.prod - Topt)/Tsd)^2),
           start = list( Topt = 7, Tsd = 1))
  summary(fit)
  se.p <- summary(fit)$sigma
  newSST <- data.frame(SST.prod = seq(6,8.6,0.1))
  yfit <- predict(fit, newdata = newSST)
  lines(newSST$SST.prod, yfit, col = "blue", lwd = 3)

# relative production fcn from fitted kernel
  Topt <- fit$m$getPars()[1]
  Tsd <-  fit$m$getPars()[2]
  rP <- function(temp,size,se=se.p) exp(-0.5 * ((temp-Topt)/Tsd)^2)  + rnorm(size,0,se) 
  
# production data to get maximum Iceland productivity
  Pobs <- read.csv('ProdMonitor.csv',header=TRUE)
  
  # Split the data by group
  split_data <- split(Pobs, Pobs$Location)

  library(dplyr)
  library(tidyr)
  
  cor_matrix <- cor(
    tidyr::pivot_wider(
      dplyr::bind_rows(
        lapply(split_data, function(x) dplyr::select(x, Year, NS)),
        .id = "Location"
      ),
      names_from  = Location,
      values_from = NS
    )[ , -1],
    use = "pairwise.complete.obs"
  )

  loc_order <- c("Akurey", "Elli?aey","G_Stgrf", "Vigur",  "Drangey", "Lundey",
                "Grimsey", "Hholmi",  "Papey", "Ihofdi" , "Dyrhol", "Eyjar")
  cor_matrix <- cor_matrix[loc_order, loc_order]
  
#  install.packages("corrplot")
  library(corrplot)
  
  corrplot(cor_matrix, method = "color", type = "upper",
           tl.col = "black", tl.cex = 0.9)  
  
    
  # Set up the plot
#  dev.new(height=4*res,width=6*res,noRStudioGD = TRUE)
  par(mar=c(5,6,2,2))
  plot(Pobs$Year, Pobs$P, type = "n", xlab = "Year",
       ylab = "Chick production (mean-centered)",cex.lab=1.5,
       las=1)  # ,ylim=c(-0.5,0.5))
  grid()
  # Add lines for each group
  color  <- rep(c("#000000", "#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),2)
  ndx <- c(4,5,7,11)
  # ndx<-1:12
  for (i in 1:12) {
    if (i %in% ndx) lines(split_data[[i]]$Year, split_data[[i]]$P # - mean(split_data[[i]]$P
                                                                         ,na.rm=TRUE #)
                          , col = color[i],
                          type='b',pch=16,lwd=3)
  }
  cnames <- names(split_data)[ndx]; cnames[1]<-'Elli?aey'
  legend("topleft", legend = cnames, col = color[ndx],
         lty = 1,pch=16,bty='n')
  abline(h=0,lwd=2,lty=1)

# use burrows for weighting  
  holes.loc <- data.frame(unique(Pobs$Location))
  holes.loc$burrow <- c(177.0,10.1,45.2,31.0,99.9,36.5,38.4,117.0,NA,425.0,
                       1125.0,NA)
  colnames(holes.loc) <- c('loc','burrow')
  holes.loc.sub <- subset(holes.loc,burrow>0)
  Pobs.sub <- Pobs[Pobs$Location %in% holes.loc.sub$loc, ]
  max.pobs <- NULL
  max.pobs[1] <- max(subset(Pobs.sub,Location=='Papey')$P,na.rm=TRUE)
  max.pobs[2] <- max(subset(Pobs.sub,Location=='Hholmi')$P,na.rm=TRUE)
  max.pobs[3] <- max(subset(Pobs.sub,Location=='Drangey')$P,na.rm=TRUE)
  max.pobs[4] <- max(subset(Pobs.sub,Location=='G_Stgrf')$P,na.rm=TRUE)
  max.pobs[5] <- max(subset(Pobs.sub,Location=='Grimsey')$P,na.rm=TRUE)
  max.pobs[6] <- max(subset(Pobs.sub,Location=='Lundey')$P,na.rm=TRUE)
  max.pobs[7] <- max(subset(Pobs.sub,Location=='Vigur')$P,na.rm=TRUE)
  max.pobs[8] <- max(subset(Pobs.sub,Location=='Akurey')$P,na.rm=TRUE)
  max.pobs[9] <- max(subset(Pobs.sub,Location=='Elli?aey')$P,na.rm=TRUE)
  max.pobs[10] <- max(subset(Pobs.sub,Location=='Eyjar')$P,na.rm=TRUE)
  (maxP <- weighted.mean(max.pobs,holes.loc.sub$burrow))
  max.pobs <- max.pobs[c(8,9,4,7,3,6,5,2,1,10)]
    
  mean.BO <- NULL
  #xyears <- subset(Pobs.sub,Year>=2011 & Year<=2020)
  xyears <- subset(Pobs.sub,Year>=2020 & Year<=2023)
  mean.BO[1] <- mean(subset(xyears,Location=='Papey')$BO,na.rm=TRUE)
  mean.BO[2] <- mean(subset(xyears,Location=='Hholmi')$BO,na.rm=TRUE)
  mean.BO[3] <- mean(subset(xyears,Location=='Drangey')$BO,na.rm=TRUE)
  mean.BO[4] <- mean(subset(xyears,Location=='G_Stgrf')$BO,na.rm=TRUE)
  mean.BO[5] <- mean(subset(xyears,Location=='Gr?msey')$BO,na.rm=TRUE)
  mean.BO[6] <- mean(subset(xyears,Location=='Lundey')$BO,na.rm=TRUE)
  mean.BO[7] <- mean(subset(xyears,Location=='Vigur')$BO,na.rm=TRUE)
  mean.BO[8] <- mean(subset(xyears,Location=='Akurey')$BO,na.rm=TRUE)
  mean.BO[9] <- mean(subset(xyears,Location=='Elli?aey')$BO,na.rm=TRUE)
  mean.BO[10] <- mean(subset(Pobs.sub,Location=='Eyjar')$BO,na.rm=TRUE)
  Nb_percolony <- (mean.BO*holes.loc.sub$burrow)*2
  BP <- sum(Nb_percolony,na.rm=TRUE)
  
  # check predicted production against 2010-2024 burrow monitoring data: Puffin Rally 2010-2024.xlsx
  B <- maxP
  Pobs <- numeric(15)
  for (i in 1:15) {
    Pobs[i] <- weighted.mean(Pobs.sub[Pobs.sub$Year==i+2009,7],holes.loc.sub$burrow,na.rm=TRUE)
  } 
  s <- c(sst.sub[sst.sub$year>=2010,2],NA)
  Ppred <- rP(s,1,0)*B
  
  pmodel <- lm(Ppred~log(Pobs))
  (s <- summary(pmodel))
  pvalue <- s$coefficients["log(Pobs)","Pr(>|t|)"]
  reg.pred <- (s$coefficients[1]+s$coefficients[2]*log(seq(0.005,0.7,0.005)))
  
#  dev.new(height=4*res,width=6*res,noRStudioGD = TRUE)
  par(mar=c(5,6,2,2))
  plot(Pobs,Ppred,xlim=c(0,0.7),ylim=c(0,0.7),las=1,pch=16,cex.lab=1.5);grid();abline(0,1)
  text(Pobs+0.02,Ppred,2010:2024,cex=1.25)
  lines(seq(0.005,0.7,0.005),reg.pred,lty=2)
  text(0.05,0.15,paste("P = ",round(pvalue,2)))
  

# final production model
temp <- seq(6,8.8,length=100)
p <- rP(temp,1,0)*B
plot(temp,p,type='l',lwd=2,ylab='Chicks/pair',xlab='SST',las=1)
grid()

    
#########################################################################
# Survival ##############################################################
#########################################################################
#sa <- 0.92 # natural survival ages 1-6+ from Erpur 2024 status report
#sa <- 0.874 # from Hafnarh?lmi with no hunting Erpur email 9//11/25
# survival from Halfdan thesis N?tt?rufr Lei?b Til_ritryna_2014 + Hafnarh?lmi
  # sa <- c(.975,.924,.935,.805,.93,.935,.935,.902,.882,.943,.971,.886,.843,.804,.935,.946,.868,.859,.881,.95,.874)
  # sa.mu <- mean(sa); sa.sd <- sd(sa); sa.se<-sa.sd/sqrt(length(sa))
  # sa.parm <- fitdist(sa,'beta',method='mle')$estimate
  # beta.par <- MOM.beta(sa.mu,(sa.sd/sqrt(21))^2)
  #sa.nodes <- gauss.quad.prob(4,dist='beta',alpha=beta.par$a,beta=beta.par$b)
  
  # curve(dbeta(x,beta.par$a,beta.par$b), from=0.5, to=1.0,lwd=3,col='blue')
  # par(mar=c(5,6,2,2))
  # hist(sa,breaks=seq(0.5,1,0.05),main='',xlab="Adult survival",freq=F,cex.lab=1.5,xlim=c(0.7,1.0),las=1) # ,add=TRUE) #,col=NULL); 
  # box()
  # curve(dbeta(x,sa.parm[1],sa.parm[2]),col='darkred',lwd=3,add=TRUE)
  # abline(v=mean(sa),lwd=2,lty=2);text(0.88,1,'mean \n(n=21)',cex=1.5)

s <- c(0.5, 0.6, 0.7) # natural survival age 0 to 1 (uncertain)
s1 <- s[3] # as in Carl's work, need this to better match expected growth in absence of harvest
  
# Grosbois, V., M. P. Harris, T. Anker-Nilssen, R. H. McCleery, D. N. Shaw, B. J. T. Morgan, and O. Gimenez. 2009. 
#... Modeling survival at multi‐population scales using mark–recapture data. Ecology 90:2922–2932.
#...  Isle of May (dependency on sand eel)
phi.f <- function(sst,size,vl=0.11) # var on logit scale
  {
  logit_phi <- 2.95 - 0.10*(sst - 7.0) - 0.20*(sst - 7.0)^2 + rnorm(size, 0, sqrt(vl))
  phi <- plogis(logit_phi)
  return(phi)
}
phi.vec <- Vectorize(phi.f)

#.... adjust for Iceland based on 0.874 survival from Erpur at Hafnarh?lmi w/mean 8.05C SST during 2020:2023
s <- mean(sst[143:146,2]) # 2020:2023 mean
sa.H <- 0.874 
pi <- sa.H
p <- phi.f(s,1,0)
prop <- pi/p
phi.fi <- function(sst, size, vl = 0.11, scale = prop) {
  logit_phi <- 2.95 - 0.10*(sst - 7.0) - 0.20*(sst - 7.0)^2 +
    rnorm(size, 0, sqrt(vl))
  phi <- plogis(logit_phi) * scale
  return(phi)
}
phii.vec <- Vectorize(phi.fi)


SST <- seq(6,9,0.1)
yf <- phi.f(SST,1,0)  
yfi <- phi.fi(SST,1,0) 
yfi.rep = phii.vec(SST,10000)
quant.fcn <- function(x) quantile(x,probs=c(0.05,0.90))
yfi.quant <- apply(yfi.rep,2,quant.fcn)

plot(SST,yf,ylab='Survival',type='l',las=1,lwd=3,ylim=c(min(yfi),1))
lines(SST,yfi,lwd=3,col='blue')
lines(SST,yfi.quant[1,],lty=2,col='blue')
lines(SST,yfi.quant[2,],lty=2,col='blue')
grid()
s <- 0.874 
points(8.05,phi.fi(8.05,1,0),pch=16,cex=1.5,col='blue')
text(8.0,0.94,pos=4,'Isle of May')
text(8.0,0.88,pos=4,'Iceland',col='blue')

# no-harvest, mean matrix
B <- maxP
temp <- mean.sst
sa <- phi.f(temp,1,0)
rp <- rP(temp,1,0)

MAY <- matrix(c(0,0,s1*pb3*B*rp/2,s1*pb4*B*rp/2,s1*pb5*B*rp/2,s1*pb6*B*rp/2,
             sa,0,0,0,0,0,
             0,sa,0,0,0,0,
             0,0,sa,0,0,0,
             0,0,0,sa,0,0,
             0,0,0,0,sa,sa),byrow=TRUE,nrow=6)
eMAY <- eigen.analysis(MAY)
eMAY$lambda1

# Iceland
sa <- phi.fi(temp,1,0)

ICE <- matrix(c(0,0,s1*pb3*B*rp/2,s1*pb4*B*rp/2,s1*pb5*B*rp/2,s1*pb6*B*rp/2,
             sa,0,0,0,0,0,
             0,sa,0,0,0,0,
             0,0,sa,0,0,0,
             0,0,0,sa,0,0,
             0,0,0,0,sa,sa),byrow=TRUE,nrow=6)
eICE <- eigen.analysis(ICE)
eICE$lambda1

  
#########################################################################
# Generation time
#########################################################################

# from popbio #############
lambda(ICE)                     # long-term growth rate
net.reproductive.rate(ICE)      # R0
(g = generation.time(ICE))            # log(R0)/log(lambda)


# life table approach ###########
s <- c(0.5, 0.6, 0.7) # natural survival age 0 to 1 (uncertain)
s1 <- s[3] # as in Carl's work, need this to better match expected growth in absence of harvest
pb3 <- pb4 <- 0.067; pb5 <- 0.75; pb6 <- 1 # probability of breeding from Erpur
B <-  maxP # max reproduction per pair (calculated in "puffin harvest mgmt.R")
f <- B * rP(mean.sst,1,0) / 2 # mean repro by breeding adult
rp <- rP(mean.sst,1,0)
sa <- phi.fi(mean.sst,1,0)
#x <- 3:50 # reproductive ages (use enough that lx drops to ~0.01)
end <- 26 #26 #50 for Isle of May; 26 for Iceland
x <- 3:end # 3% rule (Johnson et al. 2012) (depends on value of sa used)
#x <- 3:36 # max reported in captivity (Painter 2022)
(lx <- s1*sa^(x-1)) # survival until age x (1st year = 0.7; sa thereafter)
mx <- c(f*pb3,f*pb4,f*pb5,rep(f*pb6,x[length(x)]-5)) # reproductive rate; no senescence
num <- sum(x*lx*mx) # Caswell 2001: 128
den <- sum(lx*mx)
(G <- num/den)

# no-harvest, mean matrix
A <- matrix(c(0,0,s1*pb3*B*rp/2,s1*pb4*B*rp/2,s1*pb5*B*rp/2,s1*pb6*B*rp/2,
             sa,0,0,0,0,0,
             0,sa,0,0,0,0,
             0,0,sa,0,0,0,
             0,0,0,sa,0,0,
             0,0,0,0,sa,sa),byrow=TRUE,nrow=6)
e <- eigen(A)
ei <- e$values
eA <- eigen.analysis(A)

(P2 <- 2*pi/atan(Im(ei[2])/Re(ei[2]))) # period of oscillations (Caswell2001:100) 
(tx <- log(100)/log(eA$damping.ratio)) # decay time for transient dynamics to 1% of original size
   rho <- Mod(ei[2]) / Mod(ei[1]); (tx2 <-  log(0.05)/log(abs(rho)) )
pb <- c(0,0,pb3,pb4,pb5,pb6)
age <- 1:6
(AFB <- sum(age*pb)/sum(pb)) # age at first breeding
(transient.period <- AFB/log(eA$damping.ratio)) # Caswell 2001 approximation
GP <- c(g,G,P2,tx)

# Table SI Hansen 2021 article (proportions in the harvest)
pSI <- c(0.002,0.163,0.339,0.206,0.080,0.027+.024+.021+.018+.016+.014+.012)
pSI <- c(pSI[1:5],1-sum(pSI[1:5])) # to account for no harvest records of 13+ birds)
v <- pSI/eA$stable.stage 
(v <- round(v/max(v) ,2)) # age-specific vulnerabilities to harvest relative to max(v)
(RVx.unscaled <- eA$repro.value) 
c <- 1/RVx.unscaled[6]
(RVx <- RVx.unscaled*c)
(RVH <- sum(v*RVx)/sum(v)) # mean repro value of a harvested individual
sum(eA$stable.stage*RVx) # mean repro value if no differential vulnerability

#############################################################################
# lambda for a range of se temps
#############################################################################
temp <- seq(5.5,8.5,0.1)
lam <- NULL
for (i in 1:length(temp)) {
    Z <- matrix(c(0,0,s1*pb3*B*rP(temp[i],1,0)/2,s1*pb4*B*rP(temp[i],1,0)/2,s1*pb5*B*rP(temp[i],1,0)/2,s1*pb6*B*rP(temp[i],1,0)/2,
                 phi.fi(temp[i],1,0),0,0,0,0,0,
                 0,phi.fi(temp[i],1,0),0,0,0,0,
                 0,0,phi.fi(temp[i],1,0),0,0,0,
                 0,0,0,phi.fi(temp[i],1,0),0,0,
                 0,0,0,0,phi.fi(temp[i],1,0),phi.fi(temp[i],1,0)),byrow=TRUE,nrow=6)
    lam[i] <- eigen.analysis(Z)$lambda
  }    

plot(temp,lam,xlab='Sea temperature',ylab='Population growth rate',type='b',las=1,lwd=2)
grid()
abline(h=1,lty=2)
abline(v=mean.sst,lty=2)
text(7.9,1.02,pos=4,'mean SST\n(1920 - 2023)')

#########################################################################
# CALCULATE QUOTAS FOR VARIOUS NUMBERs OF BREEDING PAIRS AND MEAN SSTs
#########################################################################

################ deterministic ####################################
# p <- GP[1]
p <- 10 # choose harvest periodicity
goal <- 4000 # desired breeding individuals

Nbps <- seq(1000,3000,100) # breeding pairs
ssts <- seq(6,8.5,0.1)
Hquota <- matrix(nrow=length(Nbps),ncol=length(ssts))

for (j in 1:length(Nbps)) {
  for (i in 1:length(ssts)) {
    Nbp <- Nbps[j] # current number of breeding pairs (thousands)
    mu.sst <- ssts[i] # expected mean over next G or P2 years 
    sa <- phi.fi(mu.sst,1,0) # Icelandic survival
    rp <- rP(mu.sst,1,0)
    Z <- matrix(c(0,0,s1*pb3*B*rp/2,s1*pb4*B*rp/2,s1*pb5*B*rp/2,s1*pb6*B*rp/2,
                 sa,0,0,0,0,0,
                 0,sa,0,0,0,0,
                 0,0,sa,0,0,0,
                 0,0,0,sa,0,0,
                 0,0,0,0,sa,sa),byrow=TRUE,nrow=6)
    eZ <- eigen.analysis(Z)
    #... convert breeding pairs to population vector
    totN <- (Nbp*2)/sum(eZ$stable.stage[3]*pb3,eZ$stable.stage[4]*pb4,eZ$stable.stage[5]*pb5,eZ$stable.stage[6]*pb6)
    init.vec <- eZ$stable.stage*totN
    N <- matrix(nrow=6,ncol=round(p+1)); N[,1]<-init.vec
    for (k in 1:round(p)) { N[,k+1] <- Z%*%N[,k] }
    #... allowable harvest
    NB.end <- N[3,round(p)+1]*pb3 + N[4,round(p)+1]*pb4 + N[5,round(p)+1]*pb5 + N[6,round(p)+1]*pb6
    NB.min <- goal # number of min # of breeding individuals in thousands
    loss.max <- ifelse((NB.end-NB.min)>0,(NB.end-NB.min),0)
    RVx.unscaled <- eZ$repro.value
    c <- 1/RVx.unscaled[6]
    RVx <- RVx.unscaled*c
    RVH <- sum(v*RVx)/sum(v) # mean RV of harvested bird, weighted by vulnerability
    Harvest.total <- loss.max/RVH # allowable harvest in one generation
    Hquota[j,i] <- Harvest.total/round(p) # annual allowable
    
  }
}

# plot policy
n<-Nbps
s<-ssts
z<-matrix(Hquota,nrow=length(s),ncol=length(n),byrow=T)
# ppi <- 300
# dev.new(height=6*ppi,width=9*ppi,noRStudioGD = TRUE)
ppi<-400
hi<-600; inc<-50
Hseq <- seq(0,hi,inc)
clr <- hcl.colors((hi/inc)+1,"Spectral")
#tiff(file='MANUSCRIPT\\Figures\\deterministic policy_7years.tif',res=ppi,width = 9*ppi, height = 6*ppi)
par(mar=c(5,5,5,5))
filled.contour(x=n,y=s,z=t(z),col=clr,levels=Hseq,ylab='Expected mean SST',main='',
               xlab='Current number of breeding pairs (thousands)',cex.lab=1.25,
               key.title=title(main='Annual\nquota\n(thousands)',cex.main=0.75),
               plot.axes = { axis(1, Nbps) 
                 axis(2, seq(6,8.5,0.5)) })
#dev.off()

################ consider risk ####################################
# p <- GP[1]
p <- 10 # choose harvest periodicity
goal <- 4000 # desired breeding individuals

size <- 2000
Nbps <- seq(1000,3000,100) # breeding pairs
ssts <- seq(6,8.5,0.1)
Hquota <- matrix(nrow=length(Nbps),ncol=length(ssts))
quota <- numeric(size)
set.seed(123)
sst.sim <- arima.sim(n = 1000*round(p), list(ar = coefs[grep("ar", names(coefs))],
                                ma = coefs[grep("ma", names(coefs))]), 
                    sd = sd.best_model) + coefs[3]
blocks <- split(sst.sim, rep(1:(length(sst.sim)/round(p)), each = round(p))) # get SST in blocks of G years
CV<-NULL
for (i in 1:length(blocks)) CV[i]<- sd(blocks[[i]])/mean(blocks[[i]])
(cv <- mean(CV)) # this is ~same cv I get when using the actual data, broken into blocks

for (j in 1:length(Nbps)) {

    for (i in 1:length(ssts)) {
     Nbp <- Nbps[j] # current number of breeding pairs (thousands)

     for (l in 1:size) 
     {
     mu.sst <- rnorm(1,ssts[i],ssts[i]*cv)  
     sa <- phi.fi(mu.sst,1)
     rp <- rP(mu.sst,1)
     Z <- matrix(c(0,0,s1*pb3*B*rp/2,s1*pb4*B*rp/2,s1*pb5*B*rp/2,s1*pb6*B*rp/2,
               sa,0,0,0,0,0,
               0,sa,0,0,0,0,
               0,0,sa,0,0,0,
               0,0,0,sa,0,0,
               0,0,0,0,sa,sa),byrow=TRUE,nrow=6)
     eZ <- eigen.analysis(Z)
  #... convert breeding pairs to population vector
     totN <- (Nbp*2)/sum(eZ$stable.stage[3]*pb3,eZ$stable.stage[4]*pb4,eZ$stable.stage[5]*pb5,eZ$stable.stage[6]*pb6)
     init.vec <- eZ$stable.stage*totN
     N <- matrix(nrow=6,ncol=round(p)+1); N[,1]<-init.vec
     for (k in 1:round(p)) { N[,k+1] <- Z%*%N[,k] }
  #... allowable harvest
   NB.end <- N[3,round(p)+1]*pb3 + N[4,round(p)+1]*pb4 + N[5,round(p)+1]*pb5 + N[6,round(p)+1]*pb6
   NB.min <- goal # number of min # of breeding individuals in thousands
   loss.max <- ifelse((NB.end-NB.min)>0,(NB.end-NB.min),0)
   RVx.unscaled <- eZ$repro.value 
   c <- 1/RVx.unscaled[6]
   RVx <- RVx.unscaled*c
   RVH <- sum(v*RVx)/sum(v) # mean RV of harvested bird, weighted by vulnerability
   Harvest.total <- loss.max/RVH # allowable harvest in one generation
   quota[l] <- Harvest.total/round(p) # annual allowable
   }
      
  Hquota[j,i] <- quantile(quota,probs=0.25) 
 }
}  
Hquota.4000.25 <- Hquota

# plot policy
n<-Nbps
s<-ssts
z<-matrix(Hquota,nrow=length(s),ncol=length(n),byrow=T)

ppi<-400
hi<-600; inc<-50
Hseq <- seq(0,hi,inc)
clr <- hcl.colors((hi/inc)+1,"Spectral")
#tiff(file='MANUSCRIPT\\Figures\\Q25 policy_7years.tif',res=ppi,width = 9*ppi, height = 6*ppi)
par(mar=c(5,5,5,5))
filled.contour(x=n,y=s,z=t(z),col=clr,levels=Hseq,ylab='Expected mean SST',main='',
               xlab='Current number of breeding pairs (thousands)',cex.lab=1.25,
               key.title=title(main='Annual\nquota\n(thousands)',cex.main=0.75),
               plot.axes = { axis(1, Nbps) 
                 axis(2, seq(6,8.5,0.50)) })
#dev.off()

#########################################################################
# STOCHASTIC SIMULATION TESTING
#########################################################################
# look up function for quota
Nvals <- seq(1000,3000,100) # rows & columns of saved Hquota matrices
Tvals <- seq(6,8.5,0.1)
snap_quota <- function(N, T, Nvals, Tvals, H) {
  x <- which.min(abs(Nvals - N))
  y <- which.min(abs(Tvals - T))
  round(H[x,y]) # round to nearest thousand
}

m <- 1 # choose period
p <- round(GP[m]) # periodicity
p <- 10

set.seed(5)
n <- 10000
sst.sim <- arima.sim(n = n*p, list(ar = coefs[grep("ar", names(coefs))],
                                  ma = coefs[grep("ma", names(coefs))]), 
                    sd = sd.best_model) + coefs[3]
blocks <- split(sst.sim, rep(1:(length(sst.sim)/round(p)), each = p)) # get SST in blocks of G years

z <- 10 # generations
t <- p * z # total years
sapply(blocks[1:z],mean)

# initial quota
 Q <- Hquota.4000.25
 q <- snap_quota(2000,mean(blocks[[1]]),Nvals,Tvals,Q)
# OR...
#q <- 32

# initial population size based on SAD of mean matrix
N <- matrix(nrow=6,ncol=(t+1))
sa <- mean(phi.fi(mean.sst,2000))
rp <- mean(rP(mean.sst,2000))
M <- matrix(c(0,0, s1*pb3*B*rp/2, s1*pb4*B*rp/2, s1*pb5*B*rp/2, s1*pb6*B*rp/2,
             sa,0,0,0,0,0,
             0,sa,0,0,0,0,
             0,0,sa,0,0,0,
             0,0,0,sa,0,0,
             0,0,0,0,sa,sa),byrow=TRUE,nrow=6)
Ntotal <- 4000/sum(eigen.analysis(M)$stable.stage*pb) # begin with 4m breeding birds
N[,1] <- Ntotal*eigen.analysis(M)$stable.stage

# simulate
for (j in 1:z) { # of generations
  ndx <- (j-1)*p+1
  
  for (i in ndx:((ndx-1)+p)) # get stochastic reproduction and survival for p years
  {
    rp <- rP(blocks[[j]],1)
    sa <- phi.fi(blocks[[j]],1)
    
    for (k in 1:p) # impose quota for p years
    {
      M <- matrix(c(0,0, s1*pb3*B*rp[k]/2, s1*pb4*B*rp[k]/2, s1*pb5*B*rp[k]/2, s1*pb6*B*rp[k]/2,
                   sa[k],0,0,0,0,0,
                   0,sa[k],0,0,0,0,
                   0,0,sa[k],0,0,0,
                   0,0,0,sa[k],0,0,
                   0,0,0,0,sa[k],sa[k]),byrow=TRUE,nrow=6)
      pv <- (v*N[,i]) / sum(v*N[,i]) 
      H <- pv*q[j]
      N[,i+1] <- pmax(0, (M%*%N[,i]) - (H))
    } # k loop
  } # i loop
  
  # calculate quota for next period
  next.mean.sst <- mean(blocks[[j+1]])
  NB.end <- (N[3,ndx+p]*pb3 + N[4,ndx+p]*pb4 + N[5,ndx+p]*pb5 + N[6,ndx+p]*pb6)
  q[j+1] <- snap_quota(NB.end/2,next.mean.sst,Nvals,Tvals,Q)
  
} # j loop

summary(q)
sum(q==0)/(z+1)
q.pos <- q[q>0]; summary(q.pos)
r <- rle(q)
zeros <- r$lengths[r$values==0]
summary(zeros)
Nbp <- NULL
for (i in 1:(t+1)) Nbp[i] <- (N[3,i]*pb3 + N[4,i]*pb4 + N[5,i]*pb5 + N[6,i]*pb6)/2
plot(Nbp, type='l',lwd=2); abline(h=NB.min/2)
summary(Nbp)
sum(Nbp>=2000)/(t+1)
