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

#-----
#########################################################################
# PRELIMINARIES
#########################################################################
# sea temperatures
sst = read.csv('Sep-Aug SST.csv',header = F)
colnames(sst)=c('year','sst')
sst = sst[-nrow(sst),]
acf(sst$sst,lag.max=80)
data <- ts(sst$sst)
full_model = auto.arima(data,stepwise=FALSE)
summary(full_model) # non-stationary

sst.sub = subset(sst,year>=1920) # eliminate early very cold period
data <- ts(sst.sub$sst)
best_model = auto.arima(data,stepwise=FALSE) # now is stationary
summary(best_model)
coefs <- best_model$coef
sd.best_model = sqrt(best_model$sigma2)
mean.sst = mean(sst.sub$sst)
summary(sst.sub$sst)
res = 300
dev.new(height=4*res,width=6*res,noRStudioGD = TRUE)
par(mar=c(5,6,2,2))
hist(sst.sub$sst,breaks=seq(6,9,0.25),main='',xlab='Sea surface temperature',las=1,cex.lab=1.5); box()
abline(v=c(6.1,7.8),lty=2,lwd=2) # temps of positive population growth
abline(v=6.966655,lty=3,lwd=2)

dev.new(height=4*res,width=6*res,noRStudioGD = TRUE)
par(mar=c(5,6,2,2))
loess_model <- loess(sst ~ year, data = sst,span=0.25)
years=1878+(0:145)
plot(years,sst$sst, type='l',las=1,ylim=c(5,9),
     xlab='Year',ylab='SST',xaxp=c(1875,2025,15),lwd=1,cex.lab=1.5)
abline(h=5:9,v=seq(1875,2025,10),lty=2,col='lightgray')
lines(sst$year,predict(loess_model),lwd=3,col='black')
lines(1920:2023,best_model$fitted,type='l',col='blue',lwd=3)
dev.off()


# relative production function based on Hansen et al 2021 in the Westmans
Praw = read_xlsx('gcb15665-sup-0002-datas1.xlsx',na='NA',sheet='Data',col_names=TRUE) 
chkndx = subset(Praw, Year>=1878 & Year<=2005)
ppi = 300
dev.new(height=6*ppi,width=9*ppi,noRStudioGD = TRUE)
 par(mar=c(5,6,2,2))
 with(chkndx,plot(SST,(Prod/6.6),xlim=c(6,8.5),ylim=c(0,1.5),ylab='Relative production',
                 xlab='SST',las=1,cex.lab=1.5)) # normalize by max modeled Prod = 6.6 from Erpur 2021

#... Fit normal kernel 
  prod.norm = chkndx$Prod/6.6
  SST.prod = chkndx$SST
  fit <- nls( prod.norm ~ exp(-0.5 * ((SST.prod - Topt)/Tsd)^2),
           start = list( Topt = 7, Tsd = 1))
  summary(fit)
  newSST = data.frame(SST.prod = seq(6,8.6,0.1))
  yfit <- predict(fit, newdata = newSST)
  lines(newSST$SST.prod, yfit, col = "blue", lwd = 3)
dev.off()

  Topt = fit$m$getPars()[1]
  Tsd =  fit$m$getPars()[2]
  rP = function(temp) exp(-0.5 * ((temp-Topt)/Tsd)^2) # relative production fcn from fitted kernel
  
# production data to get maximum Iceland productivity
  Pobs = read.csv('ProdMonitor.csv',header=TRUE)
  
  # Split the data by group
  split_data <- split(Pobs, Pobs$Location)
  
  # Set up the plot
  res = 300
  dev.new(height=4*res,width=6*res,noRStudioGD = TRUE)
  par(mar=c(5,6,2,2))
  plot(Pobs$Year, Pobs$P, type = "n", xlab = "Year",
       ylab = "Chick production (mean-centered)",cex.lab=1.5,
       las=1,ylim=c(-0.5,0.5))
  grid()
  # Add lines for each group
  color  <- rep(c("#000000", "#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),2)
  ndx = c(4,5,7,11)
  # ndx=1:12
  for (i in 1:12) {
    if (i %in% ndx) lines(split_data[[i]]$Year, split_data[[i]]$P
                                                  - mean(split_data[[i]]$P,na.rm=TRUE)
                          , col = color[i],
                          type='b',pch=16,lwd=3)
  }
  cnames = names(split_data)[ndx]; cnames[1]='Elli?aey'
  legend("topleft", legend = cnames, col = color[ndx],
         lty = 1,pch=16,bty='n')
  abline(h=0,lwd=2,lty=1)

# use burrows for weighting  
  holes.loc = data.frame(unique(Pobs$Location))
  holes.loc$burrow = c(177.0,10.1,45.2,31.0,99.9,36.5,38.4,117.0,NA,425.0,
                       1125.0,NA)
  colnames(holes.loc) = c('loc','burrow')
  holes.loc.sub = subset(holes.loc,burrow>0)
  Pobs.sub <- Pobs[Pobs$Location %in% holes.loc.sub$loc, ]
  max.pobs = NULL
  max.pobs[1] = max(subset(Pobs.sub,Location=='Papey')$P,na.rm=TRUE)
  max.pobs[2] = max(subset(Pobs.sub,Location=='Hholmi')$P,na.rm=TRUE)
  max.pobs[3] = max(subset(Pobs.sub,Location=='Drangey')$P,na.rm=TRUE)
  max.pobs[4] = max(subset(Pobs.sub,Location=='G_Stgrf')$P,na.rm=TRUE)
  max.pobs[5] = max(subset(Pobs.sub,Location=='Grimsey')$P,na.rm=TRUE)
  max.pobs[6] = max(subset(Pobs.sub,Location=='Lundey')$P,na.rm=TRUE)
  max.pobs[7] = max(subset(Pobs.sub,Location=='Vigur')$P,na.rm=TRUE)
  max.pobs[8] = max(subset(Pobs.sub,Location=='Akurey')$P,na.rm=TRUE)
  max.pobs[9] = max(subset(Pobs.sub,Location=='Elli?aey')$P,na.rm=TRUE)
  max.pobs[10] = max(subset(Pobs.sub,Location=='Eyjar')$P,na.rm=TRUE)
  (maxP = weighted.mean(max.pobs,holes.loc.sub$burrow))
  
  mean.BO = NULL
  #xyears = subset(Pobs.sub,Year>=2011 & Year<=2020)
  xyears = subset(Pobs.sub,Year>=2020 & Year<=2023)
  mean.BO[1] = mean(subset(xyears,Location=='Papey')$BO,na.rm=TRUE)
  mean.BO[2] = mean(subset(xyears,Location=='Hholmi')$BO,na.rm=TRUE)
  mean.BO[3] = mean(subset(xyears,Location=='Drangey')$BO,na.rm=TRUE)
  mean.BO[4] = mean(subset(xyears,Location=='G_Stgrf')$BO,na.rm=TRUE)
  mean.BO[5] = mean(subset(xyears,Location=='Gr?msey')$BO,na.rm=TRUE)
  mean.BO[6] = mean(subset(xyears,Location=='Lundey')$BO,na.rm=TRUE)
  mean.BO[7] = mean(subset(xyears,Location=='Vigur')$BO,na.rm=TRUE)
  mean.BO[8] = mean(subset(xyears,Location=='Akurey')$BO,na.rm=TRUE)
  mean.BO[9] = mean(subset(xyears,Location=='Elli?aey')$BO,na.rm=TRUE)
  mean.BO[10] = mean(subset(Pobs.sub,Location=='Eyjar')$BO,na.rm=TRUE)
  Nb_percolony = (mean.BO*holes.loc.sub$burrow)*2
  BP = sum(Nb_percolony)
  
  # check predicted production against 2010-2024 burrow monitoring data: Puffin Rally 2010-2024.xlsx
  B = maxP
  Pobs = numeric(15)
  for (i in 1:15) {
    Pobs[i] = weighted.mean(Pobs.sub[Pobs.sub$Year==i+2009,7],holes.loc.sub$burrow,na.rm=TRUE)
  } 
  s = c(sst.sub[sst.sub$year>=2010,2],NA)
  Ppred = rP(s)*B
  res = 300
  dev.new(height=4*res,width=6*res,noRStudioGD = TRUE)
  par(mar=c(5,6,2,2))
  plot(Pobs,Ppred,xlim=c(0,0.7),ylim=c(0,0.7),las=1,pch=16,cex.lab=1.5);grid();abline(0,1)
  text(Pobs+0.02,Ppred,2010:2024,cex=1.25)
  
# demography
#sa = 0.92 # natural survival ages 1-6+ from Erpur 2024 status report
#sa = 0.874 # from Hafnarh?lmi with no hunting Erpur email 9//11/25
# survival from Halfdan thesis N?tt?rufr Lei?b Til_ritryna_2014
  sa = c(.975,.924,.935,.805,.93,.935,.935,.902,.882,.943,.971,.886,.843,.804,.935,.946,.868,.859,.881,.95,.874)
  sa.mu = mean(sa); sa.sd = sd(sa)
  sa.parm = fitdist(sa,'beta',method='mle')$estimate
  beta.par = MOM.beta(0.874,0.031^2) # Erpur from Hafnarh?lmi 
  sa.nodes = gauss.quad.prob(5,dist='beta',alpha=sa.parm[1],beta=sa.parm[2])

  #  curve(dbeta(x,beta.par$a,beta.par$b), from=0.5, to=1.0,lwd=3,col='blue')
  res = 300
  dev.new(height=4*res,width=6*res,noRStudioGD = TRUE)
  par(mar=c(5,6,2,2))
  hist(sa,breaks=seq(0.5,1,0.05),main='',xlab="Adult survival",freq=TRUE,cex.lab=1.5,xlim=c(0.7,1.0),las=1) # ,add=TRUE) #,col=NULL); 
  box()
  curve(dbeta(x,sa.parm[1],sa.parm[2]),col='darkred',lwd=3,add=TRUE)
  abline(v=mean(sa),lwd=2,lty=2);text(0.88,1,'mean \n(n=21)',cex=1.5)

sa = sa.mu
s = c(0.5, 0.6, 0.7) # natural survival age 0 to 1 (uncertain)
s1 = s[3] # as in Carl's work, need this to better match expected growth in absence of harvest
pb3 = pb4 = 0.067; pb5 = 0.75; pb6 = 1 # probability of breeding from Erpur
B =  0.729608 # max reproduction per pair (calculated in "puffin harvest mgmt.R")
f = B * rP(mean.sst) / 2 # mean repro by breeding adult

# generation time
#x = 3:50 # reproductive ages (use enough that lx drops to ~0.01)
end = 33
x = 3:end # 3% rule (Johnson et al. 2012) (depends on value of sa used)
#x = 3:36 # max reported in captivity (Painter 2022)
(lx = s1*sa^(x-1)) # survival until age x (1st year = 0.7; sa thereafter)
mx = c(f*pb3,f*pb4,f*pb5,rep(f*pb6,x[length(x)]-5)) # reproductive rate; no senescence
num = sum(x*lx*mx) # Caswell 2001: 128
den = sum(lx*mx)
(G = num/den)

# no-harvest, mean matrix
A = matrix(c(0,0,s1*pb3*B*rP(mean.sst)/2,s1*pb4*B*rP(mean.sst)/2,s1*pb5*B*rP(mean.sst)/2,s1*pb6*B*rP(mean.sst)/2,
             sa,0,0,0,0,0,
             0,sa,0,0,0,0,
             0,0,sa,0,0,0,
             0,0,0,sa,0,0,
             0,0,0,0,sa,sa),byrow=TRUE,nrow=6)
e = eigen(A)
ei = e$values
eA = eigen.analysis(A)

(P2 = 2*pi/atan(Im(ei[2])/Re(ei[2]))) # period of oscillations (Caswell2001:100) 
(tx = log(100)/log(eA$damping.ratio)) # decay time for transient dynamics to 1% of original size
pb = c(0,0,pb3,pb4,pb5,pb6)
age = 1:6
AFB = sum(age*pb)/sum(pb) # age at first breeding
(transient.period = AFB/log(eA$damping.ratio)) # Caswell 2001 approximation
GP = c(G,P2,tx)

# Table SI Hansen 2021 article (proportions in the harvest)
pSI = c(0.002,0.163,0.339,0.206,0.080,0.027+.024+.021+.018+.016+.014+.012)
pSI = c(pSI[1:5],1-sum(pSI[1:5])) # to account for no harvest records of 13+ birds)
v = pSI/eA$stable.stage 
(v = round(v/max(v) ,2)) # age-specific vulnerabilities to harvest relative to max(v)
RVx = eA$repro.value
(RVH = sum(v*RVx)/sum(v))
sum(eA$stable.stage*RVx)

#-----

#########################################################################
# CALCULATE QUOTAS FOR VARIOUS NUMBERs OF BREEDING PAIRS AND MEAN SSTs
#########################################################################
m = 2 # choose G or P2
goal = 4000 # desired breeding individuals
sa = sa.mu
# exp.future = forecast(best_model,h=round(GP[2]))
# plot(exp.future,las=1)
# mu.exp = mean(exp.future$mean)

Nbps = seq(1000,3000,100) # breeding pairs
ssts = seq(6,8.5,0.1)
Hquota = matrix(nrow=length(Nbps),ncol=length(ssts))

for (j in 1:length(Nbps)) {
  for (i in 1:length(ssts)) {
    Nbp = Nbps[j] # current number of breeding pairs (thousands)
    mu.sst = ssts[i] # expected mean over next G or P2 years 
    Z = matrix(c(0,0,s1*pb3*B*rP(mu.sst)/2,s1*pb4*B*rP(mu.sst)/2,s1*pb5*B*rP(mu.sst)/2,s1*pb6*B*rP(mu.sst)/2,
                 sa,0,0,0,0,0,
                 0,sa,0,0,0,0,
                 0,0,sa,0,0,0,
                 0,0,0,sa,0,0,
                 0,0,0,0,sa,sa),byrow=TRUE,nrow=6)
    eZ = eigen.analysis(Z)
    #... convert breeding pairs to population vector
    totN = (Nbp*2)/sum(eZ$stable.stage[3]*pb3,eZ$stable.stage[4]*pb4,eZ$stable.stage[5]*pb5,eZ$stable.stage[6]*pb6)
    init.vec = eZ$stable.stage*totN
    N = matrix(nrow=6,ncol=round(GP[m])+1); N[,1]=init.vec
    for (k in 1:round(GP[m])) { N[,k+1] = Z%*%N[,k] }
    #... allowable harvest
    NB.end = N[3,round(GP[m])+1]*pb3 + N[4,round(GP[m])+1]*pb4 + N[5,round(GP[m])+1]*pb5 + N[6,round(GP[m])+1]*pb6
    NB.min = goal # number of min # of breeding individuals in thousands
    loss.max = ifelse((NB.end-NB.min)>0,(NB.end-NB.min),0)
    RVx = eZ$repro.value
    RVH = sum(v*RVx)/sum(v) # mean RV of harvested bird, weighted by vulnerability
    Harvest.total = loss.max/RVH # allowable harvest in one generation
    Hquota[j,i] = Harvest.total/round(GP[m]) # annual allowable
    
  }
}

# plot policy
n=Nbps
s=ssts
z=matrix(Hquota,nrow=length(s),ncol=length(n),byrow=T)
#zz = interp(opt[,1],opt[,2],opt[,3])

ppi = 300
dev.new(height=6*ppi,width=9*ppi,noRStudioGD = TRUE)
par(mar=c(5,5,5,5))
#clr = mako(18)
max(Hquota)
hi = 525; inc =25
Hseq = seq(0,hi,inc)
clr = hcl.colors((hi/inc)+1,"Spectral")
filled.contour(x=n,y=s,z=t(z),col=clr,levels=Hseq,ylab='Expected mean SST',main='Min # breeding pairs = 2.0m; Periodicity P2 = 7 years',
               xlab='Current number of breeding pairs (thousands)',cex.lab=1.25,
               key.title=title(main='Annual quota\n(thousands)',cex.main=0.75),
               plot.axes = { axis(1, Nbps) 
                 axis(2, seq(6,8.5,0.25)) })

#-----

#############################################
# lambda for a range of sst and sa
#############################################
temp = seq(5.5,8.5,0.1)
san = seq(0.85,0.95,0.01)
lam = matrix(nrow=length(temp),ncol=length(san))
for (i in 1:length(temp)) {
  for (j in 1:length(san)) {
    Z = matrix(c(0,0,s1*pb3*B*rP(temp[i])/2,s1*pb4*B*rP(temp[i])/2,s1*pb5*B*rP(temp[i])/2,s1*pb6*B*rP(temp[i])/2,
                 san[j],0,0,0,0,0,
                 0,san[j],0,0,0,0,
                 0,0,san[j],0,0,0,
                 0,0,0,san[j],0,0,
                 0,0,0,0,san[j],san[j]),byrow=TRUE,nrow=6)
    lam[i,j] = eigen.analysis(Z)$lambda
  }    
}
dev.new(height=6*ppi,width=9*ppi,noRStudioGD = TRUE)
par(mar=c(5,6,2,2))
contour(x=san,y=temp,z=t(lam),las=1,xlab='Adult survival',ylab='Sea surface temperature',lwd=3,labcex=1.5,cex.lab=1.5)
abline(v=sa.mu,h=Topt,lty=2,lwd=2)
grid()
dev.off()

x=cbind(temp,lam)
rbind(c(NA,san),x)

#-----

# asymptotic properties of constant h rate
set.seed((123))
size = 10000
sst.sim = arima.sim(n = size, list(ar = coefs[grep("ar", names(coefs))],
                                   ma = coefs[grep("ma", names(coefs))]), 
                    sd = sd.best_model) + coefs[3]
h = seq(0,0.05,0.005)
lam = matrix(nrow=size,ncol=length(h))
for (j in 1:length(h))
{  
  for (i in 1:size) 
  {
    sa.new = rbeta(1,sa.parm[1],sa.parm[2])
    sa = sa.new*(1-v*h[j])
    temp = sst.sim[i]
    H = matrix(c(0,0,s1*pb3*B*rP(temp)/2,s1*pb4*B*rP(temp)/2,s1*pb5*B*rP(temp)/2,s1*pb6*B*rP(temp)/2,
                 sa[1],0,0,0,0,0,
                 0,sa[2],0,0,0,0,
                 0,0,sa[3],0,0,0,
                 0,0,0,sa[4],0,0,
                 0,0,0,0,sa[5],sa[6]),byrow=TRUE,nrow=6)
    lam[i,j] = eigen.analysis(H)$lambda1
  }
}  
#summary(lam)
#apply(lam,2,median)
neg.lam.prob = NULL
for (i in 1:ncol(lam))
{ neg.lam.prob[i] = sum(lam[,i]<1)/size }
cbind(h,neg.lam.prob)

#-----
#########################################################################
# SIMULATION TESTING
#########################################################################
# no-harvest, mean matrix
A = matrix(c(0,0,s1*pb3*B*rP(mean.sst)/2,s1*pb4*B*rP(mean.sst)/2,s1*pb5*B*rP(mean.sst)/2,s1*pb6*B*rP(mean.sst)/2,
             sa,0,0,0,0,0,
             0,sa,0,0,0,0,
             0,0,sa,0,0,0,
             0,0,0,sa,0,0,
             0,0,0,0,sa,sa),byrow=TRUE,nrow=6)
eA = eigen.analysis(A)

m = 2 # choose G or P2 or tx
p = GP[m] # periodicity
z = 5000 # how many generations (each generation gets a random initial vector N and random sa)
t = round(p) * z
initBP = 1750 # initial breeding PAIRS
NB.min = 3500 # minimum number of breeding BIRDS (Carl 1950-2000 mean=4430; a low to hi)
set.seed(123)

#... simulate long time series of sst
# this code only works for a stationary time series
sst.sim = arima.sim(n = t, list(ar = coefs[grep("ar", names(coefs))],
                                ma = coefs[grep("ma", names(coefs))]), 
                                sd = sd.best_model) + coefs[3]

blocks <- split(sst.sim, rep(1:(length(sst.sim)/round(p)), each = round(p))) # get SST in blocks of G years
# plot(sst.sim)
# time=1:t
# loess_model <- loess(sst.sim ~ time,span=0.20)
# plot(sst.sim,type='l',col='blue',lwd=0.5,las=1)
# abline(h=mean(sst.sub$sst),lty=2)
# grid()
#  maG=NULL # G-year running averages
#  for (i in round(G):t) maG[(i-(round(G)-1))] = mean(sst.sim[(i-(round(G)-1)):i])
#  lines(round(G):t,maG,lwd=3)
# lines(predict(loess_model),lwd=2,col='red')


Hquota = NBP.stop = NULL
AgeD.stop = AgeD.init = matrix(nrow=6,ncol=z)
pv = matrix(nrow=6,ncol=(round(p)))

for (j in 1:z) {
  #... calculate annual harvest quota
  next.mean.sst = mean(blocks[[j]])
  sa = sa.mu
  Z = matrix(c(0,0,s1*pb3*B*rP(next.mean.sst)/2,s1*pb4*B*rP(next.mean.sst)/2,s1*pb5*B*rP(next.mean.sst)/2,s1*pb6*B*rP(next.mean.sst)/2,
               sa,0,0,0,0,0,
               0,sa,0,0,0,0,
               0,0,sa,0,0,0,
               0,0,0,sa,0,0,
               0,0,0,0,sa,sa),byrow=TRUE,nrow=6)
  eZ = eigen.analysis(Z)

  totN = (initBP*2)/sum(eZ$stable.stage[3]*pb3,eZ$stable.stage[4]*pb4,eZ$stable.stage[5]*pb5,eZ$stable.stage[6]*pb6)
  init.vec = eZ$stable.stage*totN
  N = matrix(nrow=6,ncol=round(GP[m])+1); N[,1]=init.vec
  for (i in 1:round(p)) {
    N[,i+1] = Z%*%N[,i]
  }
  NB.end = (N[3,round(p)+1]*pb3 + N[4,round(p)+1]*pb4 + N[5,round(p)+1]*pb5 + N[6,round(p)+1]*pb6)
  loss.max = ifelse((NB.end-NB.min)>0,(NB.end-NB.min),0)
  RVx = eZ$repro.value
  RVH = sum(v*RVx)/sum(v) # mean RV of harvested bird, weighted by vulnerability
  Harvest.total = loss.max/RVH
  Hquota[j] = Harvest.total/round(p)

  # simulate dynamics
  
  # perturbation <- rnorm(length(eZ$stable.stage), mean = 0, sd = 0.05)
  # sad_perturbed <- eZ$stable.stage + perturbation
  # sad_perturbed[sad_perturbed < 0] <- 0  # truncate negatives
  # SAD.init <- sad_perturbed / sum(sad_perturbed) # renormalize
  perturb = rdirichlet(1,eZ$stable.stage*10)
  totN = (initBP*2)/sum(perturb[3]*pb3,perturb[4]*pb4,perturb[5]*pb5,perturb[6]*pb6)
  init.vec = perturb*totN
  N = matrix(nrow=6,ncol=round(p)+1); N[,1]=init.vec
  H = matrix(nrow=6,ncol=round(p))
  rp = rP(blocks[[j]])
  
  for (i in 1:round(p)) {
#    sa = sa.mu
    sa = rbeta(1,sa.parm[1],sa.parm[2])
    M = matrix(c(0,0, s1*pb3*B*rp[i]/2, s1*pb4*B*rp[i]/2, s1*pb5*B*rp[i]/2, s1*pb6*B*rp[i]/2,
                 sa,0,0,0,0,0,
                 0,sa,0,0,0,0,
                 0,0,sa,0,0,0,
                 0,0,0,sa,0,0,
                 0,0,0,0,sa,sa),byrow=TRUE,nrow=6)
    pv[,i] = (v*N[,i]) / sum(v*N[,i]) 
    H[,i] = pv[,i]*Hquota[j]
    N[,i+1] = pmax(0, (M%*%N[,i]) - H[,i])
  }
 
  NBP.stop[j] = (N[3,round(p)+1]*pb3 + N[4,round(p)+1]*pb4 + N[5,round(p)+1]*pb5 + N[6,round(p)+1]*pb6)/2
  AgeD.init[,j] = N[,1] / sum(N[,1])
  AgeD.stop[,j] = N[,round(p)+1] / sum(N[,round(p)+1])
}

 # quotas = cbind(as.numeric(lapply(blocks,mean)),Hquota)
 # plot(quotas[,1],quotas[,2])

#NB.stop final # of breeding individuals
hist(NBP.stop);abline(v=NB.min/2,lwd=2,col='darkred') 
mean(NBP.stop); 
sd(NBP.stop)/mean(NBP.stop);
summary(NBP.stop)
sum(NBP.stop>=NB.min/2)/z # frequency of planning units in which the end goal is met 

#Hquota # annual quota per generation
hist(Hquota);
mean(Hquota); 
sd(Hquota)/mean(Hquota); 
summary(Hquota); 
sum(Hquota==0)/z # frequency of closed seasons

#-----
###############################################################################################
# HOW LONG TO GOAL WITH NO HARVEST
###############################################################################################
z = 200 # how many years to inspect
y = 5000 # how many trials of z years
initBP = 1750 # initial breeding PAIRS
goal = 4000 #initBP*2
howmany = numeric(y)
NBbirds = NULL
set.seed = 123

for (j in 1:y) {
  perturb = rdirichlet(1,eA$stable.stage*10)
  totN = (initBP*2)/sum(perturb[3]*pb3,perturb[4]*pb4,perturb[5]*pb5,perturb[6]*pb6)
  init.vec = perturb*totN
  N = matrix(nrow=6,ncol=z+1); N[,1]=init.vec
  sst.sim = arima.sim(n = z, list(ar = coefs[grep("ar", names(coefs))],
                                  ma = coefs[grep("ma", names(coefs))]), 
                                  sd = sd.best_model) + coefs[3]
  
  for (i in 1:z) {
   sa = rbeta(1,sa.parm[1],sa.parm[2])
   M = matrix(c(0,0, s1*pb3*B*rP(sst.sim[i])/2, s1*pb4*B*rP(sst.sim[i])/2, s1*pb5*B*rP(sst.sim[i])/2, s1*pb6*B*rP(sst.sim[i])/2,
               sa,0,0,0,0,0,
               0,sa,0,0,0,0,
               0,0,sa,0,0,0,
               0,0,0,sa,0,0,
               0,0,0,0,sa,sa),byrow=TRUE,nrow=6)
   N[,i+1] = pmax(0, (M%*%N[,i]))

   NBbirds[i] =  (N[3,i]*pb3 + N[4,i]*pb4 + N[5,i]*pb5 + N[6,i]*pb6) 
 }
 
 howmany[j] = which(NBbirds>=goal)[1]
} 

summary(howmany)
hist(howmany,main='',xlab="Years",las=1)

#-----

#################################################################################################
# simulate P2=7 for a long-time series, retaining N from previous ending pop as new initial pop
#################################################################################################
m = 2 # choose G or P2 or tx
p = round(GP[m]) # periodicity

set.seed(123)
n = 10000
sst.sim = arima.sim(n = n*p, list(ar = coefs[grep("ar", names(coefs))],
                                ma = coefs[grep("ma", names(coefs))]), 
                    sd = sd.best_model) + coefs[3]

blocks <- split(sst.sim, rep(1:(length(sst.sim)/round(p)), each = p)) # get SST in blocks of G years

sa = sa.mu
B = maxP
A = matrix(c(0,0,s1*pb3*B*rP(mean.sst)/2,s1*pb4*B*rP(mean.sst)/2,s1*pb5*B*rP(mean.sst)/2,s1*pb6*B*rP(mean.sst)/2,
             sa,0,0,0,0,0,
             0,sa,0,0,0,0,
             0,0,sa,0,0,0,
             0,0,0,sa,0,0,
             0,0,0,0,sa,sa),byrow=TRUE,nrow=6)
eA = eigen.analysis(A)

z = n - 5000 # how many "generations" (set to less than n)
#z = 50

t = p * z # total years
RVx = eA$repro.value
RVH = sum(v*RVx)/sum(v) # mean RV of harvested bird, weighted by vulnerability
initBP = 1750 # initial breeding PAIRS
NB.min = 3500 # minimum number of breeding BIRDS (Carl 1950-2000 mean=4430; a low to hi)

N = matrix(nrow=6,ncol=(t+1))
s = matrix(nrow=6,ncol=p+1)
Hquota = NULL
totN = (initBP*2)/sum(eA$stable.stage[3]*pb3,eA$stable.stage[4]*pb4,eA$stable.stage[5]*pb5,eA$stable.stage[6]*pb6)
N[,1] = eZ$stable.stage*totN
s[,1] = N[,1]

# calculate initial quota
next.mean.sst = mean(blocks[[1]])
A = matrix(c(0,0,s1*pb3*B*rP(next.mean.sst)/2,s1*pb4*B*rP(next.mean.sst)/2,s1*pb5*B*rP(next.mean.sst)/2,s1*pb6*B*rP(next.mean.sst)/2,
             sa,0,0,0,0,0,
             0,sa,0,0,0,0,
             0,0,sa,0,0,0,
             0,0,0,sa,0,0,
             0,0,0,0,sa,sa),byrow=TRUE,nrow=6)

 for (i in 1:p) { s[,i+1] = A%*%s[,i] }  
  NB.end = (s[3,(p+1)]*pb3 + s[4,(p+1)]*pb4 + s[5,(p+1)]*pb5 + s[6,(p+1)]*pb6)
  loss.max = ifelse((NB.end-NB.min)>0,(NB.end-NB.min),0)
  Harvest.total = loss.max/RVH
  Hquota[1] = Harvest.total/p

# quota = 36 
# Hquota[1] = quota
 
# simulate
for (j in 1:z) {
    ndx = (j-1)*7+1

   # impose quota for p years
   for (i in ndx:((ndx-1)+p))
    { 
#    sa = rbeta(1,sa.parm[1],sa.parm[2]) # set random sa for block
#    sa = sa.mu                          # set mean sa for block
    rp = rP(blocks[[j]])
    for (k in 1:7) 
      {
      sa = rbeta(1,sa.parm[1],sa.parm[2])   # draw random sa for each year of block
      M = matrix(c(0,0, s1*pb3*B*rp[k]/2, s1*pb4*B*rp[k]/2, s1*pb5*B*rp[k]/2, s1*pb6*B*rp[k]/2,
                    sa,0,0,0,0,0,
                    0,sa,0,0,0,0,
                    0,0,sa,0,0,0,
                    0,0,0,sa,0,0,
                    0,0,0,0,sa,sa),byrow=TRUE,nrow=6)
       pv = (v*N[,i]) / sum(v*N[,i]) 
       H = pv*Hquota[j]
       N[,i+1] = pmax(0, (M%*%N[,i]) - (H))
      }
    }

# calculate quota for next 7 years
  s = matrix(nrow=6,ncol=p+1)
  s[,1] = N[,ndx+p]
  next.mean.sst = mean(blocks[[j+1]])
  sa = sa.mu
  X = matrix(c(0,0,s1*pb3*B*rP(next.mean.sst)/2,s1*pb4*B*rP(next.mean.sst)/2,s1*pb5*B*rP(next.mean.sst)/2,s1*pb6*B*rP(next.mean.sst)/2,
             sa,0,0,0,0,0,
             0,sa,0,0,0,0,
             0,0,sa,0,0,0,
             0,0,0,sa,0,0,
             0,0,0,0,sa,sa),byrow=TRUE,nrow=6)
  for (i in 1:p) { s[,i+1] = X%*%s[,i] }
  NB.end = (s[3,(p+1)]*pb3 + s[4,(p+1)]*pb4 + s[5,(p+1)]*pb5 + s[6,(p+1)]*pb6)
  loss.max = ifelse((NB.end-NB.min)>0,(NB.end-NB.min),0)
  Harvest.total = loss.max/RVH
  Hquota[j+1] = Harvest.total/p
  
#  Hquota[j+1] = quota
}

summary(Hquota)
sum(Hquota==0)/z
Hquota.pos = Hquota[Hquota>0]; summary(Hquota.pos)
r = rle(Hquota)
zeros = r$lengths[r$values==0]
summary(zeros)
Nbp = NULL
for (i in 1:(t+1)) Nbp[i] = (N[3,i]*pb3 + N[4,i]*pb4 + N[5,i]*pb5 + N[6,i]*pb6)/2
plot(Nbp, type='l',lwd=2); abline(h=NB.min/2)
summary(Nbp)
sum(Nbp>=1750)/(t+1)

#-----

#################################################################################################
########################### PARKING LOT #########################################################
#################################################################################################

#########################################################################
# CALCULATE QUOTAS FOR ANY SPECIFIED NUMBER OF BREEDING PAIRS AND MEAN SST
#########################################################################
exp.future = forecast(best_model,h=round(G))
plot(exp.future,las=1)
mu.exp = mean(exp.future$mean)

k = 2 # choose time frame for quotas
goal= 4430 # goal for breeding individuals 4.43 average 1950-2000 Carl's SPR
Nbp = 1500 # current number of breeding pairs (thousands) 
mu.sst = mu.exp # expected mean over next G years (can be any value desired)
Z = matrix(c(0,0,s1*pb3*B*rP(mu.sst)/2,s1*pb4*B*rP(mu.sst)/2,s1*pb5*B*rP(mu.sst)/2,s1*pb6*B*rP(mu.sst)/2,
             sa,0,0,0,0,0,
             0,sa,0,0,0,0,
             0,0,sa,0,0,0,
             0,0,0,sa,0,0,
             0,0,0,0,sa,sa),byrow=TRUE,nrow=6)
eZ = eigen.analysis(Z)
#... convert breeding pairs to population vector
totN = (Nbp*2)/sum(eZ$stable.stage[3]*pb3,eZ$stable.stage[4]*pb4,eZ$stable.stage[5]*pb5,eZ$stable.stage[6]*pb6)
init.vec = eZ$stable.stage*totN
N = matrix(nrow=6,ncol=round(GP[k])+1); N[,1]=init.vec
for (i in 1:round(GP[k])) {
  N[,i+1] = Z%*%N[,i]
}
#... allowable harvest
NB.end = N[3,round(GP[k])+1]*pb3 + N[4,round(GP[k])+1]*pb4 + N[5,round(GP[k])+1]*pb5 + N[6,round(GP[k])+1]*pb6
NB.min = goal # number of min # of breeding individuals in thousands
(loss.max = ifelse((NB.end-NB.min)>0,(NB.end-NB.min),0))
RVx = eZ$repro.value
(RVH = sum(v*RVx)/sum(v)) # mean RV of harvested bird, weighted by vulnerability
(Harvest.total = loss.max/RVH) # allowable harvest in one generation
(Harvest.total.ann = Harvest.total/round(GP[k])) # annual allowable


