#############################################################
## Financial Econometrics 2ST119
## Problem Set 1, solutions
## author: Yukai Yang
## Department of Statistics, Uppsala University
#############################################################

library(tidyverse)
library(FE)


#############################################################
## Part 1: Statistical Properties of Asset Returns
#############################################################

## A. Distributional properties of Dow Jones index returns

# 1.

# plot the log returns
plot(DJ_d$r_Dow_Jones,type='l', main='log returns of the daily Dow Jones index',
     ylab='log return',xlab='time horizon')

ggplot(DJ_d) + geom_line(aes(y=r_Dow_Jones,x=1:nrow(DJ_d)))

DJ_d %>% ggplot() + geom_line(aes(y=r_Dow_Jones,x=1:nrow(DJ_d))) +
  labs(x='time horizon',y='log return')

plot(DJ_w$r_close,type='l', main='log returns of the weekly Dow Jones index',
     ylab='log return',xlab='time horizon')

# descriptive statistics
library(psych)
describe(DJ_d$r_Dow_Jones)
describe(DJ_w$r_close)

help(describe)

# 2.

# 
df = 4
qqnorm(DJ_d$r_Dow_Jones)
qqplot(rt(length(DJ_d$r_Dow_Jones),df=5),DJ_d$r_Dow_Jones)

qqnorm(DJ_w$r_close)
qqplot(rt(length(DJ_w$r_close),df=5),DJ_d$r_Dow_Jones)


# 3.

vdata = DJ_d$r_Dow_Jones
vdata = DJ_w$r_close
vdata = (vdata - mean(vdata))/sd(vdata)


ik = 20
grids = 1:ik/ik

# (a)
vq = pnorm(vdata); hist(vq)
vn = NULL; for(val in grids) vn = c(vn,sum(vq <= val))
vn = c(vn[1],diff(vn))
test = sum((vn-length(vdata)/ik)**2/(length(vdata)/ik))
cat("test =",test," df =",ik-3," p-value =",1-pchisq(test,df=ik-3))

# (b)
df = 3
ndata = vdata*sqrt(df/(df-2))

vq = pt(ndata,df=5); hist(vq)
vn = NULL; for(val in grids) vn = c(vn,sum(vq <= val))
vn = c(vn[1],diff(vn))
test = sum((vn-length(vdata)/ik)**2/(length(vdata)/ik))
cat("test =",test," df =",ik-3," p-value =",1-pchisq(test,df=ik-3))

# (c)
alpha = 0.1514736
sigma = 4.0013995
alpha = 0.33
sigma = 2.73
ndata = vdata*sqrt((1-alpha) + alpha*sigma**2)

vq = (1-alpha)*pnorm(ndata) + alpha*pnorm(ndata,sd=sigma); hist(vq)
vn = NULL; for(val in grids) vn = c(vn,sum(vq <= val))
vn = c(vn[1],diff(vn))
test = sum((vn-length(vdata)/ik)**2/(length(vdata)/ik))
cat("test =",test," df =",ik-3," p-value =",1-pchisq(test,df=ik-3))

func <- function(vx){
  alpha = vx[1]
  sigma = vx[2]
  ndata = vdata*sqrt((1-alpha) + alpha*sigma**2)
  
  vq = (1-alpha)*pnorm(ndata) + alpha*pnorm(ndata,sd=sigma)
  vn = NULL; for(val in grids) vn = c(vn,sum(vq <= val))
  vn = c(vn[1],diff(vn))
  return(sum((vn-length(vdata)/ik)**2/(length(vdata)/ik)))
}

func(c(0.15,4))

optim(par=c(0.1,4),fn=func,method="BFGS")


#############################################################

## B. Dynamical properties of financial return series

# 1.

lret = apply(log(index_d),2,diff)
summary(lret)
matplot(lret,type='l',ylab='log returns',xlab='time horizon')

sum(is.na(lret[,'FRCAC40']))


tmp = 6
colnames(lret)[tmp]; acf(lret[!is.na(lret[,tmp]),tmp])


# 2.

LB <- function(vx,lag,ip){
  tmp = acf(vx,lag.max=lag,plot=F)$acf
  tmp = tmp[2:(lag+1)]**2
  test = sum(tmp/(length(vx)-1:lag))*length(vx)*(length(vx)+2)
  return(list(test=test, pval=1-pchisq(test,df=lag-ip)))
}

tmp = 6
LB(vx=lret[!is.na(lret[,tmp]),tmp],lag=100,ip=0)


# 3.

# choose a pair
pair = c(3,7)

# cross acf pacf
mx = lret[394:dim(lret)[1],pair]
mx = lret[,pair]
cor(mx,use="complete.obs")
cor(mx[2:dim(mx)[1],], mx[1:(dim(mx)[1]-1),],use="complete.obs")


acf(lret[394:dim(lret)[1],pair])
pacf(lret[394:dim(lret)[1],pair])


# 4.
lret2 = lret**2

tmp = 7; colnames(lret2)[tmp]; acf(lret2[!is.na(lret2[,tmp]),tmp])


#############################################################
## Part 2: Asset Return Predictability and Market Efficiency
#############################################################

## A. Testing for asset return predictability

VDR <- function(vr,iq){
  iTT = length(vr)
  im = floor(iTT/iq)
  iT = im*iq
  
  rr = vr[1:iT]
  mu = mean(rr)
  sa2 = var(rr)
  
  arr = NULL
  for(iter in 1:(iT-iq+1))
    arr = c(arr,sum(rr[iter:(iter+iq-1)]))
  
  sc2 = sum((arr-mu*iq)**2)/iq/(iT-iq+1)/(1-(1/im))
  
  VD = sc2 - sa2
  VR = sc2/sa2
  tmp = sqrt(2*(2*iq-1)*(iq-1)/3/iq)
  
  VD = VD*sqrt(iT)/sa2/tmp
  VR = (VR-1)*sqrt(iT)/tmp
  
  return(list(VD=VD, VR=VR))
}

# check that sd and var is computing unbiased sigma
tmp = rnorm(100)
sd(tmp)**2
var(tmp)
sum((tmp-mean(tmp))**2)/99


VDR(vr=DJ_d$r_Dow_Jones,iq=5)
VDR(vr=DJ_w$r_close,iq=5)


#############################################################

## B. Testing for Return Predictability in Size Portfolios


