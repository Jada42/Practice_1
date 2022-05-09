#################################################################
##################   LIBRARIES AND FUNCTIONS ####################
#################################################################
rm(list=ls())
lagpad <- function(x, k) {
  c(rep(NA, k), x)[1 : length(x)] 
}

library(psych)
library(glmnet)
library(Hmisc)
library(fmsb)
library(DataCombine)
library(qgraph)
library(pROC)
library(reshape2)
library(lmtest)
library(DataCombine)
library(DescTools)
source('beepday2consec.R', encoding = 'UTF-8')
source('lagData.R', encoding = 'UTF-8')

#See external .R files for 'beepday2consec' and 'lagData' functions
#Taken from mgm mvar internal code

#################################################################
##################  DATA SETUP AND CLEANING #####################
#################################################################

##Read in data
data=read.csv('P001-raw.csv',as.is=TRUE)

###View(data)

names<-scan("name.txt",what = "character", sep = "\n")
colnames(data)<-names

## Duplicate and lag time
lagpad <- function(x, k) {
  c(rep(NA, k), x)[1 : length(x)] 
}

data$lag=lagpad(data$start,1)

## Calculate time differences
data$tdif=as.numeric(difftime(strptime(data$start,"%m/%d/%Y %H:%M"),strptime(data$lag,"%m/%d/%Y %H:%M")))

## Replace NA
data$tdif[is.na(data$tdif)]<- 0

## Calculate cumulative sum of numeric elapsed time 
data$cumsumT=cumsum(data$tdif)

##plot data, inspect for outlying values
plot(data$drinks~data$cumsumT,type='o')
quantile(data$drinks,na.rm=T)

##Create dichotomous drinking variables
data$drinkbin=ifelse(data$drinks>0,1,0)
data$drinkH=ifelse(data$drinks>median(data$drinks,na.rm=TRUE),1,0)
length(which(data$drinkbin==1)) #22
length(which(data$drinkbin==0)) #93
length(which(data$drinkH==1)) #22
length(which(data$drinkH==0)) #93

##Trim time series to even 8 obs/day
options(width=90)
data$start
nrow(data)
datx=data[,]
datx$start

##Use for internal missing rows
new <- rep(NA, length(datx))
datx <- InsertRow(datx, NewRow=new, RowNum = 1)
datx <- InsertRow(datx, NewRow=new, RowNum = 2)
datx <- InsertRow(datx, NewRow=new, RowNum = 3)
datx$start

##Code days of the week
datx$day <- rep(1:7, each=8, length.out=nrow(datx))
## 1 = Tuesday
datx$mon=ifelse(datx$day==7,1,0)
datx$tues=ifelse(datx$day==1,1,0)
datx$wed=ifelse(datx$day==2,1,0)
datx$thur=ifelse(datx$day==3,1,0)
datx$fri=ifelse(datx$day==4,1,0)
datx$sat=ifelse(datx$day==5,1,0)
datx$sun=ifelse(datx$day==6,1,0)

##Code pings
datx$ping=seq(0,7,1)
datx$ping1=ifelse(datx$ping==0,1,0)
datx$ping2=ifelse(datx$ping==1,1,0)
datx$ping3=ifelse(datx$ping==2,1,0)
datx$ping4=ifelse(datx$ping==3,1,0)
datx$ping5=ifelse(datx$ping==4,1,0)
datx$ping6=ifelse(datx$ping==5,1,0)
datx$ping7=ifelse(datx$ping==6,1,0)
datx$ping8=ifelse(datx$ping==7,1,0)


##Temporal variables
datx$linear=scale(datx$cumsumT)
datx$quad=datx$linear^2
datx$cub=datx$linear^3
datx$cosT=cos(((2*pi)/24)*datx$cumsumT)
datx$sinT=sin(((2*pi)/24)*datx$cumsumT)
datx$cos2T=cos(((2*pi)/12)*datx$cumsumT)
datx$sin2T=sin(((2*pi)/12)*datx$cumsumT)
datx$cosW=cos(((2*pi)/168)*datx$cumsumT)
datx$sinW=sin(((2*pi)/168)*datx$cumsumT)

##Index consecutive measurements by ping and day
datx$dayvar=rep(1:(nrow(datx)/8),each=8)
datx$beepvar=datx$ping+1

##Create filter to mark cases with missing data
datx$filter=ifelse(!complete.cases(datx[,c(3:18)]),1,0)

## Remove NAs, suppress row names
daty=subset(datx,filter==0)
row.names(daty)<- NULL
##View(daty)

dayvar=daty$dayvar
beepvar=daty$beepvar

#mean-centering all self-reported continuous variables
for (i in 4:18){
daty[,i] <- (daty[,i] - mean(daty[,i]))
}


##Now we create separate lag0 (contemporaneous) and lag1 (lagged) data structures
##Functions 'beepday2consec' and 'lagData' taken from mgm mvar internal code

##Uses beepvar and dayvar info to index consecutive measurements
daty$consec <- beepday2consec(beepvar = beepvar, dayvar = dayvar)

##Object with original and lagged data
data_lagged <- lagData(data = daty, 
                       lags = 1, 
                       consec = daty$consec)

##Original & lagged data
data_response <- data_lagged$data_response
l_data_lags <- data_lagged$l_data_lags
##n_design <- nrow(data_response) #used in other mvar calculations

## delete rows that cannot be predicted
data_response <- data_response[data_lagged$included, ]
l_data_lags <- lapply(l_data_lags, function(x) x[data_lagged$included, ])

##Original dat, with appropriate rows for lagged analysis (L0 = lag zero)
datyL0=as.data.frame(data_response)

##Lagged data (L1 = lag one)
datyL1=as.data.frame(l_data_lags[[1]])

##Change column names of lagged data
colnames(datyL1)<-colnames(daty)

nrow(datx) #136
nrow(daty) #115
nrow(datyL0) #86
nrow(datyL1) #86

##Recode variables as numeric (via character)
datyL0[,c(3:57)] <- sapply(datyL0[,c(3:57)], as.character)
datyL0[,c(3:57)] <- sapply(datyL0[,c(3:57)], as.numeric)
datyL1[,c(3:57)] <- sapply(datyL1[,c(3:57)], as.character)
datyL1[,c(3:57)] <- sapply(datyL1[,c(3:57)], as.numeric)

psych::describe(datyL0)

##Inspect data
##View(datyL0)
##View(datyL1)

## Training/Testing ##
## 50% of the sample size
smp_size <- floor(0.5 * nrow(datyL1))

##Create random sets of 50% each
##set the seed to make your partition reproducible
i.results <- data.frame(i.test.drinks= rep(NA,1000),i.var.retained.drink= rep(NA,1000),i.auc=rep(NA,1000),
                      i.spec = rep(NA,1000),i.sens=rep(NA,1000),i.brier=rep(NA,1000),
                      i.var.retained.crave= rep(NA,1000),i.r2.crave=rep(NA,1000),
                      i.r2.adj.crave=rep(NA,1000),i.var.retained.want= rep(NA,1000),
                      i.r2.want=rep(NA,1000), i.r2.adj.want=rep(NA,1000))

sample.seed <- 101010
set.seed(sample.seed)
train_ind <- sample(seq_len(nrow(datyL1)), size = smp_size)

##Training and testing sets
trainL0 <- datyL0[train_ind, ]
testL0 <- datyL0[-train_ind, ]
trainL1 <- datyL1[train_ind, ]
testL1 <- datyL1[-train_ind, ]

##Visual spot check for drinking events
##Events to be predicted are located in L0 (lag-zero/contemporaneous) data
even1=seq(0,nrow(trainL0)-1,1)
even2=seq(0,nrow(testL0)-1,1)
plot(drinkbin~even1,trainL0,type='o')
plot(drinkbin~even2,testL0,type='o')
#plot(drinkbin~even,trainL1,type='o')
#plot(drinkbin~even,testL1,type='o')

i.test.drinks<- sum(testL0$drinkbin)
i.train.drinks <- sum(trainL0$drinkbin)

######################################
########  Drinking Events  ###########
######################################

##Prediction variables/Prediction matrices (feature space)
##Vars from datyL1 are lagged, vars from datyL0 are contemporaneous with DV
##Psychosocial variables (left side) + time variables (right side)
trainmat = data.matrix(cbind(trainL1[,c(4:18,26)],trainL0[,c(29:35,37:53)]))
testmat = data.matrix(cbind(testL1[,c(4:18,26)],testL0[,c(29:35,37:53)]))

alpha.set <- seq(from=0.5, to=0, by=-.05)
set.seed(8675309)
for (i in alpha.set){
  i.reg1=cv.glmnet(trainmat,trainL0$drinkbin,family='binomial',standardize=T,alpha=i)
  i.coef1=coef(i.reg1, s = "lambda.min")
  if ((sum(i.coef1[2:length(i.coef1)]!=0))>=1){ 
    i.alpha.auc <- i
    break
  } 
}
i.alpha.auc
if(is.na(i.alpha.auc)==T){
  i.var.retained.drink<- NA
  i.alpha.auc <- NA
  i.coef1 <- NA
  i.auc <- NA
  i.spec <- NA
  i.sens <- NA
  i.brier <- NA
}
if(exists("i.alpha.auc")==F){
  i.var.retained.drink<-NA
  i.alpha.auc <- NA
  i.coef1 <- NA
  i.auc <- NA
  i.spec <- NA
  i.sens <- NA
  i.brier <- NA
}

if (is.na(i.alpha.auc)==F){
  i.coef1
  i.var.retained.drink<- sum(i.coef1[2:length(i.coef1)]!=0)
  i.pred1=predict(i.reg1, newx=testmat, s = "lambda.min", type = "link")
  testL0$i.regpred1=as.numeric(i.pred1)
  i.auc<- pROC::auc(pROC::roc(testL0$drinkbin~testL0$i.regpred1))
  i.auc
  
  roc <- pROC::coords((pROC::roc(testL0$drinkbin~testL0$i.regpred1)),"best", ret=c("threshold", "specificity", "sensitivity"), transpose=FALSE)
  ifelse( nrow(roc[2])>1,i.spec <- NA, i.spec <- roc[2])
  ifelse( nrow(roc[2])>1,i.sens <- NA, i.sens <- roc[3])
  
  i.pred.prob1 <- predict(i.reg1, newx=testmat, s = "lambda.min", type = "response")
  i.brier <- mean((i.pred.prob1-as.numeric(as.character(testL0$drinkbin)))^2)
  i.brier
}

##############################
########  Craving  ###########
##############################

plot(craving~even1,trainL0,type='o')
set.seed(8675309)
for (i in alpha.set){
  
  i.reg2=cv.glmnet(trainmat,trainL0$craving,family='gaussian',standardize=T,alpha=i)
  i.coef2=coef(i.reg2, s = "lambda.min")
  if ((sum(i.coef2[2:length(i.coef2)]!=0))>=1){ 
    i.alpha.crave <- i
    break
  } 
}
i.alpha.crave
i.coef2
i.var.retained.crave<- sum(i.coef2[2:length(i.coef2)]!=0)
i.pred2=predict(i.reg2, newx=testmat, s = "lambda.min", type = "link")
testL0$i.regpred2=as.numeric(i.pred2)
i.r2.crave <- (cor.test(testL0$craving,testL0$i.regpred2)[4]$estimate)^2
i.r2.adj.crave <- 1- (((1-i.r2.crave)*(nrow(testmat)-1))/ 
                        (nrow(testmat) -((sum(i.coef2[2:length(i.coef2)]!=0)) -1)))


############################################
########  I Would Like to Drink  ###########
############################################

plot(wantdrink~even1,trainL0,type='o')
set.seed(8675309)
for (i in alpha.set){
  
  i.reg3=cv.glmnet(trainmat,trainL0$wantdrink,family='gaussian',standardize=T,alpha=i)
  i.coef3=coef(i.reg3, s = "lambda.min")
  if ((sum(i.coef3[2:length(i.coef3)]!=0))>=1){ 
    i.alpha.want <- i
    break
  } 
}
i.alpha.want
i.coef3
i.var.retained.want<- sum(i.coef3[2:length(i.coef3)]!=0)
i.pred3=predict(i.reg3, newx=testmat, s = "lambda.min", type = "link")
testL0$i.regpred3=as.numeric(i.pred3)
i.r2.want <- (cor.test(testL0$wantdrink,testL0$i.regpred3)[4]$estimate)^2
i.r2.adj.want <- 1- (((1-i.r2.want)*(nrow(testmat)-1))/ 
                       (nrow(testmat) -((sum(i.coef3[2:length(i.coef3)]!=0)) -1)))

i.results[1,] <- cbind(i.test.drinks,i.var.retained.drink,i.auc,i.spec,i.sens,i.brier,
                       i.var.retained.crave,i.r2.crave, i.r2.adj.crave,i.var.retained.want, 
                       i.r2.want, i.r2.adj.want)


#View(i.results)

################################################
########  Prep for Nomothetic Analysis  ########
################################################
n.trainL0 <- datyL0[train_ind, ]
n.testL0 <- datyL0[-train_ind, ]
n.trainL1 <- datyL1[train_ind, ]
n.testL1 <- datyL1[-train_ind, ]

ID <- rep(1,nrow(n.trainL0))
n.train<- cbind(ID,n.trainL0,n.trainL1)
ID <- rep(1,nrow(n.testL0))
n.test <- cbind(ID,n.testL0,n.testL1)

nrow(n.train) #43
nrow(n.test) #43

write.csv(n.train, "n.trainp001.csv")
write.csv(n.test, "n.testp001.csv")

##note:idiographic datasets were compiled into long format outside of r

############################################
#############Nomothetic analysis############
############################################

###############################
### DATA SETUP AND CLEANING ###
###############################
n.results <- data.frame(n.var.retained.drink= rep(NA,1000),
                        n.auc=rep(NA,1000),n.spec = rep(NA,1000),n.sens=rep(NA,1000),
                        n.brier=rep(NA,1000),n.var.retained.crave= rep(NA,1000),
                        n.r2.crave=rep(NA,1000),n.r2.adj.crave=rep(NA,1000),
                        n.var.retained.want= rep(NA,1000),n.r2.want=rep(NA,1000), 
                        n.r2.adj.want=rep(NA,1000))
##Read in data
n.train <- read.csv('n.train.csv',as.is=TRUE)
n.test <- read.csv('n.test.csv',as.is=TRUE)


####################################
###  Nomothetic Drinking Events  ###
####################################

##Prediction variables/Prediction matrices (feature space)
##Vars from n.datyL1 are lagged, vars from n.datyL0 are contemporaneous with DV
##Psychosocial variables (left side) + time variables (right side)
n.trainmat = data.matrix(cbind(n.train[,c(62:76,84)],n.train[,c(30:36,38:54)]))
n.testmat = data.matrix(cbind(n.test[,c(62:76,84)],n.test[,c(30:36,38:54)]))



set.seed(8675309)
for (i in alpha.set){
  n.reg1=cv.glmnet(n.trainmat,n.train$drinkbin,family='binomial',standardize=T,alpha=i)
  n.coef1=coef(n.reg1, s = "lambda.min")
  if ((sum(n.coef1[2:length(n.coef1)]!=0))>=1){ 
    n.alpha.auc <- i
    break
  } 
}
n.alpha.auc
if(is.na(n.alpha.auc)==T){
  n.var.retained.drink <- NA
  n.alpha.auc <- NA
  n.coef1 <- NA
  n.auc <- NA
  n.spec <- NA
  n.sens <- NA
  n.brier <- NA
}
if(exists("n.alpha.auc")==F){
  n.var.retained.drink <- NA
  n.alpha.auc <- NA
  n.coef1 <- NA
  n.auc <- NA
  n.spec <- NA
  n.sens <- NA
  n.brier <- NA
}

if (is.na(n.alpha.auc)==F) {n.coef1
  n.var.retained.drink<- sum(n.coef1[2:length(n.coef1)]!=0)
  n.pred1=predict(n.reg1, newx=n.testmat, s = "lambda.min", type = "link")
  n.test$regpred1=as.numeric(n.pred1)
  n.auc<- pROC::auc(pROC::roc(n.test$drinkbin~n.test$regpred1))
  n.auc
  
  roc <- pROC::coords((pROC::roc(n.test$drinkbin~n.test$regpred1)),"best", ret=c("threshold", "specificity", "sensitivity"), transpose=FALSE)
  ifelse( nrow(roc[2])>1,n.spec <- NA, n.spec <- roc[2])
  ifelse( nrow(roc[2])>1,n.sens <- NA, n.sens <- roc[3])
  
  n.pred.prob1 <- predict(n.reg1, newx=n.testmat, s = "lambda.min", type = "response")
  n.brier <- mean((n.pred.prob1-as.numeric(as.character(n.test$drinkbin)))^2)
  n.brier
}

###########################
### Nomothetic Craving  ###
###########################

set.seed(8675309)
for (i in alpha.set){
  
  n.reg2=cv.glmnet(n.trainmat,n.train$craving,family='gaussian',standardize=T,alpha=i)
  n.coef2=coef(n.reg2, s = "lambda.min")
  if ((sum(n.coef2[2:length(n.coef2)]!=0))>=1){ 
    n.alpha.crave <- i
    break
  } 
}
n.alpha.crave
n.coef2
n.var.retained.crave<- sum(n.coef2[2:length(n.coef2)]!=0)
n.pred2=predict(n.reg2, newx=n.testmat, s = "lambda.min", type = "link")
n.test$regpred2=as.numeric(n.pred2)
n.r2.crave <- (cor.test(n.test$craving,n.test$regpred2)[4]$estimate)^2
n.r2.adj.crave <- 1- (((1-n.r2.crave)*(nrow(n.testmat)-1))/ 
                        (nrow(n.testmat) -((sum(n.coef2[2:length(n.coef2)]!=0)) -1)))

##########################################
###  Nomothetic I Would Like to Drink  ###
##########################################

set.seed(8675309)
for (i in alpha.set){
  
  n.reg3=cv.glmnet(n.trainmat,n.train$wantdrink,family='gaussian',standardize=T,alpha=i)
  n.coef3=coef(n.reg3, s = "lambda.min")
  if ((sum(n.coef3[2:length(n.coef3)]!=0))>=1){ 
    n.alpha.want <- i
    break
  } 
}
n.alpha.want
n.coef3
n.var.retained.want<- sum(n.coef3[2:length(n.coef3)]!=0)
n.pred3=predict(n.reg3, newx=n.testmat, s = "lambda.min", type = "link")
n.test$regpred3=as.numeric(n.pred3)
n.r2.want <- (cor.test(n.test$wantdrink,n.test$regpred3)[4]$estimate)^2
n.r2.adj.want <- 1- (((1-n.r2.want)*(nrow(n.testmat)-1))/ 
                       (nrow(n.testmat) -((sum(n.coef3[2:length(n.coef3)]!=0))
                                          -1)))

n.results[1,] <- cbind(n.var.retained.crave,n.auc,n.spec,n.sens,n.brier,
                       n.var.retained.crave,n.r2.crave, n.r2.adj.crave,n.var.retained.want,
                       n.r2.want, n.r2.adj.want)


#View(n.results)



#########################################################
###idiographic application of nomothetic model results###
#########################################################

new.results <- data.frame(new.var.retained.drink= rep(NA,1000),
                          new.auc=rep(NA,1000),new.spec = rep(NA,1000),new.sens=rep(NA,1000),
                          new.brier=rep(NA,1000),new.var.retained.crave= rep(NA,1000),
                          new.r2.crave=rep(NA,1000),new.r2.adj.crave=rep(NA,1000),
                          new.var.retained.want= rep(NA,1000),new.r2.want=rep(NA,1000), 
                          new.r2.adj.want=rep(NA,1000))
##for drinking events
set.seed(8675309)
new.pred1=predict(n.reg1, newx=testmat, s = "lambda.min", type = "link")
testL0$newregpred1=as.numeric(new.pred1)
new.auc<- pROC::auc(pROC::roc(testL0$drinkbin~testL0$newregpred1))
new.auc
if(is.na(new.auc)==T){
  new.auc <- NA
  new.spec <- NA
  new.sens <- NA
  new.brier <- NA
}
if(exists("new.auc")==F){
  new.auc <- NA
  new.spec <- NA
  new.sens <- NA
  new.brier <- NA
}

if (is.na(new.auc)==F){
  roc <- pROC::coords((pROC::roc(testL0$drinkbin~testL0$newregpred1)),"best", ret=c("threshold", "specificity", "sensitivity"), transpose=FALSE)
  ifelse( nrow(roc[2])>1,new.spec <- NA, new.spec <- roc[2])
  ifelse( nrow(roc[2])>1,new.sens <- NA, new.sens <- roc[3])
  
  new.pred.prob1 <- predict(n.reg1, newx=testmat, s = "lambda.min", type = "response")
  new.brier <- mean((new.pred.prob1-as.numeric(as.character(testL0$drinkbin)))^2)
  new.brier
}
##for craving

new.pred2=predict(n.reg2, newx=testmat, s = "lambda.min", type = "link")
testL0$new.regpred2=as.numeric(new.pred2)
new.r2.crave <- (cor.test(testL0$craving,testL0$new.regpred2)[4]$estimate)^2
new.r2.adj.crave <- 1- (((1-new.r2.crave)*(nrow(testmat)-1))/ 
                          (nrow(testmat) -((sum(n.coef2[2:length(n.coef2)]!=0)) -1)))

##for wanting

new.pred3=predict(n.reg3, newx=testmat, s = "lambda.min", type = "link")
testL0$new.regpred3=as.numeric(new.pred3)
new.r2.want <- (cor.test(testL0$wantdrink,testL0$new.regpred3)[4]$estimate)^2
new.r2.adj.want <- 1- (((1-new.r2.want)*(nrow(testmat)-1))/
                         (nrow(testmat) -((sum(n.coef3[2:length(n.coef3)]!=0))-1)))

new.results[1,] <- cbind(n.var.retained.drink,new.auc,new.spec,new.sens,new.brier,
                         n.var.retained.crave,new.r2.crave, new.r2.adj.crave,
                         n.var.retained.want, new.r2.want, new.r2.adj.want)


#View(new.results)


############################################
###Multilevel version of nomothetic model###
############################################
mlm.results <- data.frame(mlm.auc=rep(NA,1000),mlm.spec = rep(NA,1000),mlm.sens=rep(NA,1000),
                          mlm.r2.crave=rep(NA,1000),mlm.r2.want=rep(NA,1000))

################################
###mlm for drinking occasions###
################################
library(lme4)

drinkmlm = glmer(drinkbin ~ stressed.1 + down.1 + pressure.1 + happy.1 + conflict.1 + craving.1 +
                   impulsive.1 + posexpect.1 + peerpercent.1 + wantdrink.1 + angry.1 + drinkbin.1 +
                   tues + wed + fri + sat + ping2 + ping6 + ping7 + ping8 + cos2T + sin2T + cosW +
                   sinW + (1|ID), family = binomial, n.train)
summary(drinkmlm)

mlm.test <- read.csv("n.testp001.csv")
drinkpred2 = predict(drinkmlm, mlm.test)
mlm.test$drinkpred2 = as.numeric(drinkpred2)
mlm.results$mlm.auc[1] <- pROC::auc(pROC::roc(mlm.test$drinkbin ~ mlm.test$drinkpred2))

roc <- pROC::coords((pROC::roc(mlm.test$drinkbin~mlm.test$drinkpred2)),"best", ret=c("threshold", "specificity", "sensitivity"), transpose=FALSE)
mlm.results$mlm.spec[1]<- as.numeric(roc[2])
mlm.results$mlm.sens[1]<- as.numeric(roc[3])

mlmpred = predict(drinkmlm, mlm.test,  type = "response") 
testL0$mlmpred = as.numeric(mlmpred)
mlm.results$mlm.brier<- DescTools::BrierScore(testL0$drinkbin, testL0$mlmpred)


#######################
### mlm for craving ###
#######################

cravemlm = lmer(craving ~ pressure.1 + enthusiastic.1 + conflict.1 + craving.1 + posexpect.1 +
                  peerpercent.1 + wantdrink.1 + drinkbin.1 + fri + sat + ping3 + ping6 +
                  ping7 + sinT + sin2T +(1|ID), n.train)
summary(cravemlm)

#R2 for intercept-only model
cravepred2 = predict(cravemlm, mlm.test)
mlm.test$cravepred2 = as.numeric(cravepred2)
mlm.results$mlm.r2.crave[1]<- (cor.test(mlm.test$craving, mlm.test$cravepred2)[4]$estimate)^2

#######################
### mlm for wanting ###
#######################

wantmlm = lmer(wantdrink ~ stressed.1 + pressure.1 + enthusiastic.1 + conflict.1 +
                 craving.1 + posexpect.1 + peerpercent.1 + wantdrink.1 + delay_grat.1 +
                 angry.1 + drinkbin.1 + mon + tues + fri + sat + sun + ping2 + ping6 +
                 ping7 + ping8 + cub + sin2T + sinW +(1|ID), n.train)
summary(wantmlm)

#R2 for intercept-only model
wantpred2 = predict(wantmlm, mlm.test)
mlm.test$wantpred2 = as.numeric(wantpred2)
mlm.results$mlm.r2.want[1]<- (cor.test(mlm.test$wantdrink, mlm.test$wantpred2)[4]$estimate)^2

full.results <- data.frame()
full.results <- cbind(i.results,n.results, new.results, mlm.results)
#View(full.results)


write.csv(full.results[1,], "full.resultsP001.csv")
