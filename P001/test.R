results <- data.frame(drinks.test= rep(NA,1000),var.retained.drink= rep(NA,1000),auc=rep(NA,1000),
                      spec = rep(NA,1000),sens=rep(NA,1000),brier=rep(NA,1000),
                      var.retained.crave= rep(NA,1000),r2.crave=rep(NA,1000),
                      r2.adj.crave=rep(NA,1000),var.retained.want= rep(NA,1000),
                      r2.want=rep(NA,1000), r2.adj.want=rep(NA,1000))
set.seed(8675309)

 for (j in 1:100){ 
 
  train_ind <- sample(seq_len(nrow(datyL1)), size = smp_size)
  
  #Training and testing sets
  trainL0 <- datyL0[train_ind, ]
  testL0 <- datyL0[-train_ind, ]
  trainL1 <- datyL1[train_ind, ]
  testL1 <- datyL1[-train_ind, ]
  
  #Visual spot check for drinking events
  #Events to be predicted are located in L0 (lag-zero/contemporaneous) data
  even=seq(0,nrow(testL0)-1,1)
  plot(drinkbin~even,trainL0,type='o')
  plot(drinkbin~even,testL0,type='o')
  #plot(drinkbin~even,trainL1,type='o')
  #plot(drinkbin~even,testL1,type='o')
  
  test.drinks<- sum(testL0$drinkbin)
  train.drinks <- sum(trainL0$drinkbin)
  if(test.drinks<2){next}
  if(train.drinks<2){next}
  
  ######################################
  ########  Drinking Events  ###########
  ######################################
  
  #Prediction variables/Prediction matrices (feature space)
  #Vars from datyL1 are lagged, vars from datyL0 are contemporaneous with DV
  #Psychosocial variables (left side) + time variables (right side)
  trainmat = data.matrix(cbind(trainL1[,c(4:18,26)],trainL0[,c(29:35,37:53)]))
  testmat = data.matrix(cbind(testL1[,c(4:18,26)],testL0[,c(29:35,37:53)]))
  
  alpha.set <- seq(from=0.5, to=0, by=-.05)
  
  for (i in alpha.set){
    
    reg1=cv.glmnet(trainmat,trainL0$drinkbin,family='binomial',standardize=T,alpha=i)
    coef1=coef(reg1, s = "lambda.min")
    if ((sum(coef1[2:length(coef1)]!=0))>=1){ 
      alpha.auc <- i
      break
    } 
  }
  alpha.auc
  coef1
  var.retained.drink<- sum(coef1[2:length(coef1)]!=0)
  pred1=predict(reg1, newx=testmat, s = "lambda.min", type = "link")
  testL0$regpred1=as.numeric(pred1)
  auc<- pROC::auc(pROC::roc(testL0$drinkbin~testL0$regpred1))
  auc
  
  roc <- pROC::coords((pROC::roc(testL0$drinkbin~testL0$regpred1)),"best", ret=c("threshold", "specificity", "sensitivity"), transpose=FALSE)
  ifelse( nrow(roc[2])>1,spec <- NA, spec <- roc[2])
  ifelse( nrow(roc[2])>1,sens <- NA, sens <- roc[3])
  
  pred.prob1 <- predict(reg1, newx=testmat, s = "lambda.min", type = "response")
  brier <- mean((pred.prob1-as.numeric(as.character(testL0$drinkbin)))^2)
  brier
  
  
  ##############################
  ########  Craving  ###########
  ##############################
  
  plot(craving~even,trainL0,type='o')
  
  for (i in alpha.set){
    
    reg2=cv.glmnet(trainmat,trainL0$craving,family='gaussian',standardize=T,alpha=i)
    coef2=coef(reg2, s = "lambda.min")
    if ((sum(coef2[2:length(coef2)]!=0))>=1){ 
      alpha.crave <- i
      break
    } 
  }
  alpha.crave
  coef2
  var.retained.crave<- sum(coef2[2:length(coef2)]!=0)
  pred2=predict(reg2, newx=testmat, s = "lambda.min", type = "link")
  testL0$regpred2=as.numeric(pred2)
  r2.crave <- (cor.test(testL0$craving,testL0$regpred2)[4]$estimate)^2
  r2.adj.crave <- 1- (((1-r2.crave)*(nrow(testmat)-1))/ (nrow(testmat) -((sum(coef2[2:length(coef2)]!=0)) -1)))
  
  
  ############################################
  ########  I Would Like to Drink  ###########
  ############################################
  
  plot(wantdrink~even,trainL0,type='o')
  
  
  for (i in alpha.set){
    
    reg3=cv.glmnet(trainmat,trainL0$wantdrink,family='gaussian',standardize=T,alpha=i)
    coef3=coef(reg3, s = "lambda.min")
    if ((sum(coef3[2:length(coef3)]!=0))>=1){ 
      alpha.want <- i
      break
    } 
  }
  alpha.want
  coef3
  var.retained.want<- sum(coef3[2:length(coef3)]!=0)
  pred3=predict(reg3, newx=testmat, s = "lambda.min", type = "link")
  testL0$regpred3=as.numeric(pred3)
  r2.want <- (cor.test(testL0$wantdrink,testL0$regpred3)[4]$estimate)^2
  r2.adj.want <- 1- (((1-r2.want)*(nrow(testmat)-1))/ (nrow(testmat) -((sum(coef2[2:length(coef3)]!=0))
                                                                       -1)))
  
  results[j,] <- cbind(test.drinks,var.retained.crave,auc,spec,sens,brier,var.retained.crave,
                       r2.crave, r2.adj.crave,var.retained.want, r2.want, r2.adj.want)
  
}
View(results)
