# convenience functions for dealing with caret models
library(caret)
library(ROCR)

#generic function for fitting models with caret
fitcaret<-function(Xtrn,ytrn,mth=c("kknn","xgbTree","svmRadialSigma","glmnet","nb"),trc=trainControl("repeatedcv",5,repeats=2)){
  #ytrn<-y[idtr]
  #ytst<-y[-idtr]
  #Xtrn<-X[idtr,]
  #Xtst<-X[-idtr,]
  mth<-match.arg(mth)
  args<-list(Xtrn,ytrn,preProcess="zv",metric="Kappa",trControl=trc,method=mth)
  if(mth=="svmRadialSigma"){ args$prob.model<-TRUE }
  if(mth=="glmnet"){ args$family<-"binomial" }
  if(mth=="nb"){
    args$preProcess<-"nzv"
    tg=expand.grid(fL=c(0,0.5,1.0),usekernel=TRUE,adjust=c(0.5,1.0))
    args$tuneGrid<-tg
    args$prior<-table(ytrn)/length(ytrn)
  }
  do.call(train,args)
}

#for any caret model, extract the true and false positive rates that can be used to make ROC
caret_roc<-function(fit,Xtst,ytst,ytrn){
  p1<-prediction(predict(fit,type="prob")[,2], ytrn)
  p2<-prediction(predict(fit,Xtst,type="prob")[,2], ytst)
  p1<-performance(p1,"tpr","fpr")
  p2<-performance(p2,"tpr","fpr")
  d1<-data.frame(tpr=p1@y.values[[1]],fpr=p1@x.values[[1]],fold="train")
  d2<-data.frame(tpr=p2@y.values[[1]],fpr=p2@x.values[[1]],fold="test")
  rbind(d1,d2)
}

roc2auc<-function(fpr,tpr){ 
  #trapezoid rule for area under the curve
  m<-(tpr[-1]+tpr[-length(tpr)])/2
  sum(m*diff(fpr)) 
}