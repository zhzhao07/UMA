
uma.predict=function(x, y, factorID = NULL, newdata, candi_models, weight,alpha,
                 method, dim)
{
  
  
  if (dim=="L"){
    x    = as.data.frame(x)           ;   x_test       = as.data.frame(newdata) ; names(x_test) = names(x)
    y    = as.matrix(y)               ;   p            = ncol(x)
    n    = length(y)                  ;   candi_models = cdm(candi_models,x)
    n_mo = nrow(candi_models)         ;   sk <- rowSums(candi_models)
    if(is.null(factorID)==FALSE)
               {
      for (i in 1:length(factorID))
               {
      x_test[,(factorID[i])]=factor(x_test[,factorID[i]])
      x[,(factorID[i])]=factor(x[,factorID[i]])
      h=summary(x[, factorID[i]])
      if (min(h)<=(nrow(x)/length(h)/2))
               {
        stop("The level of categorical data appears to be extremely unbalanced, and the number of observations at some levels is too small. It is recommended to try again after deal with categorical data reasonablely.")
               }
               }
               }
    # Check the data
    for(i in 1:p)
          {
      if(all(is.na(x[,i])))
          {
        stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to uarm.")
          }
          }
    dat          = data.frame(y, x)   ;  names(dat) = c( "y",names(x)) 
    F1                = matrix(0, nrow = nrow(x_test), ncol = n_mo)
    bname             = names(glm( y ~., data = dat)$coefficients )
    betahat           = matrix(0,length(bname), n_mo)
    rownames(betahat) = bname
    
    for (i in 1:n_mo) {
      varindex    = (candi_models[i, ] == 1)
      dati        = data.frame(y,x[, varindex])
      names(dati) = c("y",names(x)[varindex])
      LSL         = glm(as.numeric(y) ~., data=dati)
      betahat[names((LSL)$coefficients),i] = (LSL)$coefficients
      
      tedat        = data.frame(x_test[, varindex])
      names(tedat) = c(names(x)[varindex])
      F1[, i]      = predict(LSL, newdata = tedat)
                      }    
    if (method=="L1-UARM"|method=="UARM"){
      #gbm
      gbmglm       = gbm(y ~ ., data = dat, distribution = "gaussian", n.minobsinnode = 0)
      pre_gbm      = predict(gbmglm , n.trees = min(100, nrow(x_test)/2), newdata = x_test)  
      
      #l2b
      l2bglm       = glmboost(y~., data= dat ,center=F)
      pre_l2b      = predict( l2bglm, newdata = x_test )
      
      #rf
      rfglm        = randomForest(y~ ., data= dat)
      pre_rf       = predict( rfglm , newdata = x_test )
      
      #bag
      bagglm       = bagging(y ~ ., data=dat , coob=TRUE)
      pre_bag      = predict( bagglm, newdata = x_test )
      
      #bart
      pre_bart     = bart( x, as.double(y), x_test, ndpost=200, verbose=F)$yhat.test.mean
      
      
      FU1         = cbind(pre_gbm, pre_l2b, pre_rf, pre_bag, pre_bart, F1)
      pre_out     = FU1%*%weight
      #betahat=betahat
      }
    if (method=="BMA"){
      bcJ              = bicreg(x, y, strict = FALSE, OR = 20) 
      pre_out          = predict(bcJ,x_test)$mean
                      }
    if (method=="MMA"|method=="JMA"|method=="SAIC"|method=="SBIC"|method=="SFIC"|method=="ARM"|method=="L1-ARM")
    {
      pre_out =as.matrix(F1 %*% weight)
    }
    }
  
  
  if(dim=="H"){
    if (method=="MCV"){
      x_test=newdata
      q=candi_models$q
      M=candi_models$M
      MaxM=candi_models$MaxM
      Index=candi_models$Index
      Maxp=candi_models$Maxp
      beta_cv=candi_models$beta_cv
      MUs<- matrix(0,dim(x_test)[1],M)
      for(k in 1:M){
        if(k<MaxM){
          USE <- Index[(q*(k-1)+1):(q*k)]
          Zk <- x_test[,USE]
          MUs[,k] <- Zk%*%beta_cv[[k]]
        }
        if(k==MaxM){
          USE <- Index[(q*(k-1)+1):Maxp]
          Zk <- x_test[,USE]
          MUs[,k] <- Zk%*%beta_cv[[k]]
        }
      }
      pre_out=MUs%*%weight
    } else{
      x_test=as.matrix(newdata)
      p=NCOL(x)
      n=length(y)
      if (is.numeric(x) & is.null(factorID)==TRUE){
        x=x
        x_test=newdata
        data=data.frame(y,x)
        if (n<200){
          boost_lm=gbm(y~.,data=data,distribution="gaussian",n.minobsinnode = 1)
          #pre_boost0=predict(sboost_lm,newdata=data.frame(x1_rf), n.trees=n/2,type = "link")
          pre_boost=predict(boost_lm,newdata=data.frame(x_test), n.trees=n/2,type = "link")
        }else{
          boost_lm=gbm(y~.,data=data,distribution="gaussian",cv.folds = 5)
          best.iter <- gbm.perf(boost_lm, method = "cv")
          #spre_boost0=predict(sboost_lm,newdata=data.frame(x1_rf), n.trees=best.iter,type = "link")
          pre_boost=predict(boost_lm,newdata=data.frame(x_test), n.trees=best.iter,type = "link")
        }
        ##############################(3)L2Boost
        Gboost_lm=glmboost(y=as.numeric(y),x,center=FALSE)
        pre_Gboost=predict(Gboost_lm,newdata=x_test)
        ##################(5)Random forest
        randomf=randomForest(y~., data=data)
        pre_rf=predict(randomf,newdata=data.frame(x_test))
        
        
      } else {
        xrf=x
        xrf_test=newdata
        x0=x
        x0_test=newdata
        if(is.numeric(factorID)){
          cat=factorID
        } else
        { cat=which(colnames(x)==factorID)}
        
        Xm=c(0)
        XX=list()
        Xt=c(0)
        XXt=list()
        for(j in 1:dim(x0)[2]){
          if(j %in% cat){
            Xc=class.ind(as.matrix(x0[,cat[which(j==cat)]]))###dummy转换
            XX[[j]]=Xc[,-1]
            Xct=class.ind(as.matrix(x0_test[,cat[which(j==cat)]]))###dummy转换
            XXt[[j]]=Xct[,-1]
          } else{
            XX[[j]]=x0[,j]
            XXt[[j]]=x0_test[,j]
          }
          Xm=cbind(Xm,XX[[j]])###dummy covariable
          Xt=cbind(Xt,XXt[[j]])###dummy covariable
        }
        x=Xm[,-1]
        x_test=Xt[,-1]
        p=NCOL(x)
        #############3 non-parametric methods
        data_rf=data.frame(y,xrf)
        if (n<200){
          boost_lm=gbm(y~.,data=data_rf,distribution="gaussian",n.minobsinnode = 1)
          #pre_boost0=predict(sboost_lm,newdata=data.frame(x1_rf), n.trees=n/2,type = "link")
          pre_boost=predict(boost_lm,newdata=data.frame(xrf_test), n.trees=n/2,type = "link")
        }else{
          boost_lm=gbm(y~.,data=data_rf,distribution="gaussian",cv.folds = 5)
          best.iter <- gbm.perf(boost_lm, method = "cv")
          #spre_boost0=predict(sboost_lm,newdata=data.frame(x1_rf), n.trees=best.iter,type = "link")
          pre_boost=predict(boost_lm,newdata=data.frame(xrf_test), n.trees=best.iter,type = "link")
        }
        
        ##############################(3)L2Boost
        Gboost_lm=glmboost(y=as.numeric(y),x,center=FALSE)
        pre_Gboost=predict(Gboost_lm,newdata=x_test)
        ##################(5)Random forest
        data_rf=data.frame(y,xrf)
        randomf=randomForest(y~., data=data_rf)
        pre_rf=predict(randomf,newdata=data.frame(xrf_test))
      }
      
      n_mo=NROW(candi_models)##the number of candidates
      sk=rowSums(candi_models)##the dimension of candidates
      
      if (any(candi_models[1, ] == 1)) {
        betahat=matrix(0,p+1,n_mo)
        for (j in seq(n_mo)) {
          varindex=which(candi_models[j, ] == 1)
          glmfit=lm(y~x[,varindex])
          betahat[c(1,varindex+1),j]=glmfit$coefficients
        }
      } else{
        betahat=matrix(0,p+1,n_mo)
        betahat[1,1]=mean(y)
        for (j in 2:n_mo) {
          varindex=which(candi_models[j, ] == 1)
          glmfit=lm(y ~ x[,varindex])
          betahat[c(1,varindex+1),j]=glmfit$coefficients
        }
      }
      
      if (method=="UARM"){
        F1=cbind(1,x_test)%*%betahat
        pre_out=cbind(pre_boost,pre_Gboost, pre_rf,F1)%*%weight
      }
      
      if (method=="L1-UARM"){
        F1=cbind(1,x_test)%*%betahat
        pre_out=cbind(pre_boost,pre_Gboost, pre_rf,F1)%*%weight
      }
      
      
      if (method=="ARM"){
        F1=cbind(1,x_test)%*%betahat
        pre_out=F1%*%weight
      }
      
      if (method=="L1-ARM"){
        F1=cbind(1,x_test)%*%betahat
        pre_out=F1%*%weight
      }
      
      if (method=="UARM.rf"){
        F1=cbind(1,x_test)%*%betahat
        pre_out=cbind(pre_boost,pre_Gboost, pre_rf,F1)%*%weight
      }
      
      if (method=="L1-UARM.rf"){
        F1=cbind(1,x_test)%*%betahat
        pre_out=cbind(pre_boost,pre_Gboost, pre_rf,F1)%*%weight
      }
      
      if (method=="SAICp"){
        F1=cbind(1,x_test)%*%betahat
        pre_out=F1%*%weight
      }
      
      if (method=="SBICp"){
        F1=cbind(1,x_test)%*%betahat
        pre_out=F1%*%weight
      }
      
      if (method=="PMA"){
        F1=cbind(1,x_test)%*%betahat
        pre_out=F1%*%weight
      }}
  }
  
  object<-list(pre_out=pre_out)
  class(object) <- "uma.predict"
  object
}
