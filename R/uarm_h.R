###uarm_h provide the calculation of universal model averging methods for high-dimensional data.
uarm_h=function(x,y,factorID=NULL,candidate='H4',candi_models, n_train=ceiling(n/2),method='L1-UARM',
                no_rep=20, p0=0.5, psi=1, prior = TRUE) {
  # Check the data
  if(all(is.na(x)))
  {
    stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) 
         to eliminate missing data before passing X and y to calculate.")
  }
  else if(is.numeric(x) & is.null(factorID)==TRUE){
    n=length(y)
    p=NCOL(x)
    #################candidates
    if(candidate=='H4'){
      cand=SOIL(x,y,n_train=ceiling(n/2),no_rep=100,family="gaussian",method="union",weight_type="ARM",prior = TRUE)
      candi_models=cand$candidate_models
    }###4 methods
    if(candidate=='H2'){
      cand=SOIL(x,y,n_train=ceiling(n/2),no_rep=100,family="gaussian",method=="lasso",weight_type="ARM",prior = TRUE)
      candi_models=cand$candidate_models
    }###2 methods
    if(candidate=='H1'){
      thelasso.cv=cv.glmnet(x,y,family = "gaussian",alpha=1) ## first stage ridge
      bhat=as.matrix(coef(thelasso.cv,s="lambda.1se"))[-1,1] ## coef() is a sparseMatrix
      if(all(bhat==0)){
        bhat=rep(.Machine$double.eps*2,length(bhat))
      }
      adpen=(1/pmax(abs(bhat),.Machine$double.eps)) ## the adaptive lasso weight
      m2=glmnet(x,y,family = "gaussian",alpha=1,exclude=which(bhat==0),penalty.factor=adpen)
      #############################
      m2.path=as.matrix(m2$beta)
      beta.path=t(cbind(m2.path))
      candidate_models=(1-(beta.path == 0))
      candi_models=unique(candidate_models)
    }###1 method
    if(candidate=='H0'){
      candi_models=candi_models
    }

    n_mo=NROW(candi_models)##the number of candidates
    sk=rowSums(candi_models)##the dimension of candidates

    ck_compute=function(n_mo, sk, p) {
      ck=2*log(sk + 2) + ifelse(sk > 0, sk*log(exp(1)*p/sk), 0)
      return(ck)
    }
    #################betahat for candidates in regression
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
  ###########################method=UARM  
    if (method=='UARM'){
      wt_calc=function(rep_id) {
        lw_num <- matrix(0,nrow=1, 3+n_mo)
        tindex <- sample(n,n_train, replace = FALSE)
          x1=x[tindex,];x2=x[-tindex,]
          y1=y[tindex];y2=y[-tindex]
          ######5 methods
          sdata_boost=data.frame(y1,x1)
          if (n<200){
            sboost_lm=gbm(y1~.,data=sdata_boost,distribution="gaussian",n.minobsinnode = 1)
            spre_boost0=predict(sboost_lm,newdata=data.frame(x1), n.trees=n/2,type = "link")
            spre_boost1=predict(sboost_lm,newdata=data.frame(x2), n.trees=n/2,type = "link")
          }else{
            sboost_lm=gbm(y1~.,data=sdata_boost,distribution="gaussian",cv.folds = 5)
            best.iter <- gbm.perf(sboost_lm, method = "cv")
            spre_boost0=predict(sboost_lm,newdata=data.frame(x1), n.trees=best.iter,tye = "link")
            spre_boost1=predict(sboost_lm,newdata=data.frame(x2), n.trees=best.iter,type = "link")
          }
          sigma1=sqrt(sum((y1-spre_boost0)^2)/n_train)
          d1=sum((y2-spre_boost1)^2)
          lw_num[1]=(-(n-n_train))*log(sigma1)-((sigma1)^(-2))*d1/2
          ######
          ######
          sGboost_lm=glmboost(y=as.numeric(y1),x=x1,center=F)
          spre_Gboost0=predict(sGboost_lm,newdata=x1)
          spre_Gboost1=predict(sGboost_lm,newdata=x2)
          sigma2=sqrt(sum((y1-spre_Gboost0)^2)/n_train)
          d2=sum((y2-spre_Gboost1)^2)
          lw_num[2]=(-(n-n_train))*log(sigma2) - ((sigma2)^(-2))*d2/2
          ######
          data_rf=data.frame(y1,x1)
          srandomf=randomForest(y1 ~ ., data=sdata_boost)
          spre_rf0=predict(srandomf,newdata=data.frame(x1))
          spre_rf1=predict(srandomf,newdata=data.frame(x2))
          sigma3=sqrt(sum((y1-spre_rf0)^2)/n_train)
          d3=sum((y2-spre_rf1)^2)
          lw_num[3]=(-(n-n_train))*log(sigma3)-((sigma3)^(-2))*d3/2

          if (any(candi_models[1, ] == 1)) {
            for (j in seq(n_mo)) {
              varindex=(candi_models[j, ] == 1)
              glmfit =lm(y1 ~ x1[,varindex])
              sigma_k =summary(glmfit)$sigma
              if(any(is.na(glmfit$coef))) {
                lw_num[j+3] = -Inf
              } else {
                dk =sum((y2-cbind(1, x2[,varindex])%*%glmfit$coef)^2)
                lw_num[j+3]=(-(n-n_train))*log(sigma_k)-((sigma_k)^(-2))*dk/2
              }
            }
          } else {
            d1=sum((y2-mean(y1))^2)
            sigma_1=sd(y1)
            lw_num[3+1] =(-(n-n_train))*log(sigma_1)-((sigma_1)^(-2))*d1/2
            for (j in seq(2, n_mo)) {
              varindex=(candi_models[j, ] == 1)
              glmfit =lm(y1 ~ x1[,varindex])
              sigma_k=summary(glmfit)$sigma
              if (any(is.na(glmfit$coef))) {
                lw_num[j+3]=-Inf
              } else {
                dk=sum((y2 - cbind(1, x2[,varindex])%*%glmfit$coef)^2)
                lw_num[j+3] =(-(n-n_train))*log(sigma_k)-((sigma_k)^(-2))*dk/2
              }
            }
          }
        return(lw_num)
      }
      lw_num=matrix(unlist(mclapply(seq(no_rep), wt_calc)), nrow = no_rep, ncol = 3+n_mo, byrow = TRUE)
      if (prior == TRUE) {
        ck0=-log(1-p0)+log(3)
        ck1=-log(p0)+ck_compute(n_mo, sk, p) 
        ck=c(rep(ck0,3),ck1)
        lw_num=sweep(lw_num, MARGIN = 2, psi*ck, "-")
      }
      lw_num=sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max), "-")
      w_num=exp(lw_num)
      weight=colMeans(w_num/rowSums(w_num))
      weight_se=apply(w_num,2,sd)/sqrt(no_rep)
    }

  ###########################method=L1-UARM 
    if (method=='L1-UARM'){
      wt_calc <- function(rep_id) {
        lw_num <- matrix(0,nrow=1, 3+n_mo)
        tindex <- sample(n,n_train, replace = FALSE)
          x1=x[tindex,];x2=x[-tindex,]
          y1=y[tindex];y2=y[-tindex]
          ######
          sdata_boost=data.frame(y1,x1)
          if (n<200){
            sboost_lm=gbm(y1~.,data=sdata_boost,distribution="gaussian",n.minobsinnode = 1)
            spre_boost0=predict(sboost_lm,newdata=data.frame(x1), n.trees=n/2,type = "link")
            spre_boost1=predict(sboost_lm,newdata=data.frame(x2), n.trees=n/2,type = "link")
          }else{
            sboost_lm=gbm(y1~.,data=sdata_boost,distribution="gaussian",cv.folds = 5)
            best.iter <- gbm.perf(sboost_lm, method = "cv")
            spre_boost0=predict(sboost_lm,newdata=data.frame(x1), n.trees=best.iter,type = "link")
            spre_boost1=predict(sboost_lm,newdata=data.frame(x2), n.trees=best.iter,type = "link")
          }
          sigmat1=sum(abs(y1-spre_boost0))/n_train
          dt1=sum(abs(y2-spre_boost1))
          lw_num[1] <- (-(n-n_train)) *log(sigmat1)-dt1/sigmat1
          ######
          sGboost_lm=glmboost(y=as.numeric(y1),x=x1,center=F)
          spre_Gboost0=predict(sGboost_lm,newdata=x1)
          spre_Gboost1=predict(sGboost_lm,newdata=x2)
          sigmat2=sum(abs(y1-spre_Gboost0))/n_train
          dt2=sum(abs(y2-spre_Gboost1))
          lw_num[2] <- (-(n-n_train))*log(sigmat2)-dt2/sigmat2
          ######
          data_rf=data.frame(y1,x1)
          srandomf=randomForest(y1 ~ ., data=sdata_boost)
          spre_rf0=predict(srandomf,newdata=data.frame(x1))
          spre_rf1=predict(srandomf,newdata=data.frame(x2))
          sigmat3=sum(abs(y1-spre_rf0))/n_train
          dt3=sum(abs(y2-spre_rf1))
          lw_num[3] <- (-(n-n_train))*log(sigmat3)-dt3/sigmat3

          if (any(candi_models[1, ] == 1)) {
            for (j in 1:n_mo) {
              varindex <- which(candi_models[j, ] == 1)
              glmfit <- lm(y1 ~ x1[,varindex])
              dk1=sum(abs(y1-cbind(1,x1[,varindex])%*%glmfit$coefficients))/(n_train)
              if(any(is.na(glmfit$coef))) {
                lw_num[j+3] <- -Inf
              } else {
                Dk1 <- sum(abs(y2-cbind(1,x2[,varindex])%*%glmfit$coefficients))
                lw_num[j+3] <- (-(n-n_train))*log(dk1)-Dk1/dk1
              }
            }
          } else {
            Dk0 <- sum(abs(y2 - mean(y1)))
            dk0<-sum(abs(y1-mean(y1)))/(n_train)
            lw_num[3+1] <- (-(n-n_train))*log(dk0)-Dk0/dk0
            for (j in 2:n_mo) {
              varindex <- which(candi_models[j, ] == 1)
              glmfit <- lm(y1 ~ x1[,varindex])
              dk1=sum(abs(y1-cbind(1,x1[,varindex])%*%glmfit$coefficients))/(n_train)
              if (any(is.na(glmfit$coef))) {
                lw_num[j+3] <- -Inf
              } else {
                Dk1 <-sum(abs(y2-cbind(1,x2[,varindex])%*%glmfit$coefficients))
                lw_num[j+3] <- (-(n-n_train))*log(dk1)-Dk1/dk1
              }
            }
          }
        return(lw_num)
      }
      lw_num= matrix(unlist(mclapply(seq(no_rep), wt_calc)), nrow = no_rep, ncol = 3+n_mo, byrow = TRUE)
      if (prior==TRUE) {
        ck0=-log(1-p0)+log(3)
        ck1=-log(p0)+ck_compute(n_mo, sk, p)
        ck=c(rep(ck0,3),ck1)
        lw_num=sweep(lw_num, MARGIN = 2, psi*ck, "-")
      }
      lw_num=sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max), "-")
      w_num=exp(lw_num)
      weight=colMeans(w_num/rowSums(w_num))
      weight_se=apply(w_num,2,sd)/sqrt(no_rep)
    }

  ###########################method=UARM.rf
    if (method=='UARM.rf'){
      lscvrf=function(x, y,candi_models, n_train,Rep) {
        p=NCOL(x)
        n=length(y)
        n_mo=NROW(candi_models)
        sk=rowSums(candi_models)
        sigrf1=c()
        for(j in 1:Rep){
          tindex=sample(n, n_train, replace = FALSE)
          x1=x[tindex,];x2=x[-tindex,]
          y1=y[tindex];y2=y[-tindex]
          #####
          sdata_boost=data.frame(y1,x1)
          srandomf=randomForest(y1 ~ ., data=sdata_boost, importance=TRUE,proximity=TRUE)
          spre_rf0=predict(srandomf,newdata=data.frame(x1))
          spre_rf1=predict(srandomf,newdata=data.frame(x2))
          sigrf1[j]=sqrt(sum((y2-spre_rf1)^2)/n_train)
        }
        Sigma1=mean(sigrf1)
        list(Sigma1=Sigma1)
      }
      iscv=lscvrf(x=X,y=y,candi_models,n_train,Rep=20)
      sig1=iscv$Sigma1 ####for same sigmahat estimate
      wt_calc=function(rep_id) {
        lw_num <- matrix(0,nrow=1, 3+n_mo)
        tindex <- sample(n,n_train, replace = FALSE)
          x1=x[tindex,];x2=x[-tindex,]
          y1=y[tindex];y2=y[-tindex]
          ######3 methods
          sdata_boost=data.frame(y1,x1)
          if (n<200){
            sboost_lm=gbm(y1~.,data=sdata_boost,distribution="gaussian",n.minobsinnode = 1)
            #spre_boost0=predict(sboost_lm,newdata=data.frame(x1), n.trees=n/2,type = "link")
            spre_boost1=predict(sboost_lm,newdata=data.frame(x2), n.trees=n/2,type = "link")
          }else{
            sboost_lm=gbm(y1~.,data=sdata_boost,distribution="gaussian",cv.folds = 5)
            best.iter <- gbm.perf(sboost_lm, method = "cv")
            #spre_boost0=predict(sboost_lm,newdata=data.frame(x1), n.trees=best.iter,type = "link")
            spre_boost1=predict(sboost_lm,newdata=data.frame(x2), n.trees=best.iter,type = "link")
          }
          sigma1=sig1
          d1=sum((y2-spre_boost1)^2)
          lw_num[1]=(-(n-n_train))*log(sigma1)-((sigma1)^(-2))*d1/2
          ######
          ######
          sGboost_lm=glmboost(y=as.numeric(y1),x=x1,center=F)
          #spre_Gboost0=predict(sGboost_lm,newdata=x1)
          spre_Gboost1=predict(sGboost_lm,newdata=x2)
          sigma2=sig1
          d2=sum((y2-spre_Gboost1)^2)
          lw_num[2]=(-(n-n_train))*log(sigma2) - ((sigma2)^(-2))*d2/2
          ######
          data_rf=data.frame(y1,x1)
          srandomf=randomForest(y1 ~ ., data=sdata_boost)
          #spre_rf0=predict(srandomf,newdata=data.frame(x1))
          spre_rf1=predict(srandomf,newdata=data.frame(x2))
          sigma3=sig1
          d3=sum((y2-spre_rf1)^2)
          lw_num[3]=(-(n-n_train))*log(sigma3)-((sigma3)^(-2))*d3/2

          if (any(candi_models[1, ] == 1)) {
            for (j in seq(n_mo)) {
              varindex=(candi_models[j, ] == 1)
              glmfit =lm(y1 ~ x1[,varindex])
              sigma_k =sig1
              if(any(is.na(glmfit$coef))) {
                lw_num[j+3] = -Inf
              } else {
                dk =sum((y2-cbind(1, x2[,varindex])%*%glmfit$coef)^2)
                lw_num[j+3]=(-(n-n_train))*log(sigma_k)-((sigma_k)^(-2))*dk/2
              }
            }
          } else {
            d1=sum((y2-mean(y1))^2)
            sigma_1=sig1
            lw_num[3+1] =(-(n-n_train))*log(sigma_1)-((sigma_1)^(-2))*d1/2
            for (j in seq(2, n_mo)) {
              varindex=(candi_models[j, ] == 1)
              glmfit =lm(y1 ~ x1[,varindex])
              sigma_k=sig1
              if (any(is.na(glmfit$coef))) {
                lw_num[j+3]=-Inf
              } else {
                dk=sum((y2 - cbind(1, x2[,varindex])%*%glmfit$coef)^2)
                lw_num[j+3] =(-(n-n_train))*log(sigma_k)-((sigma_k)^(-2))*dk/2
              }
            }
          }
        return(lw_num)
      }
      lw_num <- matrix(unlist(mclapply(seq(no_rep), wt_calc)), nrow = no_rep, ncol = 3+n_mo, byrow = TRUE)
      if (prior == TRUE) {
        ck0=-log(1-p0)+log(3)
        ck1=-log(p0)+ck_compute(n_mo, sk, p)
        ck=c(rep(ck0,3),ck1)
        lw_num <- sweep(lw_num, MARGIN = 2, psi*ck, "-")
      }
      lw_num <- sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max), "-")
      w_num <- exp(lw_num)
      weight <- colMeans(w_num/rowSums(w_num))
      weight_se=apply(w_num,2,sd)/sqrt(no_rep)
    }

  ###########################method=L1-UARM.rf 
    if (method=='L1-UARM.rf'){
      lscvrf=function(x, y,candi_models, n_train,Rep) {
        p=NCOL(x)
        n=length(y)
        n_mo=NROW(candi_models)
        sk=rowSums(candi_models)
        sigrf2=c()
        for(j in 1:Rep){
          tindex=sample(n, n_train, replace = FALSE)
          x1=x[tindex,];x2=x[-tindex,]
          y1=y[tindex];y2=y[-tindex]
          #####
          sdata_boost=data.frame(y1,x1)
          srandomf=randomForest(y1 ~ ., data=sdata_boost, importance=TRUE,proximity=TRUE)
          spre_rf0=predict(srandomf,newdata=data.frame(x1))
          spre_rf1=predict(srandomf,newdata=data.frame(x2))
          sigrf2[j]=sum(abs(y2-spre_rf1))/n_train
        }
        Sigma2=mean(sigrf2)
        list(Sigma2=Sigma2)
      }
      iscv=lscvrf(x=X,y=y,candi_models,n_train,Rep=20)
      sig2=iscv$Sigma2
      wt_calc <- function(rep_id) {
        lw_num <- matrix(0,nrow=1, 3+n_mo)
        tindex <- sample(n,n_train, replace = FALSE)

          x1=x[tindex,];x2=x[-tindex,]
          y1=y[tindex];y2=y[-tindex]
          ######
          sdata_boost=data.frame(y1,x1)
          if (n<200){
            sboost_lm=gbm(y1~.,data=sdata_boost,distribution="gaussian",n.minobsinnode = 1)
            #spre_boost0=predict(sboost_lm,newdata=data.frame(x1), n.trees=n/2,type = "link")
            spre_boost1=predict(sboost_lm,newdata=data.frame(x2), n.trees=n/2,type = "link")
          }else{
            sboost_lm=gbm(y1~.,data=sdata_boost,distribution="gaussian",cv.folds = 5)
            best.iter <- gbm.perf(sboost_lm, method = "cv")
            #spre_boost0=predict(sboost_lm,newdata=data.frame(x1), n.trees=best.iter,type = "link")
            spre_boost1=predict(sboost_lm,newdata=data.frame(x2), n.trees=best.iter,type = "link")
          }
          sigmat1=sig2
          dt1=sum(abs(y2-spre_boost1))
          lw_num[1] <- (-(n-n_train)) *log(sigmat1)-dt1/sigmat1
          ######
          sGboost_lm=glmboost(y=as.numeric(y1),x=x1,center=F)
          spre_Gboost0=predict(sGboost_lm,newdata=x1)
          spre_Gboost1=predict(sGboost_lm,newdata=x2)
          sigmat2=sig2
          dt2=sum(abs(y2-spre_Gboost1))
          lw_num[2] <- (-(n-n_train))*log(sigmat2)-dt2/sigmat2
          ######
          data_rf=data.frame(y1,x1)
          srandomf=randomForest(y1 ~ ., data=sdata_boost)
          spre_rf0=predict(srandomf,newdata=data.frame(x1))
          spre_rf1=predict(srandomf,newdata=data.frame(x2))
          sigmat3=sig2
          dt3=sum(abs(y2-spre_rf1))
          lw_num[3] <- (-(n-n_train))*log(sigmat3)-dt3/sigmat3

          if (any(candi_models[1, ] == 1)) {
            for (j in 1:n_mo) {
              varindex <- which(candi_models[j, ] == 1)
              glmfit <- lm(y1 ~ x1[,varindex])
              dk1=sig2
              if(any(is.na(glmfit$coef))) {
                lw_num[j+3] <- -Inf
              } else {
                Dk1 <- sum(abs(y2-cbind(1,x2[,varindex])%*%glmfit$coefficients))
                lw_num[j+3] <- (-(n-n_train))*log(dk1)-Dk1/dk1
              }
            }
          } else {
            Dk0 <- sum(abs(y2 - mean(y1)))
            dk0<-sig2
            lw_num[3+1] <- (-(n-n_train))*log(dk0)-Dk0/dk0
            for (j in 2:n_mo) {
              varindex <- which(candi_models[j, ] == 1)
              glmfit <- lm(y1 ~ x1[,varindex])
              dk1=sig2
              if (any(is.na(glmfit$coef))) {
                lw_num[j+3] <- -Inf
              } else {
                Dk1 <-sum(abs(y2-cbind(1,x2[,varindex])%*%glmfit$coefficients))
                lw_num[j+3] <- (-(n-n_train))*log(dk1)-Dk1/dk1
              }
            }
          }
        return(lw_num)
      }
      lw_num <- matrix(unlist(mclapply(seq(no_rep), wt_calc)), nrow = no_rep, ncol = 3+n_mo, byrow = TRUE)
      if (prior == TRUE) {
        ck0=-log(1-p0)+log(3)
        ck1=-log(p0)+ck_compute(n_mo, sk, p)
        ck=c(rep(ck0,3),ck1)
        lw_num <- sweep(lw_num, MARGIN = 2, psi*ck, "-")
      }
      lw_num <- sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max), "-")
      w_num <- exp(lw_num)
      weight=colMeans(w_num/rowSums(w_num))
      weight_se=apply(w_num,2,sd)/sqrt(no_rep)
    }
    
    names(weight) = c("GBM", "L2B", "RF", paste("LM", 1:(length(weight) - 3), sep = ""))
    list(weight = weight,weight_se=weight_se,candi_models=candi_models)
  }else {
    candidate='CH3'
    x_rf=x
    n=length(y)
    if(is.numeric(factorID)){
      cat=factorID
    } else
    { cat=which(colnames(x)==factorID)
    }
  
    lengthp=c()
    for(i in 1:length(cat)){
      lengthp[i]=length(summary(x[,cat[i]]))-1
    }
    leng=sum(lengthp)

    for(i in 1:length(cat)){
      h=summary(as.factor(x[,cat[i]]))

      if (min(h)<=(n/length(h)/2))
      {
        stop("The level of categorical data appears to be extremely unbalanced, and the number of observations at some levels is too small.it is recommended to try again after deal with categorical data reasonablely.")
      }
    }

    Xm=c(0)
    XX=list()
    for(j in 1:dim(x)[2]){
      if(j %in% cat){
        Xc=class.ind(as.matrix(x[,cat[which(j==cat)]]))###dummy转换
        XX[[j]]=Xc[,-1]
      } else{
        XX[[j]]=x[,j]
      }
      Xm=cbind(Xm,XX[[j]])###dummy covariable
    }
    x0=Xm[,-1]
    index_group0=list()
    index_group=c()
    for (i in 1:length(XX)){
      index_group0[[i]]=rep(i,dim(as.matrix(XX[[i]]))[2])
      index_group=c(index_group,index_group0[[i]])
    }
    if(candidate=='CH3'){
    m1=grpreg(x0,y,group=index_group,penalty="grLasso")
    m2=grpreg(x0, y,group=index_group,penalty="grMCP")
    m3=grpreg(x0,y,group=index_group,penalty="grSCAD")
    #############################
    #############################
    m1.path <- as.matrix(m1$beta)
    m2.path <- as.matrix(m2$beta)
    m3.path <- as.matrix(m3$beta)
    beta.path <- t(cbind(m1.path, m2.path,m3.path))
    candidate_models <- (1 - (beta.path == 0))
    candidate_models=unique(candidate_models)
    sk0=rowSums(candidate_models)##the dimension of candidates
    de=which(sk0<=n/2-(leng+1))
    candi_models=candidate_models[c(de),-1]
     }
    x=x0
    p=NCOL(x)
    n_mo=NROW(candi_models)##the number of candidates
    sk=rowSums(candi_models)##the dimension of candidates

    ck_compute=function(n_mo, sk, p) {
      ck=2*log(sk + 2) + ifelse(sk > 0, sk*log(exp(1)*p/sk), 0)
      return(ck)
    }

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

    ###########################method=UARM  
    if (method=='UARM'){
      wt_calc=function(rep_id) {
        lw_num <- matrix(0,nrow=1, 3+n_mo)
        tindex <- sample(n,n_train, replace = FALSE)
        x1=x[tindex,];x2=x[-tindex,]
        y1=y[tindex];y2=y[-tindex]
        ######5 methods
        sdata_boost=data.frame(y1,x1)
        if (n<200){
          sboost_lm=gbm(y1~.,data=sdata_boost,distribution="gaussian",n.minobsinnode = 1)
          spre_boost0=predict(sboost_lm,newdata=data.frame(x1), n.trees=n/2,type = "link")
          spre_boost1=predict(sboost_lm,newdata=data.frame(x2), n.trees=n/2,type = "link")
        }else{
          sboost_lm=gbm(y1~.,data=sdata_boost,distribution="gaussian",cv.folds = 5)
          best.iter <- gbm.perf(sboost_lm, method = "cv")
          spre_boost0=predict(sboost_lm,newdata=data.frame(x1), n.trees=best.iter,tye = "link")
          spre_boost1=predict(sboost_lm,newdata=data.frame(x2), n.trees=best.iter,type = "link")
        }
        sigma1=sqrt(sum((y1-spre_boost0)^2)/n_train)
        d1=sum((y2-spre_boost1)^2)
        lw_num[1]=(-(n-n_train))*log(sigma1)-((sigma1)^(-2))*d1/2
        ######
        ######
        sGboost_lm=glmboost(y=as.numeric(y1),x=x1,center=F)
        spre_Gboost0=predict(sGboost_lm,newdata=x1)
        spre_Gboost1=predict(sGboost_lm,newdata=x2)
        sigma2=sqrt(sum((y1-spre_Gboost0)^2)/n_train)
        d2=sum((y2-spre_Gboost1)^2)
        lw_num[2]=(-(n-n_train))*log(sigma2) - ((sigma2)^(-2))*d2/2
        ######
        data_rf=data.frame(y1,x1)
        srandomf=randomForest(y1 ~ ., data=sdata_boost)
        spre_rf0=predict(srandomf,newdata=data.frame(x1))
        spre_rf1=predict(srandomf,newdata=data.frame(x2))
        sigma3=sqrt(sum((y1-spre_rf0)^2)/n_train)
        d3=sum((y2-spre_rf1)^2)
        lw_num[3]=(-(n-n_train))*log(sigma3)-((sigma3)^(-2))*d3/2
        
        if (any(candi_models[1, ] == 1)) {
          for (j in seq(n_mo)) {
            varindex=(candi_models[j, ] == 1)
            glmfit =lm(y1 ~ x1[,varindex])
            sigma_k =summary(glmfit)$sigma
            if(any(is.na(glmfit$coef))) {
              lw_num[j+3] = -Inf
            } else {
              dk =sum((y2-cbind(1, x2[,varindex])%*%glmfit$coef)^2)
              lw_num[j+3]=(-(n-n_train))*log(sigma_k)-((sigma_k)^(-2))*dk/2
            }
          }
        } else {
          d1=sum((y2-mean(y1))^2)
          sigma_1=sd(y1)
          lw_num[3+1] =(-(n-n_train))*log(sigma_1)-((sigma_1)^(-2))*d1/2
          for (j in seq(2, n_mo)) {
            varindex=(candi_models[j, ] == 1)
            glmfit =lm(y1 ~ x1[,varindex])
            sigma_k=summary(glmfit)$sigma
            if (any(is.na(glmfit$coef))) {
              lw_num[j+3]=-Inf
            } else {
              dk=sum((y2 - cbind(1, x2[,varindex])%*%glmfit$coef)^2)
              lw_num[j+3] =(-(n-n_train))*log(sigma_k)-((sigma_k)^(-2))*dk/2
            }
          }
        }
        return(lw_num)
      }
      lw_num=matrix(unlist(mclapply(seq(no_rep), wt_calc)), nrow = no_rep, ncol = 3+n_mo, byrow = TRUE)
      if (prior == TRUE) {
        ck0=-log(1-p0)+log(3)
        ck1=-log(p0)+ck_compute(n_mo, sk, p) 
        ck=c(rep(ck0,3),ck1)
        lw_num=sweep(lw_num, MARGIN = 2, psi*ck, "-")
      }
      lw_num=sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max), "-")
      w_num=exp(lw_num)
      weight=colMeans(w_num/rowSums(w_num))
      weight_se=apply(w_num,2,sd)/sqrt(no_rep)
    }
    
    ###########################method=L1-UARM 
    if (method=='L1-UARM'){
      wt_calc <- function(rep_id) {
        lw_num <- matrix(0,nrow=1, 3+n_mo)
        tindex <- sample(n,n_train, replace = FALSE)
        x1=x[tindex,];x2=x[-tindex,]
        y1=y[tindex];y2=y[-tindex]
        ######
        sdata_boost=data.frame(y1,x1)
        if (n<200){
          sboost_lm=gbm(y1~.,data=sdata_boost,distribution="gaussian",n.minobsinnode = 1)
          spre_boost0=predict(sboost_lm,newdata=data.frame(x1), n.trees=n/2,type = "link")
          spre_boost1=predict(sboost_lm,newdata=data.frame(x2), n.trees=n/2,type = "link")
        }else{
          sboost_lm=gbm(y1~.,data=sdata_boost,distribution="gaussian",cv.folds = 5)
          best.iter <- gbm.perf(sboost_lm, method = "cv")
          spre_boost0=predict(sboost_lm,newdata=data.frame(x1), n.trees=best.iter,type = "link")
          spre_boost1=predict(sboost_lm,newdata=data.frame(x2), n.trees=best.iter,type = "link")
        }
        sigmat1=sum(abs(y1-spre_boost0))/n_train
        dt1=sum(abs(y2-spre_boost1))
        lw_num[1] <- (-(n-n_train)) *log(sigmat1)-dt1/sigmat1
        ######
        sGboost_lm=glmboost(y=as.numeric(y1),x=x1,center=F)
        spre_Gboost0=predict(sGboost_lm,newdata=x1)
        spre_Gboost1=predict(sGboost_lm,newdata=x2)
        sigmat2=sum(abs(y1-spre_Gboost0))/n_train
        dt2=sum(abs(y2-spre_Gboost1))
        lw_num[2] <- (-(n-n_train))*log(sigmat2)-dt2/sigmat2
        ######
        data_rf=data.frame(y1,x1)
        srandomf=randomForest(y1 ~ ., data=sdata_boost)
        spre_rf0=predict(srandomf,newdata=data.frame(x1))
        spre_rf1=predict(srandomf,newdata=data.frame(x2))
        sigmat3=sum(abs(y1-spre_rf0))/n_train
        dt3=sum(abs(y2-spre_rf1))
        lw_num[3] <- (-(n-n_train))*log(sigmat3)-dt3/sigmat3
        
        if (any(candi_models[1, ] == 1)) {
          for (j in 1:n_mo) {
            varindex <- which(candi_models[j, ] == 1)
            glmfit <- lm(y1 ~ x1[,varindex])
            dk1=sum(abs(y1-cbind(1,x1[,varindex])%*%glmfit$coefficients))/(n_train)
            if(any(is.na(glmfit$coef))) {
              lw_num[j+3] <- -Inf
            } else {
              Dk1 <- sum(abs(y2-cbind(1,x2[,varindex])%*%glmfit$coefficients))
              lw_num[j+3] <- (-(n-n_train))*log(dk1)-Dk1/dk1
            }
          }
        } else {
          Dk0 <- sum(abs(y2 - mean(y1)))
          dk0<-sum(abs(y1-mean(y1)))/(n_train)
          lw_num[3+1] <- (-(n-n_train))*log(dk0)-Dk0/dk0
          for (j in 2:n_mo) {
            varindex <- which(candi_models[j, ] == 1)
            glmfit <- lm(y1 ~ x1[,varindex])
            dk1=sum(abs(y1-cbind(1,x1[,varindex])%*%glmfit$coefficients))/(n_train)
            if (any(is.na(glmfit$coef))) {
              lw_num[j+3] <- -Inf
            } else {
              Dk1 <-sum(abs(y2-cbind(1,x2[,varindex])%*%glmfit$coefficients))
              lw_num[j+3] <- (-(n-n_train))*log(dk1)-Dk1/dk1
            }
          }
        }
        return(lw_num)
      }
      lw_num= matrix(unlist(mclapply(seq(no_rep), wt_calc)), nrow = no_rep, ncol = 3+n_mo, byrow = TRUE)
      if (prior==TRUE) {
        ck0=-log(1-p0)+log(3)
        ck1=-log(p0)+ck_compute(n_mo, sk, p)
        ck=c(rep(ck0,3),ck1)
        lw_num=sweep(lw_num, MARGIN = 2, psi*ck, "-")
      }
      lw_num=sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max), "-")
      w_num=exp(lw_num)
      weight=colMeans(w_num/rowSums(w_num))
      weight_se=apply(w_num,2,sd)/sqrt(no_rep)
    }
    
    ###########################method=UARM.rf
    if (method=='UARM.rf'){
      lscvrf=function(x, y,candi_models, n_train,Rep) {
        p=NCOL(x)
        n=length(y)
        n_mo=NROW(candi_models)
        sk=rowSums(candi_models)
        sigrf1=c()
        for(j in 1:Rep){
          tindex=sample(n, n_train, replace = FALSE)
          x1=x[tindex,];x2=x[-tindex,]
          y1=y[tindex];y2=y[-tindex]
          #####
          sdata_boost=data.frame(y1,x1)
          srandomf=randomForest(y1 ~ ., data=sdata_boost, importance=TRUE,proximity=TRUE)
          spre_rf0=predict(srandomf,newdata=data.frame(x1))
          spre_rf1=predict(srandomf,newdata=data.frame(x2))
          sigrf1[j]=sqrt(sum((y2-spre_rf1)^2)/n_train)
        }
        Sigma1=mean(sigrf1)
        list(Sigma1=Sigma1)
      }
      iscv=lscvrf(x=X,y=y,candi_models,n_train,Rep=20)
      sig1=iscv$Sigma1 ####for same sigmahat estimate
      wt_calc=function(rep_id) {
        lw_num <- matrix(0,nrow=1, 3+n_mo)
        tindex <- sample(n,n_train, replace = FALSE)
        x1=x[tindex,];x2=x[-tindex,]
        y1=y[tindex];y2=y[-tindex]
        ######3 methods
        sdata_boost=data.frame(y1,x1)
        if (n<200){
          sboost_lm=gbm(y1~.,data=sdata_boost,distribution="gaussian",n.minobsinnode = 1)
          #spre_boost0=predict(sboost_lm,newdata=data.frame(x1), n.trees=n/2,type = "link")
          spre_boost1=predict(sboost_lm,newdata=data.frame(x2), n.trees=n/2,type = "link")
        }else{
          sboost_lm=gbm(y1~.,data=sdata_boost,distribution="gaussian",cv.folds = 5)
          best.iter <- gbm.perf(sboost_lm, method = "cv")
          #spre_boost0=predict(sboost_lm,newdata=data.frame(x1), n.trees=best.iter,type = "link")
          spre_boost1=predict(sboost_lm,newdata=data.frame(x2), n.trees=best.iter,type = "link")
        }
        sigma1=sig1
        d1=sum((y2-spre_boost1)^2)
        lw_num[1]=(-(n-n_train))*log(sigma1)-((sigma1)^(-2))*d1/2
        ######
        ######
        sGboost_lm=glmboost(y=as.numeric(y1),x=x1,center=F)
        #spre_Gboost0=predict(sGboost_lm,newdata=x1)
        spre_Gboost1=predict(sGboost_lm,newdata=x2)
        sigma2=sig1
        d2=sum((y2-spre_Gboost1)^2)
        lw_num[2]=(-(n-n_train))*log(sigma2) - ((sigma2)^(-2))*d2/2
        ######
        data_rf=data.frame(y1,x1)
        srandomf=randomForest(y1 ~ ., data=sdata_boost)
        #spre_rf0=predict(srandomf,newdata=data.frame(x1))
        spre_rf1=predict(srandomf,newdata=data.frame(x2))
        sigma3=sig1
        d3=sum((y2-spre_rf1)^2)
        lw_num[3]=(-(n-n_train))*log(sigma3)-((sigma3)^(-2))*d3/2
        
        if (any(candi_models[1, ] == 1)) {
          for (j in seq(n_mo)) {
            varindex=(candi_models[j, ] == 1)
            glmfit =lm(y1 ~ x1[,varindex])
            sigma_k =sig1
            if(any(is.na(glmfit$coef))) {
              lw_num[j+3] = -Inf
            } else {
              dk =sum((y2-cbind(1, x2[,varindex])%*%glmfit$coef)^2)
              lw_num[j+3]=(-(n-n_train))*log(sigma_k)-((sigma_k)^(-2))*dk/2
            }
          }
        } else {
          d1=sum((y2-mean(y1))^2)
          sigma_1=sig1
          lw_num[3+1] =(-(n-n_train))*log(sigma_1)-((sigma_1)^(-2))*d1/2
          for (j in seq(2, n_mo)) {
            varindex=(candi_models[j, ] == 1)
            glmfit =lm(y1 ~ x1[,varindex])
            sigma_k=sig1
            if (any(is.na(glmfit$coef))) {
              lw_num[j+3]=-Inf
            } else {
              dk=sum((y2 - cbind(1, x2[,varindex])%*%glmfit$coef)^2)
              lw_num[j+3] =(-(n-n_train))*log(sigma_k)-((sigma_k)^(-2))*dk/2
            }
          }
        }
        return(lw_num)
      }
      lw_num <- matrix(unlist(mclapply(seq(no_rep), wt_calc)), nrow = no_rep, ncol = 3+n_mo, byrow = TRUE)
      if (prior == TRUE) {
        ck0=-log(1-p0)+log(3)
        ck1=-log(p0)+ck_compute(n_mo, sk, p)
        ck=c(rep(ck0,3),ck1)
        lw_num <- sweep(lw_num, MARGIN = 2, psi*ck, "-")
      }
      lw_num <- sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max), "-")
      w_num <- exp(lw_num)
      weight <- colMeans(w_num/rowSums(w_num))
      weight_se=apply(w_num,2,sd)/sqrt(no_rep)
    }
    
    ###########################method=L1-UARM.rf 
    if (method=='L1-UARM.rf'){
      lscvrf=function(x, y,candi_models, n_train,Rep) {
        p=NCOL(x)
        n=length(y)
        n_mo=NROW(candi_models)
        sk=rowSums(candi_models)
        sigrf2=c()
        for(j in 1:Rep){
          tindex=sample(n, n_train, replace = FALSE)
          x1=x[tindex,];x2=x[-tindex,]
          y1=y[tindex];y2=y[-tindex]
          #####
          sdata_boost=data.frame(y1,x1)
          srandomf=randomForest(y1 ~ ., data=sdata_boost, importance=TRUE,proximity=TRUE)
          spre_rf0=predict(srandomf,newdata=data.frame(x1))
          spre_rf1=predict(srandomf,newdata=data.frame(x2))
          sigrf2[j]=sum(abs(y2-spre_rf1))/n_train
        }
        Sigma2=mean(sigrf2)
        list(Sigma2=Sigma2)
      }
      iscv=lscvrf(x=X,y=y,candi_models,n_train,Rep=20)
      sig2=iscv$Sigma2
      wt_calc <- function(rep_id) {
        lw_num <- matrix(0,nrow=1, 3+n_mo)
        tindex <- sample(n,n_train, replace = FALSE)
        
        x1=x[tindex,];x2=x[-tindex,]
        y1=y[tindex];y2=y[-tindex]
        ######
        sdata_boost=data.frame(y1,x1)
        if (n<200){
          sboost_lm=gbm(y1~.,data=sdata_boost,distribution="gaussian",n.minobsinnode = 1)
          #spre_boost0=predict(sboost_lm,newdata=data.frame(x1), n.trees=n/2,type = "link")
          spre_boost1=predict(sboost_lm,newdata=data.frame(x2), n.trees=n/2,type = "link")
        }else{
          sboost_lm=gbm(y1~.,data=sdata_boost,distribution="gaussian",cv.folds = 5)
          best.iter <- gbm.perf(sboost_lm, method = "cv")
          #spre_boost0=predict(sboost_lm,newdata=data.frame(x1), n.trees=best.iter,type = "link")
          spre_boost1=predict(sboost_lm,newdata=data.frame(x2), n.trees=best.iter,type = "link")
        }
        sigmat1=sig2
        dt1=sum(abs(y2-spre_boost1))
        lw_num[1] <- (-(n-n_train)) *log(sigmat1)-dt1/sigmat1
        ######
        sGboost_lm=glmboost(y=as.numeric(y1),x=x1,center=F)
        spre_Gboost0=predict(sGboost_lm,newdata=x1)
        spre_Gboost1=predict(sGboost_lm,newdata=x2)
        sigmat2=sig2
        dt2=sum(abs(y2-spre_Gboost1))
        lw_num[2] <- (-(n-n_train))*log(sigmat2)-dt2/sigmat2
        ######
        data_rf=data.frame(y1,x1)
        srandomf=randomForest(y1 ~ ., data=sdata_boost)
        spre_rf0=predict(srandomf,newdata=data.frame(x1))
        spre_rf1=predict(srandomf,newdata=data.frame(x2))
        sigmat3=sig2
        dt3=sum(abs(y2-spre_rf1))
        lw_num[3] <- (-(n-n_train))*log(sigmat3)-dt3/sigmat3
        
        if (any(candi_models[1, ] == 1)) {
          for (j in 1:n_mo) {
            varindex <- which(candi_models[j, ] == 1)
            glmfit <- lm(y1 ~ x1[,varindex])
            dk1=sig2
            if(any(is.na(glmfit$coef))) {
              lw_num[j+3] <- -Inf
            } else {
              Dk1 <- sum(abs(y2-cbind(1,x2[,varindex])%*%glmfit$coefficients))
              lw_num[j+3] <- (-(n-n_train))*log(dk1)-Dk1/dk1
            }
          }
        } else {
          Dk0 <- sum(abs(y2 - mean(y1)))
          dk0<-sig2
          lw_num[3+1] <- (-(n-n_train))*log(dk0)-Dk0/dk0
          for (j in 2:n_mo) {
            varindex <- which(candi_models[j, ] == 1)
            glmfit <- lm(y1 ~ x1[,varindex])
            dk1=sig2
            if (any(is.na(glmfit$coef))) {
              lw_num[j+3] <- -Inf
            } else {
              Dk1 <-sum(abs(y2-cbind(1,x2[,varindex])%*%glmfit$coefficients))
              lw_num[j+3] <- (-(n-n_train))*log(dk1)-Dk1/dk1
            }
          }
        }
        return(lw_num)
      }
      lw_num <- matrix(unlist(mclapply(seq(no_rep), wt_calc)), nrow = no_rep, ncol = 3+n_mo, byrow = TRUE)
      if (prior == TRUE) {
        ck0=-log(1-p0)+log(3)
        ck1=-log(p0)+ck_compute(n_mo, sk, p)
        ck=c(rep(ck0,3),ck1)
        lw_num <- sweep(lw_num, MARGIN = 2, psi*ck, "-")
      }
      lw_num <- sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max), "-")
      w_num <- exp(lw_num)
      weight=colMeans(w_num/rowSums(w_num))
      weight_se=apply(w_num,2,sd)/sqrt(no_rep)
    }
  names(weight) = c("GBM", "L2B", "RF", paste("LM", 1:(length(weight) - 3), sep = ""))
  list(weight = weight,weight_se=weight_se,candi_models=candi_models)
  }
}

