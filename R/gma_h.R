###gma_h provide the calculation of model averaging methods for high-dimensional data.
gma_h=function(x, y,factorID=NULL,candidate='H4',candi_models,method='L1-ARM', psi=1,n_train=ceiling(n/2),
               no_rep=20,lambda=log(n),alpha=0.05,prior=TRUE) {
  wbetahat=c()
  # Check the data
  if(all(is.na(x)))
    {
      stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) 
                to eliminate missing data before passing X and y to calculate.")
    }
    else if(is.numeric(x) & is.null(factorID)==TRUE){
      p=NCOL(x)
      n=length(y)
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
      m2=glmnet(x,y,family="gaussian",alpha=1,exclude=which(bhat==0),penalty.factor=adpen)
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
############################method=SAICp
    if (method=='SAICp'){
      ik <- rep(NA, n_mo)
      if (any(candi_models[1, ] == 1)) {
        for (i in 1:n_mo) {
          LSL <- lm(y ~ x[, candi_models[i, ] == 1])
          rss <- crossprod(summary(LSL)$res, summary(LSL)$res)
          ik[i] <- n * log(rss/(n)) + sk[i] * 2
        }
      } else {
        rss <- sum((y - mean(y))^2)
        ik[1] <- n * log(rss/(n)) + sk[1] * 2
        for (i in 2:n_mo) {
          LSL <- lm(y ~ x[, candi_models[i, ] == 1])
          rss <- crossprod(summary(LSL)$res, summary(LSL)$res)
          ik[i] <- n * log(rss/(n)) + sk[i] * 2
        }
      }
      if (prior == TRUE) {
        ck <- ck_compute(n_mo, sk, p)
        ik <- ik + 2 * psi * ck
      }
      ik <- ik - min(ik)
      weight <- exp(-ik/2)/sum(exp(-ik/2))
      wbetahat=betahat%*%weight
    }

    ############################SBICp
    if (method=='SBICp'){
      ik <- rep(NA, n_mo)
      if (any(candi_models[1, ] == 1)) {
        for (i in 1:n_mo) {
          LSL <- lm(y ~ x[, candi_models[i, ] == 1])
          rss <- crossprod(summary(LSL)$res, summary(LSL)$res)
          ik[i] <- n * log(rss/(n)) + sk[i] * log(n)
        }
      } else {
        rss <- sum((y - mean(y))^2)
        ik[1] <- n * log(rss/(n)) + sk[1] * log(n)
        for (i in 2:n_mo) {
          LSL <- lm(y ~ x[, candi_models[i, ] == 1])
          rss <- crossprod(summary(LSL)$res, summary(LSL)$res)
          ik[i] <- n * log(rss/(n)) + sk[i] * log(n)
        }
      }
      if (prior == TRUE) {
        ck <- ck_compute(n_mo, sk, p)
        ik <- ik + 2 * psi * ck
      }
      ik <- ik - min(ik)
      weight <- exp(-ik/2)/sum(exp(-ik/2))
      wbetahat=betahat%*%weight
    }

    ########################ARM
    if (method=='ARM'){
      wt_calc=function(rep_id) {
        lw_num <- matrix(0,nrow=1,n_mo)
        tindex <- sample(n,n_train, replace = FALSE)
          x1=x[tindex,];x2=x[-tindex,]
          y1=y[tindex];y2=y[-tindex]

          if (any(candi_models[1, ] == 1)) {
            for (j in seq(n_mo)) {
              varindex=(candi_models[j, ] == 1)
              glmfit =lm(y1 ~ x1[,varindex])
              sigma_k =summary(glmfit)$sigma
              if(any(is.na(glmfit$coef))) {
                lw_num[j] = -Inf
              } else {
                dk =sum((y2-cbind(1, x2[,varindex])%*%glmfit$coef)^2)
                lw_num[j]=(-(n-n_train))*log(sigma_k)-((sigma_k)^(-2))*dk/2
              }
            }
          } else {
            d1=sum((y2-mean(y1))^2)
            sigma_1=sd(y1)
            lw_num[1] =(-(n-n_train))*log(sigma_1)-((sigma_1)^(-2))*d1/2
            for (j in seq(2, n_mo)) {
              varindex=(candi_models[j, ] == 1)
              glmfit =lm(y1 ~ x1[,varindex])
              sigma_k=summary(glmfit)$sigma
              if (any(is.na(glmfit$coef))) {
                lw_num[j]=-Inf
              } else {
                dk=sum((y2 - cbind(1, x2[,varindex])%*%glmfit$coef)^2)
                lw_num[j] =(-(n-n_train))*log(sigma_k)-((sigma_k)^(-2))*dk/2
              }
            }
          }
        return(lw_num)
      }
      lw_num=matrix(unlist(mclapply(seq(no_rep), wt_calc)),nrow = no_rep,ncol=n_mo, byrow=TRUE)
      if (prior == TRUE) {
        ck=ck_compute(n_mo, sk, p)
        lw_num=sweep(lw_num, MARGIN = 2, psi * ck, "-")
      }
      lw_num=sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max), "-")
      w_num=exp(lw_num)
      weight=colMeans(w_num/rowSums(w_num))
      weight_se=apply(w_num,2,sd)/sqrt(no_rep)
      wbetahat=betahat%*%weight
    }
    
    ########################l1-ARM
    if (method=='L1-ARM'){
      wt_calc <- function(rep_id) {
        lw_num <- matrix(0,nrow=1, n_mo)
        tindex <- sample(n,n_train, replace = FALSE)
          x1=x[tindex,];x2=x[-tindex,]
          y1=y[tindex];y2=y[-tindex]
          ######
          if (any(candi_models[1, ] == 1)) {
            for (j in 1:n_mo) {
              varindex <- which(candi_models[j, ] == 1)
              glmfit <- lm(y1 ~ x1[,varindex])
              dk1=sum(abs(y1-cbind(1,x1[,varindex])%*%glmfit$coefficients))/(n_train)
              if(any(is.na(glmfit$coef))) {
                lw_num[j] <- -Inf
              } else {
                Dk1 <- sum(abs(y2-cbind(1,x2[,varindex])%*%glmfit$coefficients))
                lw_num[j] <- (-(n-n_train))*log(dk1)-Dk1/dk1
              }
            }
          } else {
            Dk0 <- sum(abs(y2 - mean(y1)))
            dk0<-sum(abs(y1-mean(y1)))/(n_train)
            lw_num[1] <- (-(n-n_train))*log(dk0)-Dk0/dk0
            for (j in 2:n_mo) {
              varindex <- which(candi_models[j, ] == 1)
              glmfit <- lm(y1 ~ x1[,varindex])
              dk1=sum(abs(y1-cbind(1,x1[,varindex])%*%glmfit$coefficients))/(n_train)
              if (any(is.na(glmfit$coef))) {
                lw_num[j] <- -Inf
              } else {
                Dk1 <-sum(abs(y2-cbind(1,x2[,varindex])%*%glmfit$coefficients))
                lw_num[j] <- (-(n-n_train))*log(dk1)-Dk1/dk1
              }
            }
          }
        return(lw_num)
      }
      lw_num=matrix(unlist(mclapply(seq(no_rep), wt_calc)), nrow = no_rep, ncol =n_mo, byrow = TRUE)
      if (prior == TRUE) {
        ck = ck_compute(n_mo, sk, p)
        lw_num =sweep(lw_num, MARGIN = 2, psi * ck, "-")
      }
      lw_num=sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max), "-")
      w_num=exp(lw_num)
      weight=colMeans(w_num/rowSums(w_num))
      weight_se=apply(w_num,2,sd)/sqrt(no_rep)
      wbetahat=betahat%*%weight
    }

    ########################l1ARM
    if (method=='PMA'){
      Psigma0=c()
      muhat=matrix(0,n,n_mo)
      if (any(candi_models[1, ] == 1)) {
        for(j in 1:n_mo){
          dd=c(which(candi_models[j,]==1))
          xp=as.matrix(x[,dd])
          np=ncol(xp)
          lmp=lm(y~xp)
          muhat[,j]=cbind(1,xp)%*%lmp$coefficients
          Psigma0[j]=sum((y-cbind(1,xp)%*%lmp$coefficients)^2)/(n-np)
        }
      }else {
        Psigma0[1]=(sd(y))^2
        muhat[,1]=mean(y)
        for(j in 2:n_mo){
          dd=c(which(candi_models[j,]==1))
          xp=as.matrix(x[,dd])
          np=ncol(xp)
          lmp=lm(y~xp)
          muhat[,j]=cbind(1,xp)%*%lmp$coefficients
          Psigma0[j]=sum((y-cbind(1,xp)%*%lmp$coefficients)^2)/(n-np)
        }
      }
      #####################################weight compute
      YY=matrix(y,n,1)%*%rep(1,n_mo)#每一列都是Y
      ehat=YY-muhat
      a1 <- t(ehat) %*% ehat
      if (qr(a1)$rank<ncol(ehat)) {
        a1 <- a1 + diag(n_mo)*1e-10
      }
      psk=Psigma0[which.max(sk)]*sk
      a2<-matrix(as.vector(-lambda*c(psk)),n_mo,1) #嵌套
      a3 <- t(rbind(matrix(1,nrow=1,ncol=n_mo),diag(n_mo),-diag(n_mo)))
      a4 <- rbind(1,matrix(0,nrow=n_mo,ncol=1),matrix(-1,nrow=n_mo,ncol=1))
      w0 <- matrix(1,nrow=n_mo,ncol=1)/n_mo #%权重初???
      QP1 <- solve.QP(a1,a2,a3,a4,1)
      w_pma <- QP1$solution
      w_pma <- w_pma*(w_pma>0)
      w_pma=w_pma/sum(w_pma)
      ########################final estimators
      weight=w_pma
      wbetahat=betahat%*%weight
    }

    ########################MCV
    if (method=='MCV'){
      COR <- rep(0,len=p)
      for(i in 1:p){COR[i] <- cor.test(x[,i],y,"two.sided")$p.value}
      Index <- cbind(1:p,COR)
      Index <- Index[order(Index[,2],decreasing=F),1:2]
      Maxp <- sum(Index[,2]<=alpha)
      if (Maxp==0){
        stop("there are no variables have been chosen under the correlation test,
             it is recommended to increase the alpha before trying again.")
      }
      Q <- seq(0.05,0.40,len=8)*n   #grid points  n*h

      #Model search
      minCV.score <- 10^10
      for(Qind in 1:length(Q)){
        q <- Q[Qind]
        if(trunc(Maxp/q)-abs(Maxp/q)==0){MaxM <- trunc(Maxp/q)}
        if(trunc(Maxp/q)-abs(Maxp/q)!=0){MaxM <- trunc(Maxp/q)+1}
        for(M in 2:MaxM){
          MU <- matrix(0,n,M)
          for(k in 1:M){
            if(k<MaxM){
              USE <- Index[(q*(k-1)+1):(q*k)]
              Zk <- x[,USE]
              nk=dim(as.matrix(Zk))[2]
              Hk <- Zk%*%solve(t(Zk)%*%Zk+diag(0.001,nk,nk))%*%t(Zk)
              Dk <- diag(1/as.vector(diag(diag(1,n)-Hk)))
              Hk <- Dk%*%(Hk-diag(1,n))+diag(1,n)
              MU[,k] <- Hk%*%y
            }
            if(k==MaxM){
              USE <- Index[(q*(k-1)+1):Maxp]
              Zk <- x[,USE]
              nk=dim(as.matrix(Zk))[2]
              Hk <- Zk%*%solve(t(Zk)%*%Zk+diag(0.001,nk,nk))%*%t(Zk)
              Dk <- diag(1/as.vector(diag(diag(1,n)-Hk)))
              Hk <- Dk%*%(Hk-diag(1,n))+diag(1,n)
              MU[,k] <- Hk%*%y
            }

            CV <- function(w){
              sum( (MU%*%w-y)^2 )
            }
            w <- rep(0.1,len=M)
            fit <- optim(w,fn=CV,method="L-BFGS-B",lower=rep(0,len=M),
                         upper=rep(1,len=M))#
            w <- fit$par
            CV.score <- fit$value

            if(CV.score<=minCV.score){
              #minCV.score <- CV.score
              Opt.M <- M
              Opt.w <- w
              Opt.q <- q
            }

          }
        }
      }
      #Construction of the final estimator
      M <- Opt.M
      q <- Opt.q
      wcv<- Opt.w
      if(trunc(Maxp/q)-abs(Maxp/q)==0){MaxM <- trunc(Maxp/q)}
      if(trunc(Maxp/q)-abs(Maxp/q)!=0){MaxM <- trunc(Maxp/q)+1}

      MU <- matrix(0,n,M)
      beta_cv=list(0)
      for(k in 1:M){
        if(k<MaxM){
          USE <- Index[(q*(k-1)+1):(q*k)]
          Zk <- x[,USE]
          nk=dim(as.matrix(Zk))[2]
          Hk <- Zk%*%solve(t(Zk)%*%Zk+diag(0.001,nk,nk))%*%t(Zk)
          MU[,k] <- Hk%*%y
          beta_cv[[k]]=solve(t(Zk)%*%Zk+diag(0.001,nk,nk))%*%t(Zk)%*%y
        }
        if(k==MaxM){
          USE <- Index[(q*(k-1)+1):Maxp]
          Zk <- x[,USE]
          nk=dim(as.matrix(Zk))[2]
          Hk <- Zk%*%solve(t(Zk)%*%Zk+diag(0.001,nk,nk))%*%t(Zk)
          MU[,k] <- Hk%*%y
          beta_cv[[k]]=solve(t(Zk)%*%Zk+diag(0.001,nk,nk))%*%t(Zk)%*%y
        }
      }
       weight=wcv
       M=M;q=q;MaxM=MaxM;Index=Index;Maxp=Maxp
      }
    if (method=='MCV'){
      list(weight=weight,candi_models=list(beta_cv=beta_cv,Index=Index,Maxp=Maxp,MaxM=MaxM,q=q,M=M))
    }else if (method=='PMA'|method=="SAICp"|method=="SBICp"){
      list(weight=weight,candi_models=candi_models,betahat=betahat,wbetahat=wbetahat)
    }else{
      list(weight=weight,weight_se=weight_se,candi_models=candi_models,betahat=betahat,wbetahat=wbetahat)
    }
  } else {
    n=length(y)
    candidate='CH3'
    #####factorID
    if(is.numeric(factorID)){
      cat=factorID
    } else {
    cat=which(colnames(x)==factorID)
    }
    
    lengthp=c()
    for(i in 1:length(cat)){
      lengthp[i]=length(summary(x[,cat[i]]))-1
    }
    leng=sum(lengthp) ##the sum of categorical

    for(i in 1:length(cat)){
      h=summary(as.factor(x[,cat[i]]))
      if (min(h)<=((n/length(h))/2))
      {
        stop("The level of categorical data appears to be extremely unbalanced, 
             and the number of observations at some levels is too small.it is recommended to 
             try again after deal with categorical data reasonablely.")
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
    x0=Xm[,-1] ###DUMMY MATRIX
    index_group0=list()
    index_group=c()
    for (i in 1:length(XX)){
      index_group0[[i]]=rep(i,dim(as.matrix(XX[[i]]))[2])
      index_group=c(index_group,index_group0[[i]])
    }
    
    if(candidate=='CH3'){
      m1=grpreg(x0,y,group=index_group,penalty="grLasso")
      m2=grpreg(x0,y,group=index_group,penalty="grMCP")
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
      candi_models=candidate_models[c(de),-1]##delete intercept 
    }
    x=x0  ##the version contains dummy variables
    p=NCOL(x)
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
  ##############################SAICp
  if (method=='SAICp'){
    ik <- rep(NA, n_mo)
    if (any(candi_models[1, ] == 1)) {
      for (i in 1:n_mo) {
        LSL <- lm(y ~ x[, candi_models[i, ] == 1])
        rss <- crossprod(summary(LSL)$res, summary(LSL)$res)
        ik[i] <- n * log(rss/(n)) + sk[i] * 2
      }
    } else {
      rss <- sum((y - mean(y))^2)
      ik[1] <- n * log(rss/(n)) + sk[1] * 2
      for (i in 2:n_mo) {
        LSL <- lm(y ~ x[, candi_models[i, ] == 1])
        rss <- crossprod(summary(LSL)$res, summary(LSL)$res)
        ik[i] <- n * log(rss/(n)) + sk[i] * 2
      }
    }
    if (prior == TRUE) {
      ck <- ck_compute(n_mo, sk, p)
      ik <- ik + 2 * psi * ck
    }
    ik <- ik - min(ik)
    weight <- exp(-ik/2)/sum(exp(-ik/2))
    wbetahat=betahat%*%weight
  }

  ############################SBICp
  if (method=='SBICp'){
    ik <- rep(NA, n_mo)
    if (any(candi_models[1, ] == 1)) {
      for (i in 1:n_mo) {
        LSL <- lm(y ~ x[, candi_models[i, ] == 1])
        rss <- crossprod(summary(LSL)$res, summary(LSL)$res)
        ik[i] <- n * log(rss/(n)) + sk[i] * log(n)
      }
    } else {
      rss <- sum((y - mean(y))^2)
      ik[1] <- n * log(rss/(n)) + sk[1] * log(n)
      for (i in 2:n_mo) {
        LSL <- lm(y ~ x[, candi_models[i, ] == 1])
        rss <- crossprod(summary(LSL)$res, summary(LSL)$res)
        ik[i] <- n * log(rss/(n)) + sk[i] * log(n)
      }
    }
    if (prior == TRUE) {
      ck <- ck_compute(n_mo, sk, p)
      ik <- ik + 2 * psi * ck
    }
    ik <- ik - min(ik)
    weight <- exp(-ik/2)/sum(exp(-ik/2))
    wbetahat=betahat%*%weight
  }

    ########################ARM
    if (method=='ARM'){
      wt_calc=function(rep_id) {
        lw_num <- matrix(0,nrow=1,n_mo)
        tindex <- sample(n,n_train, replace = FALSE)
        x1=x[tindex,];x2=x[-tindex,]
        y1=y[tindex];y2=y[-tindex]
        
        if (any(candi_models[1, ] == 1)) {
          for (j in seq(n_mo)) {
            varindex=(candi_models[j, ] == 1)
            glmfit =lm(y1 ~ x1[,varindex])
            sigma_k =summary(glmfit)$sigma
            if(any(is.na(glmfit$coef))) {
              lw_num[j] = -Inf
            } else {
              dk =sum((y2-cbind(1, x2[,varindex])%*%glmfit$coef)^2)
              lw_num[j]=(-(n-n_train))*log(sigma_k)-((sigma_k)^(-2))*dk/2
            }
          }
        } else {
          d1=sum((y2-mean(y1))^2)
          sigma_1=sd(y1)
          lw_num[1] =(-(n-n_train))*log(sigma_1)-((sigma_1)^(-2))*d1/2
          for (j in seq(2, n_mo)) {
            varindex=(candi_models[j, ] == 1)
            glmfit =lm(y1 ~ x1[,varindex])
            sigma_k=summary(glmfit)$sigma
            if (any(is.na(glmfit$coef))) {
              lw_num[j]=-Inf
            } else {
              dk=sum((y2 - cbind(1, x2[,varindex])%*%glmfit$coef)^2)
              lw_num[j] =(-(n-n_train))*log(sigma_k)-((sigma_k)^(-2))*dk/2
            }
          }
        }
        return(lw_num)
      }
      lw_num=matrix(unlist(mclapply(seq(no_rep), wt_calc)),nrow = no_rep,ncol=n_mo, byrow=TRUE)
      if (prior == TRUE) {
        ck=ck_compute(n_mo, sk, p)
        lw_num=sweep(lw_num, MARGIN = 2, psi * ck, "-")
      }
      lw_num=sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max), "-")
      w_num=exp(lw_num)
      weight=colMeans(w_num/rowSums(w_num))
      weight_se=apply(w_num,2,sd)/sqrt(no_rep)
      wbetahat=betahat%*%weight
    }
    
    ########################l1-ARM
    if (method=='L1-ARM'){
      wt_calc <- function(rep_id) {
        lw_num <- matrix(0,nrow=1, n_mo)
        tindex <- sample(n,n_train, replace = FALSE)
        x1=x[tindex,];x2=x[-tindex,]
        y1=y[tindex];y2=y[-tindex]
        ######
        if (any(candi_models[1, ] == 1)) {
          for (j in 1:n_mo) {
            varindex <- which(candi_models[j, ] == 1)
            glmfit <- lm(y1 ~ x1[,varindex])
            dk1=sum(abs(y1-cbind(1,x1[,varindex])%*%glmfit$coefficients))/(n_train)
            if(any(is.na(glmfit$coef))) {
              lw_num[j] <- -Inf
            } else {
              Dk1 <- sum(abs(y2-cbind(1,x2[,varindex])%*%glmfit$coefficients))
              lw_num[j] <- (-(n-n_train))*log(dk1)-Dk1/dk1
            }
          }
        } else {
          Dk0 <- sum(abs(y2 - mean(y1)))
          dk0<-sum(abs(y1-mean(y1)))/(n_train)
          lw_num[1] <- (-(n-n_train))*log(dk0)-Dk0/dk0
          for (j in 2:n_mo) {
            varindex <- which(candi_models[j, ] == 1)
            glmfit <- lm(y1 ~ x1[,varindex])
            dk1=sum(abs(y1-cbind(1,x1[,varindex])%*%glmfit$coefficients))/(n_train)
            if (any(is.na(glmfit$coef))) {
              lw_num[j] <- -Inf
            } else {
              Dk1 <-sum(abs(y2-cbind(1,x2[,varindex])%*%glmfit$coefficients))
              lw_num[j] <- (-(n-n_train))*log(dk1)-Dk1/dk1
            }
          }
        }
        return(lw_num)
      }
      lw_num=matrix(unlist(mclapply(seq(no_rep), wt_calc)), nrow = no_rep, ncol =n_mo, byrow = TRUE)
      if (prior == TRUE) {
        ck = ck_compute(n_mo, sk, p)
        lw_num =sweep(lw_num, MARGIN = 2, psi * ck, "-")
      }
      lw_num=sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max), "-")
      w_num=exp(lw_num)
      weight=colMeans(w_num/rowSums(w_num))
      weight_se=apply(w_num,2,sd)/sqrt(no_rep)
      wbetahat=betahat%*%weight
    }

  ########################PMA
  if (method=='PMA'){
    Psigma0=c()
    muhat=matrix(0,n,n_mo)
    if (any(candi_models[1, ] == 1)) {
      for(j in 1:n_mo){
        dd=c(which(candi_models[j,]==1))
        xp=as.matrix(x[,dd])
        np=ncol(xp)
        lmp=lm(y~xp)
        muhat[,j]=cbind(1,xp)%*%lmp$coefficients
        Psigma0[j]=sum((y-cbind(1,xp)%*%lmp$coefficients)^2)/(n-np)
      }
    }else {
      Psigma0[1]=(sd(y))^2
      muhat[,1]=mean(y)
      for(j in 2:n_mo){
        dd=c(which(candi_models[j,]==1))
        xp=as.matrix(x[,dd])
        np=ncol(xp)
        lmp=lm(y~xp)
        muhat[,j]=cbind(1,xp)%*%lmp$coefficients
        Psigma0[j]=sum((y-cbind(1,xp)%*%lmp$coefficients)^2)/(n-np)
      }
    }
    #####################################weight compute
    YY=matrix(y,n,1)%*%rep(1,n_mo)#每一列都是Y
    ehat=YY-muhat
    a1 <- t(ehat) %*% ehat
    if (qr(a1)$rank<ncol(ehat)) {
      a1 <- a1 + diag(n_mo)*1e-10
    }
    psk=Psigma0[which.max(sk)]*sk
    a2<-matrix(as.vector(-lambda*c(psk)),n_mo,1) #嵌套
    a3 <- t(rbind(matrix(1,nrow=1,ncol=n_mo),diag(n_mo),-diag(n_mo)))
    a4 <- rbind(1,matrix(0,nrow=n_mo,ncol=1),matrix(-1,nrow=n_mo,ncol=1))
    w0 <- matrix(1,nrow=n_mo,ncol=1)/n_mo #%权重初???
    QP1 <- solve.QP(a1,a2,a3,a4,1)
    w_pma <- QP1$solution
    w_pma <- w_pma*(w_pma>0)
    w_pma=w_pma/sum(w_pma)
    ########################final estimators
    weight=w_pma
    wbetahat=betahat%*%weight
  }
  ########################
if (method=='PMA'|method=="SAICp"|method=="SBICp"){
    list(weight=weight,candi_models=candi_models,betahat=betahat,wbetahat=wbetahat)
  }else{
    list(weight=weight,weight_se=weight_se,candi_models=candi_models,betahat=betahat,wbetahat=wbetahat)
  }
 }
}

