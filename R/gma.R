cdm <- function(candi_models,x)
{
  p <- ncol(x); 
  s=as.matrix(candi_models)
  if ((nrow(s)==1) && (ncol(s)==1)){
    if (candi_models == 1){
      s <- matrix(1,nrow=p,ncol=p)
      s[upper.tri(s)] <- 0
      s <- rbind(0,s)
    }
    if (candi_models == 2){
      s <- matrix(0,nrow=2^p,ncol=p)
      s0 <- matrix(c(1,rep(0,p-1)),1,p)
      s1 <- matrix(c(rep(0,p)),1,p)
      for (i in 2:2^p){
        s1 <- s0 + s1
        for (j in 1:p){
          if (s1[1,j] == 2){
            s1[1,j+1] <- s1[1,j+1]+1
            s1[1,j] <- 0
          }
        }
        s[i,] <- s1
      }
    }
  }
  candi_models=s
  return(s)
}

wFIC.normal <- function(Y, X, Z, XZeval=cbind(X,Z), weightsvec=rep(1,nrow(XZeval)), 
                        subsetsmatrix, partial.delta.theta, partial.delta.gamma)
{
  weightsvec <- as.vector(weightsvec)
  
  qq      = ncol(as.matrix(Z))
  pp      = ncol(as.matrix(X)) 
  nn      = nrow(as.matrix(X))
  
  XZeval  = matrix(XZeval,nrow = nn)
  X.eval  = XZeval[,1:pp]; X.eval<-matrix(X.eval,nrow = nn)
  Z.eval  = XZeval[,-(1:pp)]; Z.eval<-matrix(Z.eval,nrow = nn)
  
  full.fit =  lm(Y~X+Z-1)  # exclude intercept since already in the X-part
  
  sigma.hat.sq<-as.vector(t(full.fit$residuals)%*%full.fit$residuals/(nn-pp-qq))
  
  J <- matrix(0,nrow<-(pp+qq),ncol<-(pp+qq)) 
  #J[1,1]<-2
  J[(1:(pp)),(1:(pp))] <- t(X)%*%X
  J[(pp+1):(pp+qq),(1:(pp))] <- t(Z)%*%X
  J[(1:(pp)),(pp+1):(pp+qq)] <- t(X)%*%Z
  J[(pp+1):(pp+qq),(pp+1):(pp+qq)] <- t(Z)%*%Z
  
  J <- J/(nn*sigma.hat.sq)
  
  J11 <- J[(pp+1):(pp+qq),(pp+1):(pp+qq)] 
  J00 <- J[1:(pp),1:(pp)]
  J10 <- J[(pp+1):(pp+qq),1:(pp)]
  
  invJ <- solve(J)
  K <- invJ[-(1:(pp)),-(1:(pp))]
  
  ## estimate (beta,gamma) in biggest model
  betagamma.hat <- full.fit$coefficients
  deltahat <- nn^{1/2}*betagamma.hat[-(1:pp)]
  
  n.eval <- nrow(XZeval)
  
  wFIC.model<-rep(NA,nrow(subsetsmatrix))
  
  XZeval <- as.matrix(XZeval)
  X.eval <- XZeval[,1:pp]; X.eval<-matrix(X.eval,ncol=pp)
  Z.eval <- XZeval[,-(1:pp)]; Z.eval<-matrix(Z.eval,ncol=qq)
  
  #CAR create the A matrix, as defined on p. 507 of Claeskens and Hjort (5008a)      
  omega.omega.t.sum<-matrix(0,nrow<-qq,ncol<-qq)
  weights.sum<-0
  
  for (k in 1:n.eval) {
    if (weightsvec[k] > 0) {
      this.x <- X.eval[k,]
      this.z <- Z.eval[k,]
      omega <- J10%*%solve(J00)%*%partial.delta.theta(this.x,this.z) - 
        partial.delta.gamma(this.x,this.z)
      omega.omega.t.sum <- omega.omega.t.sum + omega%*%t(omega) * weightsvec[k]
      weights.sum <- weights.sum+weightsvec[k]
    }
  }
  
  A <- omega.omega.t.sum/weights.sum 
  
  Iq <- diag(rep(1,qq))
  
  for (j in 1:nrow(subsetsmatrix)) {
    indicS <- subsetsmatrix[j,]*(1:qq)
    variables <- subsetsmatrix[j,]
    if (sum(variables) > 0)
    {
      qq <- length(variables)
      indic <- variables*(1:qq)
      
      Id <-  diag (rep(1,qq))
      pi.S <- matrix(Id[indic,],ncol=qq)
      
      K.S <- solve(pi.S %*% solve(K) %*% t(pi.S))
      G.S <- t(pi.S)%*%K.S%*%pi.S%*%solve(K)   ###p.506 of Claeskens and Hjort (5008a)
      
      sq.bias.estimate <- sum(
        diag((Iq-G.S)%*%(deltahat%*%t(deltahat)-K)%*%t((Iq-G.S))%*%A))
      wFIC.model[j] <- max(sq.bias.estimate,0) + sum(diag(A%*%t(pi.S)%*%K.S%*%pi.S))
    } else {
      wFIC.model[j] <- max(sum(diag((deltahat%*%t(deltahat)-K)%*%A)), 0)
    }
    wFIC=  exp(-0.5*wFIC.model/c(t(omega)%*%K%*%omega))/sum(exp(-0.5*wFIC.model/c(t(omega)%*%K%*%omega)))
  } 
  # return(wFIC.model)
  return(wFIC) 
}

ARMS <-function(x, y,method = "L1-ARM",candi_models=candi_models, n_train=ceiling(n/2), no_rep=50) 
  {
    x            = as.data.frame(x)   ;    p       = ncol(x)            ;   n  = length(y)
    candi_models = cdm(candi_models,x);    n_mo    = nrow(candi_models) ;   sk = rowSums(candi_models)

    wt_calc <- function(rep_id)                     {
    lwnmutr  = lwnmute = rep(NA, n_mo); 
    tridx    = sample(n, n_train, replace = FALSE)
      for (j in seq(n_mo))                      {
        varindex     = (candi_models[j, ] == 1)
        trdat        = data.frame(y[tridx]  , x[tridx, varindex] )  ;  names(trdat) = c( "y",names(x)[varindex] )   
        tedat        = data.frame(y[-tridx] , x[-tridx, varindex])  ;  names(tedat) = c( "y",names(x)[varindex] )     
        glmtr        = glm(y ~., data = trdat)                      ;  #glmte        = glm(y ~., data = tedat)
         if(  any( is.na(glmtr$coef) )  )  {
         lwnmutr[j]   = lwnumte[j] = -Inf
                                           } 
  else {
         if(method=="L1-ARM") { 
         sigtr_k      = mean(abs(glmtr$res))                   
         #sigte_k      = mean(abs(glmte$res))
         trdk         = sum( abs(tedat$y - predict(glmtr,newdata = tedat) ) )
         #tedk         = sum( abs(trdat$y - predict(glmte,newdata = trdat) ) )
         lwnmutr[j]   = - n_train    * log(sigtr_k) - ((sigtr_k)^(-1)) * trdk
         #lwnmute[j]   = -(n-n_train) * log(sigte_k) - ((sigte_k)^(-1)) * tedk
                              }
         if(method=="ARM")    {
         sigtr_k      = sqrt(sum((glmtr$res)^2)/(n_train-sk[j]-1))
         #sigte_k      = sqrt(sum((glmte$res)^2)/(n-n_train-sk[j]-1))
         trdk         = sum( (tedat$y - predict(glmtr, newdata = tedat) )^2 )
         #tedk         = sum( (trdat$y - predict(glmte, newdata = trdat) )^2 )
         lwnmutr[j]   = - n_train    * log(sigtr_k) - ((sigtr_k)^(-2)) * trdk/2
         #lwnmute[j]   = -(n-n_train) * log(sigte_k) - ((sigte_k)^(-2)) * tedk/2
                              }              
      }
                                                  }
      lw_num          = lwnmutr#(lwnmutr+lwnmute)/2                                        
      return(lw_num)
                                                      }
    lw_num    = matrix(unlist(mclapply(seq(no_rep), wt_calc)), nrow = no_rep, ncol = n_mo, byrow = TRUE)
    lw_num    = sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max), "-")
    w_num     = exp(lw_num)
    weight    = colMeans(w_num/rowSums(w_num))
    weight_se = (apply(w_num, 2, sd))/sqrt(no_rep)
    
    list(weight = weight,weight_se = weight_se )
  }



gma<-function (x, y,factorID=NULL, method = "L1-ARM", candi_models = 2, n_train=ceiling(n/2), no_rep=50) 
{
  x = as.data.frame(x) 
  
  if(is.null(factorID)==FALSE)
  {for (i in 1:length(factorID))
  {
    x[,(factorID[i])] = factor(x[,factorID[i]])
  }
  }
  
  p =  ncol(x);     n =  length(y)     ;    candi_models = cdm(candi_models,x)
  n_mo = nrow(candi_models)            ;    ik = sk = rowSums(candi_models)
  ee = eej =  matrix(0, nrow = nrow(x), ncol = n_mo)
  bname=names(glm( y ~., data=x)$coefficients)
  betahat=matrix(0,length(bname),n_mo) ;    rownames(betahat)=bname
  
  for (i in 1:n_mo) {
      dati        = data.frame(y,x[, candi_models[i, ] == 1])
      names(dati) = c("y",names(x)[candi_models[i, ] == 1])
      LSL         = glm(y~., data=dati)
      sk[i]       = length((LSL)$coefficients)
      betahat[names((LSL)$coefficients),i] = (LSL)$coefficients
      ee[, i]     = residuals(LSL)
      hi          = hatvalues(LSL)
      eej[, i]    = ee[, i] * ( 1 / (1 - hi) )
        if(method == "SBIC") ik[i] <- BIC(LSL) #n * log(rss/n) + sk[i] * log(n)
        if(method == "SAIC") ik[i] <- AIC(LSL) #n * log(rss/n) + sk[i] * 2 
                    }
  if (method == "SAIC" | method == "SBIC") {
    weight    = exp(-0.5 * (ik - min(ik)))/sum(exp(-0.5 *  (ik - min(ik))))
    weight    = matrix(weight,n_mo,1)
    wbetahat  = betahat %*% weight
    return(list(weight = weight, wbetahat=wbetahat, betahat=betahat, candi_models=candi_models))
                                           }
  if (method == "MMA" | method == "JMA") {
    if (method == "MMA") {
      a1        = t(ee) %*% ee
      mK        = which.max((sk))
      sighat    = (t(ee[, mK]) %*% ee[, mK])/(n - candi_models[mK])
      a2        = matrix(c(-sighat) * sk)
                         }
    if (method == "JMA") {
      a1        = t(eej) %*% eej
      a2        = matrix(0, nrow = n_mo, ncol = 1)                 
                         }
    if (qr(a1)$rank < ncol(ee)) 
      a1 = a1 + diag(n_mo) * 1e-10
    if (qr(a1)$rank < ncol(ee)) 
      a1 = a1 + diag(n_mo) * 1e-08
    if (qr(a1)$rank < ncol(ee)) 
    a1       = a1 + diag(n_mo) * 1e-06
    a3       = t(rbind(matrix(1, nrow = 1, ncol = n_mo), diag(n_mo),-diag(n_mo)))
    a4       = rbind(1, matrix(0, nrow = n_mo, ncol = 1), matrix(-1, nrow = n_mo, ncol = 1))
    w0       = matrix(1, nrow = n_mo, ncol = 1)/n_mo
    QP       = solve.QP(a1, a2, a3, a4, 1)
    w        = QP$solution
    w        = as.matrix(w)
    w        = w * (w > 0)
    weight   = w/sum(w0)
    weight    = matrix(weight,n_mo,1)
    wbetahat = betahat %*% weight
    return(list(weight = weight, wbetahat=wbetahat, betahat=betahat, candi_models=candi_models))
                                          }
  if (method == "L1-ARM" | method == "ARM")    {
    if(is.null(factorID)==FALSE){
    for(i in 1:length(factorID)){
    h=summary(x[, factorID[i]])
    if (min(h)<=(n/length(h)/2))
      {
        stop("The level of categorical data appears to be extremely unbalanced, and the number of observations at some levels is too small. It is recommended to try again after deal with categorical data reasonablely.")
      }
                               }
                               }
  if (method == "L1-ARM") {
    ar = ARMS(x, y,method = "L1-ARM",candi_models=candi_models, n_train=ceiling(n/2), no_rep=50)
                          }
  if (method == "ARM") {
    ar = ARMS(x, y,method = "ARM",candi_models=candi_models, n_train=ceiling(n/2), no_rep=50)
                       }
    weight      = ar$weight
    weight      = matrix(weight,n_mo,1)
    wbetahat    = betahat %*% weight
    weight_se   = ar$weight_se
    return(list(weight = weight, weight_se = weight_se, wbetahat=wbetahat, betahat=betahat, candi_models=candi_models))
                                              }
  if (method == "SFIC") {
  wf = wFIC.normal(Y = y, X = matrix(1, nrow = length(y)), 
                   Z = as.matrix(x), subsetsmatrix = candi_models, partial.delta.theta = function(x, z) 
                     {c(1)}, partial.delta.gamma = function(x, z) { z})
  weight    = wf
  weight    = matrix(weight,n_mo,1)
  wbetahat  = betahat %*% weight
  return(list(weight = weight, wbetahat=wbetahat, betahat=betahat, candi_models=candi_models))
                        }
}


