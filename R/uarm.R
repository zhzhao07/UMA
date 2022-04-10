ck_compute <- function(n_mo, sk, p) {
  ck = 2*log(sk + 2) + ifelse(sk > 0, sk * log(exp(1)*p/sk), 0)
  return(ck)
}

uarm <- function(x, y, factorID = NULL, candi_models=2, n_train = ceiling(n/2), no_rep = 50, psi = 0.1, 
              method ="L1-UARM", prior = TRUE, p0 = 0.5)
{
  x            = as.data.frame(x)     ;  p <- NCOL(x)                   ;   n  = length(y)
  candi_models = cdm(candi_models,x)  ;  n_mo <- nrow(candi_models)     ;   sk = rowSums(candi_models)
  dat          = data.frame(y, x)     ;  names(dat) = c( "y",names(x)) 
  if(is.null(factorID)==FALSE)
    {
    for (i in 1:length(factorID))
    {
    x[,(factorID[i])] = factor(x[,factorID[i]])
    h=summary(x[, factorID[i]])
    if (min(h)<=(nrow(x)/length(h)/2))
    {
      stop("The level of categorical data appears to be extremely unbalanced, and the number of observations at some levels is too small. It is recommended to try again after deal with categorical data reasonablely.")
    }
    }
    }
  for (i in 1:p)            {
    if (all(is.na(x[, i]))) {
      stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to uarm.")
                            }
                            }

      lscvrf  <- function(x, y, n_train, Rep) {
    sigrf1=sigrf2=c()
    for(j in 1:Rep){
      tindex    = sample(n, n_train, replace = FALSE)
      trdat     = dat[tindex,]   ;    tedat       = dat[-tindex,]
      srf       = randomForest(y~ ., data = trdat)
      spre_rf   = predict(srf, newdata = tedat)
      sigrf1[j] = sum( abs(tedat$y-spre_rf) )/n_train
      sigrf2[j] = sqrt(sum((tedat$y-spre_rf)^2)/n_train)
                  }
    Sigma1   = mean(sigrf1) ;  Sigma2 = mean(sigrf2)
    list(Sigma1 = Sigma1, Sigma2 = Sigma2)
                                       }
      iscv     = lscvrf(x=x,y=y,n_train = ceiling(n/2),Rep=50)
      Sig1     = iscv$Sigma1 ;    Sig2  = iscv$Sigma2
      if (method == "L1-UARM") {            
      wt_calc <- function(rep_id)                     {
        lwnmutr  = lwnmute = rep(NA, n_mo + 5); 
        tridx    = sample(n, n_train, replace = FALSE)
        for (j in seq(n_mo))                      {
          varindex     = (candi_models[j, ] == 1)
          trdat        = data.frame(y[tridx]  , x[tridx, varindex] );  names(trdat) = c( "y",names(x)[varindex] )   
          tedat        = data.frame(y[-tridx] , x[-tridx, varindex]);  names(tedat) = c( "y",names(x)[varindex] )     
          glmtr        = glm(y ~., data = trdat)                    ;  #glmte        = glm(y ~., data = tedat)
          sigtr_k      = mean(abs(glmtr$res))                       ;  #sigte_k      = mean(abs(glmte$res))
          trdk         = sum( abs(tedat$y - predict(glmtr,newdata = tedat) ) )
          #tedk         = sum( abs(trdat$y - predict(glmte,newdata = trdat) ) )
          lwnmutr[j+5] = - n_train    * log(sigtr_k) - ((sigtr_k)^(-1)) * trdk
          #lwnmute[j+5] = -(n-n_train) * log(sigte_k) - ((sigte_k)^(-1)) * tedk
                                                    }
          tradat        = dat[tridx,];       teadat        = dat[-tridx,]     
          
          #GMB
          gbmglmtr      = gbm(y ~ ., data = tradat, distribution = "gaussian", n.minobsinnode = 0)
          #gbmglmte      = gbm(y ~ ., data = teadat, distribution = "gaussian", n.minobsinnode = 0)
          trdkgbm       = sum (abs (teadat$y - predict(gbmglmtr,  newdata = teadat, n.trees = min(100,nrow(teadat)/2)) ) )
          #tedkgbm       = sum (abs (tradat$y - predict(gbmglmte,  newdata = tradat, n.trees = min(100,nrow(tradat)/2)) ) )
          lwnmutr[1]    = -n_train     * log(Sig1) - trdkgbm/Sig1
          #lwnmute[1]    = -(n-n_train) * log(Sig1) - tedkgbm/Sig1
          #L2B
          l2bglmtr      = glmboost(y ~ .,data= tradat, center = FALSE)
          #l2bglmte      = glmboost(y ~ .,data= teadat, center = FALSE)
          trdkl2b       = sum (abs (teadat$y - predict(l2bglmtr, newdata = teadat) ) )
          #tedkl2b       = sum (abs (tradat$y - predict(l2bglmte, newdata = tradat) ) )
          lwnmutr[2]    = -n_train     * log(Sig1) - trdkl2b/Sig1
          #lwnmute[2]    = -(n-n_train) * log(Sig1) - tedkl2b/Sig1
          #rf
          rfglmtr       = randomForest(y ~ .,data= tradat)
          #rfglmte       = randomForest(y ~ .,data= teadat)
          trdkrf        = sum (abs (teadat$y - predict(rfglmtr, newdata = teadat) ) )
          #tedkrf        = sum (abs (tradat$y - predict(rfglmte, newdata = tradat) ) )
          lwnmutr[3]    = -n_train     * log(Sig1) - trdkrf/Sig1
          #lwnmute[3]    = -(n-n_train) * log(Sig1) - tedkrf/Sig1
          #bag
          bagglmtr      = bagging(y ~ .,data= tradat, coob = TRUE)
          #bagglmte      = bagging(y ~ .,data= teadat, coob = TRUE)
          trdkbag       = sum (abs (teadat$y - predict(bagglmtr, newdata = teadat) ) )
          #tedkbag       = sum (abs (tradat$y - predict(bagglmte, newdata = tradat) ) )
          lwnmutr[4]    = -n_train     * log(Sig1) - trdkbag/Sig1
          #lwnmute[4]    = -(n-n_train) * log(Sig1) - tedkbag/Sig1
          #bart
          bart_tr       = bart(tradat[,-1], as.double(tradat[,1]), teadat[,-1], ndpost = 500, verbose = F)$yhat.test.mean
          #bart_te       = bart(teadat[,-1], as.double(teadat[,1]), tradat[,-1], ndpost = 500, verbose = F)$yhat.test.mean
          trdkbart      = sum (abs (teadat$y - bart_tr ) )
          #tedkbart      = sum (abs (tradat$y - bart_te ) )
          lwnmutr[5]    = -n_train     * log(Sig1) - trdkbart/Sig1
          #lwnmute[5]    = -(n-n_train) * log(Sig1) - tedkbart/Sig1
          lw_num        = lwnmutr# (lwnmutr+lwnmute)/2 
          return (lw_num)
      } 
      }
    
    
      
      if (method == "UARM") {            
      wt_calc <- function(rep_id)                     {
          lwnmutr  = lwnmute = rep(0, n_mo + 5); 
          tridx    = sample(n, n_train, replace = FALSE)
          for (j in seq(n_mo))                      {
            varindex     = (candi_models[j, ] == 1)
            trdat        = data.frame(y[tridx]  , x[tridx, varindex] )   ;      names(trdat) = c( "y",names(x)[varindex] )   
            tedat        = data.frame(y[-tridx] , x[-tridx, varindex])   ;      names(tedat) = c( "y",names(x)[varindex] )     
            glmtr        = glm(y ~., data = trdat)                       ;    #  glmte        = glm(y ~., data = tedat)
            sigtr_k      = sqrt(sum((glmtr$res)^2)/(n_train-sk[j]-1))    ;    #  sigte_k      = sqrt(sum((glmte$res)^2)/(n-n_train-sk[j]-1))
            trdk         = sum( (tedat$y - predict(glmtr,newdata = tedat) )^2 )
            #tedk         = sum( (trdat$y - predict(glmte,newdata = trdat) )^2 )
            lwnmutr[j+5] = - n_train    * log(sigtr_k) - (sigtr_k^(-2)) * trdk/2
            #lwnmute[j+5] = -(n-n_train) * log(sigte_k) - (sigte_k^(-2)) * tedk/2
                                                     }
            tradat        = dat[tridx,];       teadat        = dat[-tridx,]     
            #GMB
            gbmglmtr      = gbm(y ~ ., data = tradat, distribution = "gaussian", n.minobsinnode = 0)
           # gbmglmte      = gbm(y ~ ., data = teadat, distribution = "gaussian", n.minobsinnode = 0)
            trdkgbm       = sum ( (teadat$y - predict(gbmglmtr, n.trees = min(100,nrow(teadat)/2), newdata = teadat) )^2 )
          #  tedkgbm       = sum ( (tradat$y - predict(gbmglmte, n.trees = min(100,nrow(tradat)/2), newdata = tradat) )^2 )
            lwnmutr[1]    = -n_train     * log(Sig2) - trdkgbm/(2*Sig2^2)
          #  lwnmute[1]    = -(n-n_train) * log(Sig2) - tedkgbm/(2*Sig2^2)
            #L2B
            l2bglmtr      = glmboost(y ~ .,data= tradat, center = FALSE)
          #  l2bglmte      = glmboost(y ~ .,data= teadat, center = FALSE)
            trdkl2b       = sum ( (teadat$y - predict(l2bglmtr, newdata = teadat) )^2 )
          #  tedkl2b       = sum ( (tradat$y - predict(l2bglmte, newdata = tradat) )^2 )
            lwnmutr[2]    = -n_train     * log(Sig2) - trdkl2b/(2*Sig2^2)
          #  lwnmute[2]    = -(n-n_train) * log(Sig2) - tedkl2b/(2*Sig2^2)
            #rf
            rfglmtr       = randomForest(y ~ .,data= tradat)
          #  rfglmte       = randomForest(y ~ .,data= teadat)
            trdkrf        = sum ( (teadat$y - predict(rfglmtr, newdata = teadat) )^2 )
          #  tedkrf        = sum ( (tradat$y - predict(rfglmte, newdata = tradat) )^2 )
            lwnmutr[3]    = -n_train     * log(Sig2) - trdkrf/(2*Sig2^2)
           # lwnmute[3]    = -(n-n_train) * log(Sig2) - tedkrf/(2*Sig2^2)
            #bag
            bagglmtr      = bagging(y ~ .,data= tradat, coob = TRUE)
           # bagglmte      = bagging(y ~ .,data= teadat, coob = TRUE)
            trdkbag       = sum ( (teadat$y - predict(bagglmtr, newdata = teadat) )^2 )
           # tedkbag       = sum ( (tradat$y - predict(bagglmte, newdata = tradat) )^2 )
            lwnmutr[4]    = -n_train     * log(Sig2) - trdkbag/(2*Sig2^2)
           # lwnmute[4]    = -(n-n_train) * log(Sig2) - tedkbag/(2*Sig2^2)
            #bart
            bart_tr       = bart(tradat[,-1], as.double(tradat[,1]), teadat[,-1], ndpost = 500, verbose = F)$yhat.test.mean
           # bart_te       = bart(teadat[,-1], as.double(teadat[,1]), tradat[,-1], ndpost = 500, verbose = F)$yhat.test.mean
            trdkbart      = sum ( (teadat$y - bart_tr )^2 )
           # tedkbart      = sum ( (tradat$y - bart_te )^2 )
            lwnmutr[5]    = -n_train     * log(Sig2) - trdkbart/(2*Sig2^2)
           # lwnmute[5]    = -(n-n_train) * log(Sig2) - tedkbart/(2*Sig2^2)
            
          
            lw_num          = lwnmutr# (lwnmutr+lwnmute)/2 
          return(lw_num)
                                                       } 
                            }
      lw_num  = matrix(unlist(mclapply(seq(no_rep), wt_calc)), nrow = no_rep, ncol = 5 + n_mo, byrow = TRUE)
      if (prior == TRUE)
                           {
      ck0     = log(5) - log(1 - p0)
      ck1     = -log(p0) + ck_compute(n_mo, sk, p)
      ck      = c(rep(ck0, 5), ck1)
      lw_num  = sweep(lw_num, MARGIN = 2, psi * ck, "-")
                           }
      lw_num           = sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max),"-")
      w_num            = exp(lw_num)
      weight           = colMeans(w_num/rowSums(w_num))
      weight           = as.matrix(weight,ncol=1)
      rownames(weight) = c("GBM","L2B","RF","BAG","BART", paste("LM",1:(length(weight)-5), sep = ""))
      weight_se        = as.matrix(apply(w_num, 2, sd),ncol=1)/sqrt(no_rep)
      object           = list(weight = weight, weight_se = weight_se, candi_models = candi_models)
      class(object)    = "UARM"
      object
}
