\name{uma.predict}
\alias{uma.predict}
\title{Prediction for Universal Adaptive Regression by Mixing}
\description{
The predictions based on different MA methods, including SAIC (SAICp), SBIC (SBICp), SFIC, ARM, L1-ARM, UARM, L1-UARM, MMA, JMA, PMA, BMA and MCV.}
\usage{
uma.predict(x,y,factorID=NULL,newdata,candi_models,weight,method,dim)
}
\arguments{
  \item{x}{Matrix of predictors.}

  \item{y}{Response variable.}

  \item{factorID}{Indication on whether there are categorical variables among the predictors.} { If factorID= NULL, the predictors are all continuous or have the identifiable categorical variables; If factorID=\code{'colnames'} or the location numbers of categorical variables, the name or location of variables provided by the user are treated as categorical variables in the linear model. The default factorID is NULL.}

  \item{candi_models}{The candidate models under specific method, you can be calculated by gma, gma_h, uarm, and uarm_h functions, as shown in the examples.}


  \item{newdata}{An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.}

  \item{weight}{The weights of candidate models under specific methods, you can be calculated by gma, gma_h, uarm, and uarm_h functions, as shown in the examples}

  \item{method}{The method= \code{'UARM'} is the Universal Adaptive Regression by Mixing method; the method= \code{'L1-UARM'} is the L1 Universal Adaptive Regression by Mixing method; the method= \code{'SAIC'} is the Smooth-AIC method; the method= \code{'SBIC'} is the Smooth-BIC method; the method= \code{'SAICp'} is the Smooth-AIC method with the penalty term; the method= \code{'SBICp'} is the Smooth-BIC method with the penalty term; the method= \code{'SFIC'} is the Smooth-FIC method; the method= \code{'ARM'} is the Adaptive Regression by Mixing method; the method= \code{'L1-ARM'} is the L1 Adaptive Regression by Mixing method; the method= \code{'MMA'} is the Mallows Model Averaging (MMA); the method= \code{'JMA'} is the Jackknife Model Averaging (JMA); the method= \code{'UARM.rf'} is the Universal Adaptive Regression by Mixing method using the random forest to estimate the standard deviation of random error for candidate models; the method= \code{'L1-UARM.rf'} is the L1 Universal Adaptive Regression by Mixing method using the random forest to estimate the standard deviation of random error for candidate models; the method= \code{'PMA'} is the Parsimonious Model Averaging; the method= \code{'MCV'} is the Cross-validation for Model Averaging (MCV).}


 \item{dim}{High-dimensional or low-dimensional methods are used for prediction. If dim =\code{'H'}, high-dimensional methods are used; otherwise, low-dimensional methods are used.}
}

\details{
See the paper provided in Reference section.
}

\value{
A \code{'uma.predict'} object is retured. The components is:

\item{pre_out}{The prediction by given method.}

}



\examples{
### low dimension

# generate simulation data
n<-50
p<-8
beta<-c(3,1.5,0,0,2,0,0,0)
b0<-1
x<-matrix(rnorm(n*p,0,1),nrow=n,ncol=p)
e<-rnorm(n,0,3)
y<-x\%*\%beta+b0+e

# user supplied candidate models
candi_models<-rbind(c(0,0,0,0,0,0,0,1),
                    c(0,1,0,0,0,0,0,1),
                    c(0,1,1,1,0,0,0,1),
                    c(0,1,1,0,0,0,0,1),
                    c(1,1,0,1,1,0,0,0),
                    c(1,1,0,0,1,0,0,0))

# compute weight for candidate models using L1-UARM
weightL<-uarm(x,y,factorID=NULL,candi_models=candi_models,n_train=ceiling(n/2),
         no_rep=50,psi=1,method='L1-UARM',prior=TRUE,p0=0.5)$weight

# compute the prediction by method L1-UARM
luma.predict<-uma.predict(x,y,factorID=NULL,newdata=x,candi_models=candi_models,
              weight=weightL,method='L1-UARM',dim='L')$pre_out

# early COVID-19 data in China
data(covid19)
y<-covid19[,1]
x<-covid19[,-1]
n<-length(y)

# compute the predicts for L1-UARM, MMA and SFIC
# user supplied all subsets candidate models
Cl1uarmw<-uarm(x,y,factorID=NULL,candi_models=2,n_train=ceiling(n/2),no_rep=50,
          psi=1,method='L1-UARM',prior=TRUE,p0=0.5)$weight
Cmmaw<-gma(x,y,factorID=NULL,method='MMA',candi_models=2)$weight
Csficw<-gma(x,y,factorID=NULL,method='SFIC',candi_models=2)$weight

# compute the prediction by methods L1-UARM, MMA, SFIC and BMA
cl1uarm.predict<-uma.predict(x,y,factorID=NULL,newdata=x,candi_models=2,
                             weight=Cl1uarmw,method='L1-UARM',dim='L')$pre_out
cmma.predict<-uma.predict(x,y,factorID=NULL,newdata=x,candi_models=2,
                          weight=Cmmaw,method='MMA',dim='L')$pre_out
csfic.predict<-uma.predict(x,y,factorID=NULL,newdata=x,candi_models=2,
                           weight=Cmmaw,method='SFIC',dim='L')$pre_out

#The BMA prediction does not depend on candidate models
cbma.predict<-uma.predict(x,y,factorID=NULL,newdata=x,candi_models=2,
                          method='BMA',dim='L')$pre_out



###high dimension
library(mvtnorm) 
n1=100;n2=1000
p=200
sigma0=1
######
b=rep(0,len=p) #1*p,beta
for(j in 1:12){
  b[j]=2/j
}
# cov setting
Sig = matrix(0,p,p)
rho = 0.5
for(i in 1:p)
{
  for(j in 1:p)
  {
    Sig[i,j] = rho^abs(i-j)
  }
}
# new data
X=matrix(rmvnorm(n1,matrix(0,ncol=1,nrow=p),Sig),nrow=n1)
X_test=matrix(rmvnorm(n2,matrix(0,ncol=1,nrow=p),Sig),nrow=n2)
mu0=X\%*\%b
mu_test=X_test\%*\%b
y=mu0+rnorm(n1,0,sigma0)##normal distribution
##########the prediction on each methods
g1=gma_h(x=X,y,factorID=NULL,candidate='H4',method='SBICp',psi=1, prior=TRUE)
pre_out1=uma.predict(x=X,y,factorID=NULL,newdata=X_test,method='SBICp',weight=g1$weight,
                    candi_models=g1$candi_models,dim='H')$pre_out

g2=gma_h(x=X,y,factorID=NULL,candidate='H4',method='PMA',lambda=log(n1))
pre_out2=uma.predict(x=X,y,factorID=NULL,newdata=X_test,method='PMA',weight=g2$weight,
                     candi_models=g2$candi_models,dim='H')$pre_out

g3=gma_h(x=X,y,factorID=NULL,method='MCV',alpha = 0.05)
pre_out3=uma.predict(x=X,y,factorID=NULL,newdata=X_test,method='MCV',weight=g3$weight,
                     candi_models=g3$candi_models,dim='H')$pre_out

g4=gma_h(x=X, y, factorID=NULL,candidate='H4',method='ARM',n_train=n1/2, no_rep=50, psi=1) 
pre_out4=uma.predict(x=X,y,factorID=NULL,newdata=X_test,method='ARM',weight=g4$weight,
                     candi_models=g4$candi_models,dim='H')$pre_out

g5=gma_h(x=X, y, factorID=NULL,candidate='H4',method='L1-ARM',n_train=n1/2, no_rep=50, psi=1)
pre_out5=uma.predict(x=X,y,factorID=NULL,newdata=X_test,method='L1-ARM',weight=g5$weight,
                     candi_models=g5$candi_models,dim='H')$pre_out

g6=uarm_h(x=X,y,factorID=NULL,candidate='H4',n_train=n1/2,method='UARM',
            no_rep=50,p0=0.5,psi=1, prior = TRUE)
pre_out6=uma.predict(x=X,y,factorID=NULL,newdata=X_test,method='UARM',weight=g6$weight,
                candi_models=g6$candi_models,dim='H')$pre_out

g7=uarm_h(x=X,y,factorID=NULL,candidate='H4',n_train=n1/2,method='L1-UARM',
          no_rep=50,p0=0.5,psi=1, prior = TRUE) 
pre_out7=uma.predict(x=X,y,factorID=NULL,newdata=X_test,method='L1-UARM',weight=g7$weight,
                candi_models=g7$candi_models,dim='H')$pre_out

####the performance of different methods
pre_out=cbind(pre_out1,pre_out2,pre_out3,pre_out4,pre_out5,pre_out6,pre_out7)
se=(pre_out-matrix(mu_test,n2,1)\%*\%rep(1,7))^2
colnames(se)=c('SBICp','PMA','MCV','ARM','L1-ARM','UARM','L1-UARM')
Pre=apply(se,2,mean)
Pre_se=apply(se,2,se)

}

