\name{gma}
\alias{gma}
\title{general model averaging for low-dimensional inputs}
\description{gma provides model averaging for linear regrssion with low dimensional inputs (no more than 20 covariates). The MA methods included are SAIC, SBIC, SFIC, ARM, L1-ARM, MMA and JMA.}
\usage{
  gma(x,y,factorID=NULL,method='L1-ARM',candi_models=2,n_train=ceiling(n/2),
      no_rep=50)
}
\arguments{
  \item{x}{Matrix of predictors.}

  \item{y}{Response variable.}

  \item{factorID}{Indication on whether there are categorical variables among the predictors. If factorID= NULL, the predictors are all continuous or have the identifiable categorical variables; If factorID=\code{'colnames'} or the location numbers of categorical variables, the name or location of variables provided by the user are treated as categorical variables in the linear model. The default factorID is NULL.}

  \item{method}{The method for calculating weights. The method= \code{'SAIC'} is the Smooth-AIC method; the method= \code{'SBIC'} is the Smooth-BIC method; the method= \code{'SFIC'} is the Smooth-FIC method; the method= \code{'ARM'} is the Adaptive Regression by Mixing method; the method= \code{'L1-ARM'} is the L1 Adaptive Regression by Mixing method; the method= \code{'MMA'} is the Mallow's Model Averaging (MMA); the method= \code{'JMA'} is the Jackknife Model Averaging (JMA). The default is \code{'L1-ARM'}.}

  \item{candi_models}{Set to 1 for nested subset models in the order given in predictors; set to 2 for all combinations of subsets; input an m*p matrix, where m is the number of models to be combined, and each row of which is a 0/1 indicator vector representing whether each variable is included/excluded in the model. The default is 2.}

  \item{n_train}{Size of training set when the weight function is \code{L1-ARM} or \code{ARM} with \code{prior=TRUE}. The default value is \code{n_train=ceiling(n/2).}}

  \item{no_rep}{Number of replications when the weight function is \code{L1-ARM} and \code{ARM}. The default value is \code{no_rep=50}.}
}

\details{
  See the paper provided in Reference section.
}

\value{
  A \code{'gma'} object is retured. The components are:
  \item{weight}{The weight for each candidate model.}
  \item{weight_se}{The standard error of the weights of the candidate models over the data-splittings under the method= \code{'ARM'} or method=\code{'L1-ARM'}.}
  \item{wbetahat}{The weighted estimation of the coefficients.}
  \item{betahat}{The coefficients matrix estimated by candidate models.}
  \item{candi_models}{The candidate models.}
  }

\examples{
# generate simulation data
n<-50
p<-8
beta<-c(3,1.5,0,0,2,0,0,0)
b0<-1
x<-matrix(rnorm(n*p,0,1),nrow=n,ncol=p)
e<-rnorm(n,0,3)
y<-x\%*\%beta+b0+e

# compute weight for candidate models using L1-ARM, JMA and SAIC with nested subsets candidate models
lw<-gma(x,y,factorID=NULL,method='L1-ARM',candi_models=1)$weight
jw<-gma(x,y,factorID=NULL,method='JMA',candi_models=1)$weight
saw<-gma(x,y,factorID=NULL,method='SAIC',candi_models=1)$weight

# output the candidate models used for method L1-ARM
candi_models<-gma(x,y,factorID=NULL,method='L1-ARM',candi_models=1)$candi_models

# simulation with categorical variables
n<-100
x1<-rnorm(n)
x2<-rnorm(n)
x3<-rnorm(n)
x4<-factor(sample(1:5,n,replace=T),levels=c(1:5))
X<-data.frame(x1,x2,x3,x4)
Z<-as.matrix(model.matrix(~.-1,data=as.data.frame(X)))[,-4]
mu<-Z\%*\%c(0.1,0.3,0.5,1,-2,4,-3)
y<-mu+rnorm(n,0,3)

# compute weight for candidate models using MMA with nested subsets candidate models
mmaw <- gma(X, y, factorID = 'x4', method = 'MMA', candi_models = 1)$weight

# early COVID-19 data in China
data(covid19)
y<-covid19[,1]
x<-covid19[,-1]
n<-length(y)

# the weighted estimation using L1-ARM, MMA and SFIC with all subsets candidate models
Cl<-gma(x,y,factorID=NULL,method='L1-ARM',candi_models=2)$wbetahat
Cm<-gma(x,y,factorID=NULL,method='MMA',candi_models=2)$wbetahat
Csf<-gma(x,y,factorID=NULL,method='SFIC',candi_models=2)$wbetahat
}
