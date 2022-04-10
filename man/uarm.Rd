\name{uarm}
\alias{uarm}
\title{Universal Adaptive Regression by Mixing with low-dimensional inputs}
\description{
Universal Adaptive Regression by Mixing (UARM) provides an adaptive model averaging with both linear models and nonparamatric methods considered as candidates. The nonparamatric methods include Generalized Boosted Regression modeling (GBM), L2Boosting (L2B), Random Forests (RF), Bagging (BAG), and Bayesian Additive Regression Trees (BART) on low-dimensional inputs.}
\usage{
uarm(x,y,factorID=NULL,candi_models,n_train=ceiling(n/2),no_rep=20,psi=0.1,
    method='L1-UARM',prior=TRUE,p0=0.5)
}
\arguments{
  \item{x}{Matrix of predictors.}

  \item{y}{Response variable.}

  \item{factorID}{Indication on whether there are categorical variables among the predictors. If factorID= NULL, the predictors are all continuous or have the identifiable categorical variables; If factorID=\code{'colnames'} or the location numbers of categorical variables, the name or location of variables provided by the user are treated as categorical variables in the linear model. The default factorID is NULL.}

  \item{candi_models}{Set to 1 for nested subset models in the order given in predictors; set to 2 for all combinations of subsets; input an m*p matrix, where m is the number of models to be combined, and each row of which is a 0/1 indicator vector representing whether each variable is included/excluded in the model. The default is 2.}

  \item{n_train}{Size of training set when the weight function is \code{L1-UARM} or \code{UARM}. The default value is \code{n_train=ceiling(n/2).}}

  \item{no_rep}{Number of replications when the weight function is \code{L1-UARM} and \code{UARM}. The default value is \code{no_rep=50}.}

  \item{psi}{A positive number to control the influence of the prior weight on the models. The default value is 0.1.}

  \item{prior}{Whether to use prior in the weighting function. The default is \code{TRUE}.}

  \item{method}{The method for calculating weights. Users can choose between \code{'L1-UARM'} and \code{'UARM'}. The default is \code{'L1-UARM'}.}

  \item{p0}{Prior probabilities of parametric methods. The default value is 0.5.}

}

\details{
See the paper provided in Reference section.
}

\value{
A \code{'uarm'} object is retured. The components are:

\item{weight}{The weight for each candidate model.}
\item{weight_se}{The standard error of the weights of the candidate models over the data-splittings.}

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



# compute weight and weight_se for candidate models using L1-UARM
#all subsets candidate models
Lu1<-uarm(x,y,factorID=NULL,candi_models=2,n_train=ceiling(n/2),no_rep=50,psi=0.1,
     method = 'L1-UARM', prior = TRUE, p0 = 0.5)
Lw1<-Lu1$weight
Ls1<-Lu1$weight_se

# compute weight and weight_se for candidate models using UARM
Lu2<-uarm(x,y,factorID=NULL,candi_models=2,n_train=ceiling(n/2),no_rep=50,psi=0.1,
     method = 'UARM', prior = TRUE, p0 = 0.5)
Luw2<-Lu2$weight
Ls2<-Lu2$weight_se

# output the candidate models
candi_models<-Lu2$candi_models

# simulation with categorical variables
n=100
x1<-rnorm(n)
x2<-rnorm(n)
x3<-rnorm(n)
x4<-factor(sample(1:5,n,replace=T),levels=c(1:5))
X<-data.frame(x1,x2,x3,x4)
Z<-as.matrix(model.matrix(~.-1,data=as.data.frame(X)))[,-4]
mu<-Z\%*\%c(0.1,0.3,0.5,1,-2,4,-3)
y<-mu+rnorm(n,0,3)

# compute weight for candidate models using MMA with nested subsets candidate models
Lu3<-uarm(X,y,factorID='x4',candi_models=2,n_train=ceiling(n/2),no_rep=50,psi=0.1,
     method='UARM',prior=TRUE,p0=0.5)
Luw3<-Lu3$weight
Ls3<-Lu3$weight_se

# early COVID-19 data in China
data(covid19)
y<-covid19[,1]
x<-covid19[,-1]
n<-length(y)

# compute weight and weight_se for model using L1-UARM
# all subsets candidate models
LC<-uarm(x,y,factorID=NULL,candi_models=2,n_train=ceiling(n/2),no_rep=50,psi=0.1,
    method='L1-UARM',prior=TRUE,p0=0.5)
LCw<-LC$weight
LCs<-LC$weight_se

# compute weight and weight_se for candidate models using UARM
LC2<-uarm(x,y,factorID=NULL,candi_models=2,n_train=ceiling(n/2),no_rep=50,psi=0.1,
     method='UARM',prior=TRUE,p0=0.5)
LCw2<-LC2$weight
LCs2<-LC2$weight_se
}


