\name{uarm_h}
\alias{uarm_h}
\title{Universal Adaptive Regression by Mixing with high-dimensional inputs}
\description{
High-dimensional Universal Adaptive Regression by Mixing (uarm_h) provides an adaptive model averaging with both linear and nonparamatric methods considered as candidates (Generalized Boosted Regression modeling (GBM), L2Boosting (L2B) and Random Forests (RF)) on high-dimensional inputs.}

\usage{
uarm_h(x, y, factorId=NULL,candidate='H4', method='L1-UARM', n_train=ceiling(n/2),
      p0=0.5, no_rep=50, psi=1, prior = TRUE)
}
\arguments{
  \item{x}{Matrix of predictors.}

  \item{y}{Response variable.}

  \item{factorID}{Indication on whether there are categorical variables among the predictors.} { If factorID= NULL, the predictors are all continuous or have the identifiable categorical variables; If factorID=`colnames' or the location numbers of categorical variables, the name or location of variables provided by the user are treated as categorical variables in the linear model. The default factorID is NULL.}

  \item{candidate}{Method for preparing candidate models.}{ The method of candidate selection differs depending on whether it contains categorical variables. If the predictors are all continuous variables, the candidate= `H4', the candidate models are on solution paths of 4 common methods, which are lasso, adaptive lasso, SCAD and MCP;  the candidate= `H2', the candidate models are on solution paths of 2 common methods, which are lasso, adaptive lasso; the candidate=`H1', the candidate models are on solution paths of the lasso. Otherwise, the candidate= `CH3', the candidate models are on solution paths of 3 group selection methods in categorical regression, which are group lasso, group MCP and group SCAD, treating the categorical variable as the indivual groups. If candidate= `H0', the candidate model should be input by the user. When the method is MCV, candidate are not required.}

  \item{method}{The method for calculating weights. The method= `UARM' is the Universal Adaptive Regression by Mixing method ; the method= `L1-UARM' is the L1 Universal Adaptive Regression by Mixing method; the method= `UARM.rf' is the Universal Adaptive Regression by Mixing method using the random forest to estimate the standard deviation of random error for candidate models; the method= `L1-UARM.rf' is the L1 Universal Adaptive Regression by Mixing method using the random forest to estimate the standard deviation of random error for candidate models.}

  \item{n_train}{Size of training set when the weight function is \code{UARM}  or  \code{L1-UARM}. The default value is \code{n_train=ceiling(n/2).}}

  \item{no_rep}{Number of replications when the weight function is \code{UARM} and \code{L1-UARM}. The default value is \code{no_rep=20}.}

  \item{psi}{A positive number to control the improvement of the prior weight. The default value is \code{psi=1.}}

  \item{p0}{Prior probability of parametric methods. The default value is \code{p0=0.5.}}

  \item{candi_models}{This component is used by the user to input the candidate model matrix.}  {It is a matrix of candidate models, each row of which is a 0/1 indicator vector representing whether each variable is included/excluded in the model. For details see example section. }

  \item{prior}{Whether to use prior in the weighting function. The default is \code{TRUE}.}
}

\details{
  See the paper provided in reference section.
}

\value{
  A "uarm_h" object is retured. The components are:
  \item{weight}{The weight for each candidate model.}
  \item{weight_se}{The standard error of candidate models over the each data-splitting.}
  \item{candi_models}{The candidate models.}
  }



\examples{
#####linear case
library(mvtnorm) 
n=100
p=200
b=rep(0,len=p)
b[1:10]=2*c(0.878,0.636,-0.257,0.320,0.181,0.174,0.795,-0.523,0.228,-0.727)
X=matrix(rnorm(n*p,0,1),nrow=n,ncol=p)
e=rnorm(n,0,0.5^2)
y=X\%*\%b + e
#compute the weight and prediction
uarm=uarm_h(x=X,y,factorID=NULL,candidate='H4',n_train=n/2,method='UARM',
            p0=0.5,no_rep=50,psi=1, prior = TRUE)
wy=uarm$weight
uarm$weight_se
l1uarm=uarm_h(x=X,y,factorID=NULL,candidate='H4',n_train=n/2,method='L1-UARM',
              p0=0.5,no_rep=50,psi=1, prior = TRUE)
wy1=l1uarm$weight
l1uarm$weight_se

##compare the nonparametric weights and parametric weights,
##the nonparametric weights almost zero. (the first three weights are nonparametric).
weight_U=cbind(sum(uarm$weight[1:3]),sum(uarm$weight[-c(1:3)]))
weight_L1U=cbind(sum(l1uarm$weight[1:3]),sum(l1uarm$weight[-c(1:3)]))


# Nonlinear Case
library(mvtnorm) 
n=100
p=400
sigma0=1
###beta and beta0;decay
b=rep(0,p)
alpha=1
for (j in 1:10){
  b[j]=1/j
}
b=b/sum(b)
b0=5.2 ##nonpara degree
###################cov setting
Sig= matrix(0,p,p)
rho = 0.5
for(i in 1:p)
{
  for(j in 1:p)
  {
    Sig[i,j]=rho^abs(i-j)
  }
}
#################train data
X=matrix(rmvnorm(n,matrix(0,ncol=1,nrow=p),Sig),nrow=n)
mu0=X\%*\%b+b0*X[,1]*X[,2]
y=mu0+rnorm(n,0,sigma0)##normal distribution
#compute the weight and prediction
l1uarm_rf=uarm_h(x=X,y,factorID=NULL,candidate='H4',n_train=n/2,
                 method='L1-UARM.rf',p0=0.5,no_rep=50,psi=1, prior = TRUE)
weight=l1uarm_rf$weight
#####comparr the nonparametric and parametric weights, 
#####the parametric weights are small.
weight_U=cbind(sum(l1uarm_rf$weight[1:3]),sum(l1uarm_rf$weight[-c(1:3)]))
weight_se=l1uarm_rf$weight_se


#######categorical regression
library(nnet)
library(CatReg)
library(mvtnorm) 
n=100
p=200#
sigma0=1
#######1*p,beta
b <- rep(0,len=p)
##decay
for (j in 1:7){
  b[j]=1/j
}
###################cov setting
Sig= matrix(0,p,p)
rho = 0.5
for(i in 1:p)
{
  for(j in 1:p)
  {
    Sig[i,j]=rho^abs(i-j)
  }
}
x0=UniformDesignMatrix(n, 1, 3)
X_c=class.ind(as.matrix(x0))
X0=matrix(rmvnorm(n,matrix(0,ncol=1,nrow=p),Sig),nrow=n)
mu0=X_c[,-1]\%*\%c(2,4)+X0\%*\%b
y=mu0+rnorm(n,0,sigma0)##normal distribution
X=cbind(x0,X0)###combine dummy covariable
###categorical regression
isarm=uarm_h(x=X,y,factorID=1,candidate='CH3',method='UARM',n_train=n/2,
            no_rep=50,psi=1,prior=TRUE)
iluarm=uarm_h(x=X,y,factorID=1,candidate='CH3',method='L1-UARM',n_train=n/2,
             no_rep=50,psi=1,prior=TRUE)
weight=cbind(isarm$weight,iluarm$weight)


}
