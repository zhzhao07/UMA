\name{gma_h}
\alias{gma_h}
\title{general model averaging for high-dimensional inputs}
\description{gma_h provides model averaging for linear regression with high-dimensional inputs. The Model Averaging methods included are SAICp, SBICp, ARM, L1-ARM, PMA, MCV.}
\usage{
  gma_h(x, y, factorID=NULL,candidate='H4',method='L1-ARM', psi=1,
        n_train=ceiling(n/2),no_rep=20, lambda=log(n), alpha=0.05, prior = TRUE)
}
\arguments{
  \item{x}{Matrix of predictors.}

  \item{y}{Response variable.}

  \item{factorID}{Indication on whether there are categorical variables among the predictors.} { If factorID= NULL, the predictors are all continuous or have the identifiable categorical variables; If factorID=`colnames' or the location numbers of categorical variables, the name or location of variables provided by the user are treated as categorical variables in the linear model. The default factorID is NULL.}

   \item{candidate}{Method for preparing candidate models.}{ The method of candidate selection differs depending on whether it contains categorical variables. If the predictors are all continuous variables, the candidate= `H4', the candidate models are on solution paths of 4 common methods, which are lasso, adaptive lasso, SCAD and MCP;  the candidate= `H2', the candidate models are on solution paths of 2 common methods, which are lasso, adaptive lasso; the candidate=`H1', the candidate models are on solution paths of the lasso. Otherwise, the candidate= `CH3', the candidate models are on solution paths of 3 group selection methods in categorical regression, which are group lasso, group MCP and group SCAD, treating the categorical variable as the indivual groups. If candidate= `H0', the candidate model should be input by the user. When the method is MCV, candidate are not required.}

  \item{candi_models}{This component is used by the user to input the candidate model matrix. It is a matrix of candidate models, each row of which is a 0/1 indicator vector representing whether each variable is included/excluded in the model. For details see example section.}


  \item{method}{The method for calculating weights.} {The method= `SAICp' is the Smooth-AIC method with the penalty term; the method= `SBICp' is the Smooth-BIC method with the penalty term; the method= `ARM' is the Adaptive Regression by Mixing method; the method= `L1-ARM' is the L1 Adaptive Regression by Mixing method; the method= `PMA' is the Parsimonious Model Averaging (PMA); the method= `MCV' is the Cross-validation for Model Averaging (MCV).}

  \item{n_train}{Size of training set when the weight function is \code{ARM}  or  \code{L1-ARM}. The default value is \code{n_train=ceiling(n/2).}}

  \item{no_rep}{Number of replications when the weight function is \code{ARM} or \code{L1-ARM}. The default value is \code{no_rep=20}.}

  \item{lambda}{It is the tunning parameter in \code{PMA}. The default is \code{log(n)}.}

  \item{alpha}{Threshold value for marginal correlation test in MCV. The default is \code{0.05}.}

  \item{psi}{A positive number to control the influence of the prior weight on the models. The default value is 1.}

  \item{prior}{Whether to use prior in the weighting function when the method=
               (`ARM',`L1-ARM',`SAIC',`SBIC'). The default is \code{TRUE}.}
}

\details{
  See the paper provided in Reference section.
}

\value{
  A \code{'gma_h'} object is retured. The components are:
  \item{weight}{The weight for each candidate model.}
\item{weight_se}{The standard error of the weights of the candidate models over the data-splittings under the method= `ARM' or method=`L1-ARM'.}
  \item{wbetahat}{The weighted estimation of the coefficients.}
  \item{betahat}{The coefficients matrix estimated by candidate models.}
  \item{candi_models}{The candidate models. NOTE: The \code{weights} and \code{candi_models} of MCV method, see the article for details. its outputs are mainly for prediction.}
  }



\examples{
library(mvtnorm) 
n=100
p=200
alpha=1
b <- rep(0,len=p)
for (j in 1:10){
  b[j]=(j^(-alpha-0.5))*sqrt(2*alpha)
}
b=5*b/sum(b)
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
mu0=X\%*\%b
y=mu0+rt(n,df=3)##t distribution
#the calculating for different method.
###################################################
g1=gma_h(x=X, y,factorID=NULL, candidate='H4',method='SAICp',psi=1) 
g2=gma_h(x=X, y, factorID=NULL,candidate='H4',method='ARM',n_train=n/2, no_rep=50, psi=1) 
g3=gma_h(x=X, y, factorID=NULL,candidate='H4',method='L1-ARM',n_train=n/2, no_rep=50, psi=1)
g4=gma_h(x=X, y, factorID=NULL,candidate='H4',method='PMA',lambda=log(n))
#####weight
weight=cbind(g1$weight,g2$weight,g3$weight,g4$weight)
weight_se=cbind(g2$weight_se,g3$weight_se)
#####coefficients estimation and prediction
wbetahat=cbind(g1$wbetahat,g2$wbetahat,g3$wbetahat,g4$wbetahat)
Xs=cbind(1,X)
pre=cbind(Xs\%*\%g1$wbetahat,Xs\%*\%g2$wbetahat,Xs\%*\%g3$wbetahat,Xs\%*\%g4$wbetahat)
se=(pre-matrix(mu0,n,1)\%*\%rep(1,4))^2
colnames(se)=c('SAICp','ARM','L1-ARM','PMA')
boxplot(se)
apply(se,2,mean)

##############categorical case
library(nnet)
library(CatReg)
library(mvtnorm) 
n=100
p=200
sigma0=1
######
b <- rep(0,len=p) #1*p,beta
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
x0= UniformDesignMatrix(n, 1, 3) ##categorical variable
X_c=class.ind(as.matrix(x0))
X0=matrix(rmvnorm(n,matrix(0,ncol=1,nrow=p),Sig),nrow=n)
mu0=X_c[,-1]\%*\%c(2,4)+X0\%*\%b
y=mu0+rnorm(n,0,sigma0) ##normal distribution
X=cbind(x0,X0)### combine dummy covariable

###categorical regression
isarm=gma_h(x=X,y,factorID=1,candidate='CH3',method='ARM',n_train=n/2,
            no_rep=50,psi=1,prior=TRUE)
pma=gma_h(x=X,y,factorID=1,candidate='CH3',method='PMA',lambda=log(n))

weight=cbind(isarm$weight,pma$weight)
wbeta=cbind(isarm$wbetahat,pma$wbetahat)

}
