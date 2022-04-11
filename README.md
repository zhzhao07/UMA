# UMA
This package provides adaptive model averaging (MA) with both linear and nonparamatric methods. It also allows the use of other averaging methods such as smoothed information criteria and Mallow's MA.   
Authors: Li Wen <wlwendy1008@163.com>, Zhihao Zhao <zhzhao@cueb.edu.cn>, Yuhong Yang <yyang@stat.umn.edu>.  
### Installation
To install this package in R, run the following commands:  
`library(devtools)`  
`devtools::install_github("zhzhao07/UMA")`  

### Example usage:
Below is an example of using the function gma, uarm, and uma.predict:  
```#generate simulation data  
`library(UMA)`  
`n <- 50`  
`p <- 8`  
`beta <- c(3,1.5,0,0,2,0,0,0)`  
`b0 <- 1`  
`x <- matrix(rnorm(n*p,0,1),nrow=n,ncol=p)`  
`e <- rnorm(n,0,3)`  
`y <- x%*%beta+b0+e`  


#the weighted estimation using L1-ARM, MMA and SFIC with all subsets candidate models  
`Cl  <- gma(x,y,factorID=NULL,method='L1-ARM',candi_models=2)$wbetahat`  
`Cm  <- gma(x,y,factorID=NULL,method='MMA',candi_models=2)$wbetahat`  
`Csf <- gma(x,y,factorID=NULL,method='SFIC',candi_models=2)$wbetahat`  

#compute weight and weight_se for model using L1-UARM with all subsets candidate models  
`LC  <- uarm(x,y,factorID=NULL,candi_models=2,n_train=ceiling(n/2),no_rep=50,psi=0.1,`  
`method='L1-UARM',prior=TRUE,p0=0.5)`    
`LCw <- LC$weight`    
`LCs <- LC$weight_se`  

#compute weight and weight_se for candidate models using UARM with all subsets candidate models  
`LC2  <- uarm(x,y,factorID=NULL,candi_models=2,n_train=ceiling(n/2),no_rep=50,psi=0.1,`  
`        method='UARM',prior=TRUE,p0=0.5)`    
`LCw2 <- LC2$weight`    
`LCs2 <- LC2$weight_se`  

#compute weight and weight_se for model using MMA and SFIC with all subsets candidate models  
`mmaw <- gma(x,y,factorID=NULL,method='MMA',candi_models=2)$weight`    
`sficw <- gma(x,y,factorID=NULL,method='SFIC',candi_models=2)$weight`  

#compute the prediction by methods L1-UARM, MMA, SFIC and BMA  
`l1uarm.predict <- uma.predict(x,y,factorID=NULL,newdata=x,candi_models=2,`  
`            weight=LCw,method='L1-UARM',dim='L')$pre_out`    
`mma.predict <- uma.predict(x,y,factorID=NULL,newdata=x,candi_models=2,`  
`            weight=mmaw,method='MMA',dim='L')$pre_out`  
`sfic.predict <- uma.predict(x,y,factorID=NULL,newdata=x,candi_models=2,`  
`            weight=sficw,method='SFIC',dim='L')$pre_out`  

#the BMA prediction does not depend on candidate models  
`bma.predict  <- uma.predict(x,y,factorID=NULL,newdata=x,candi_models=2,`  
`            method='BMA',dim='L')$pre_out`
