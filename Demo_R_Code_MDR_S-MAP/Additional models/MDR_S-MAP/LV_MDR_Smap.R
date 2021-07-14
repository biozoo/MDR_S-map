# R-code of  Multiview Distance Regularized (MDR) S-map
# Article: Reconstructing high dimensional time-varying interaction networks for natural dynamical systems 
# Authors: Chun-Wei Chang, Takeshi Miki, Masayuki Ushio, Po-Ju Ke, Hsiao-Pei Lu, Fuh-Kwo Shiah, and Chih-hao Hsieh
# Published in: 

# R code for
Root_Path <- 'D:\\data\\Meta_analysis\\AdpSmap\\Code'
Sub_Path <- '\\Demo_R_Code_MDR_S-MAP\\Additional models\\'
OPath <- paste(Root_Path,Sub_Path,sep='')
setwd(OPath)

library(vegan)
library(rEDM)
library(doParallel)
library(parallel)
library(foreach)
library(Kendall)
library(MASS)
library(dplyr)
library(glmnet)
source(paste(Root_Path,'\\Demo_R_Code_MDR_S-MAP\\Demo_MDR_function.R',sep=''))
seed <- 49563
set.seed(seed)

#############################################################
#### Data Preparation 
# Read dataset
da.range <- 901:1000 # Subsample for data analysis
out.sample <- T # T/F for out-of-sample forecast
nout <- 2  # number of out-of-sample

(da.name <- 'LV')
do=read.csv(paste('./Data/','LV_dynamics20210611.csv',sep=''),header=T,stringsAsFactors = F)
do <- do[do[,1]%in%da.range,]
dot <- do[,1]
do <- do[,-1]
colnames(do) <- paste('X',1:ncol(do),sep='_')
ndo <- nrow(do)
nin <- ndo-nout # library sample size

# Excluding rare species with no clear temporal dynamics 
dto <- apply(do,1,sum) # total abundance
dpo <- do*repmat(dto^-1,1,ncol(do))
# Using a slightly more strict criteria (mean abundance >1.2*10^-3) to select species numbers comparable with other models  
pcri <- 0;bcri <- 1.2*10^-3; # criteria for selecting species based on proportion of presnece (pcri) and mean relative abundance (bri) 
doind2 <- (apply(dpo,2,mean,na.rm=T)>(bcri))&((apply(do>0,2,sum,na.rm=T)/nrow(do))>pcri) # index for selected species 
exsp2 <- setdiff(1:ncol(do),which(doind2))   # index for rare species 
do <- do[,-exsp2]                            # Dataset excluded rare species
(nsp <- ncol(do))                            # number of selected species

# In-sample
do.mean <- apply(do[1:nin,],2,mean,na.rm=T)  # mean abundance in in-sample
do.sd <- apply(do[1:nin,],2,sd,na.rm=T)      # SD of abundance in in-sample
d <- do[1:(nin-1),]                          # In-sample dataset at time t
d_tp1 <- do[2:(nin),]                        # In-sample dataset at time t+1
ds <- (d-repmat(do.mean,nrow(d),1))*repmat(do.sd,nrow(d),1)^-1 # Normalized in-sample dataset at time t
ds_tp1 <- (d_tp1-repmat(do.mean,nrow(d_tp1),1))*repmat(do.sd,nrow(d_tp1),1)^-1 # Normalized in-sample dataset at time t+1

# Out-sample
if(out.sample&nout!=0){
  d.test <- do[nin:(ndo-1),]                 # Out-of-sample dataset at time t 
  dt_tp1 <- do[(nin+1):ndo,]                 # Out-of-sample dataset at time t+1
  ds.test <- (d.test-repmat(do.mean,nrow(d.test),1))*repmat(do.sd,nrow(d.test),1)^-1 # Normalized out-of-sample dataset at time t
  dst_tp1 <- (dt_tp1-repmat(do.mean,nrow(dt_tp1),1))*repmat(do.sd,nrow(dt_tp1),1)^-1 # Normalized out-of-sample dataset at time t+1
}else{d.test <- dt_tp1 <- ds.test <- NULL}

# Compiled data at time t 
ds.all <- rbind(ds,ds.test)

#############################################################
# Find the optimal embedding dimension & nonlinearity parameter for each variable 
# based on univariate simplex projection and S-map, respectively

# Univariate simplex projection
Emax <- 10
cri <- 'rmse' # model selection 
Ed <- NULL
forecast_skill_simplex <- NULL
for(i in 1:ncol(ds)){
  spx.i <- simplex(ds[,i],E=2:Emax)
  Ed <- c(Ed,spx.i[which.min(spx.i[,cri])[1],'E'])
  forecast_skill_simplex <- c(forecast_skill_simplex,spx.i[which.min(spx.i[,cri])[1],'rho'])
}
Ed # The optimal embedding dimension for each variable
forecast_skill_simplex # Forecast skills for each variable based on simplex projection

######################################################################
# Find causal variables by CCM analysis for multiview embedding
# Warning: It is time consuming for calculating the causation for each node
# CCM causality test for all node pairs 
do.CCM <- F 
if(do.CCM){ 
  ccm.out <- ccm.fast.demo(ds,Epair=T,cri=cri,Emax=Emax)
  ccm.sig <- ccm.out[['ccm.sig']]
  ccm.rho <- ccm.out[['ccm.rho']]
  if(SaveFile){
    # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
    write.csv(ccm.sig,paste('./Output/','ccm_sig_',da.name,'_nin',nin,'_demo_NEW.csv',sep=''),row.names=F)
    write.csv(ccm.rho,paste('./Output/','ccm_rho_',da.name,'_nin',nin,'_demo_NEW.csv',sep=''),row.names=F)
  }
}

ccm.sig <- read.csv(paste('./Output/','ccm_sig_',da.name,'_nin',nin,'_demo.csv',sep=''),header=T,stringsAsFactors = F)
ccm.rho <- read.csv(paste('./Output/','ccm_rho_',da.name,'_nin',nin,'_demo.csv',sep=''),header=T,stringsAsFactors = F)

######################################################################
# Perform multiview embedding analysis for each node
# Warning: It is time consuming for running multview embedding for each nodes
do.multiview <- F
if(do.multiview){
  esele_lag <- esim.lag.demo(ds,ccm.rho,ccm.sig,Ed,kmax=10000,kn=100,max_lag=3,Emax=Emax)
  # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
  if(SaveFile){write.csv('./Output/',esele_lag,paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep=''),row.names=F)}
}

esele <- read.csv(paste('./Output/','eseleLag_',da.name,'_nin',nin,'_demo.csv',sep=''),header=T)


####################################################
## The computation of multiview distance
dmatrix.mv <- mvdist.demo(ds,ds.all,esele)
dmatrix.train.mvx <- dmatrix.mv[['dmatrix.train.mvx']]
dmatrix.test.mvx <- dmatrix.mv[['dmatrix.test.mvx']]

######## Leave-one-out cross-validation for finding the optimal parameters for MDR S-map analysis
### Warning: The cross-validation is the most time-consuming step in MDR S-map requiring massive computations and .  
### Thus, we recommend dividing job into smaller parts (sub.da>1)  or used parallel computation (parall=T, ncore>=1)
### The example showing below divided the parameter space into five parts and ran independently (sub.da=5).
do.MDR.CV <- F
### The parameter cv.unit determines the precision of selected parameters and strongly influences computation time.
### In our cases, we used cv.unit=0.025 to obtain more precise estimations
### This parameter may be adjusted to 0.05 or even 0.1, depending on how sensitive the results to parameter precision. 
cv.unit <- 0.025
alpha.so <- seq(0, 1, cv.unit);            # Sequence of alpha
sub.da <- 5                                # Divide the computation job into five parts 
afsp <- eqsplit(1:length(alpha.so),sub.da) # Divide the parameter space based on alpha parameter
alf <- 1                                   # Run CV in the first parameter subset 

# Cross-validation of MDR analysis    
if(do.MDR.CV){
  alpha.s <- alpha.so[afsp[alf,1]:afsp[alf,2]] # Subset parameter pace
  cv.ind <- cv.MDR.demo(ds, ds_tp1, dmatrix.list=dmatrix.train.mvx, 
                        parall=F, ncore=1, keep_intra=T,alpha.seq=alpha.s)
  # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
  if(SaveFile){write.csv(cv.ind,paste('./Output/',da.name,'_nin',nin,'_cvunit',cv.unit,'_alph',alpha.s[1]*100,'_cvout_Nmvx_Rallx.csv',sep=''),row.names=F)}
}
# Repeat the CV under different parameter subsets by changing alf=2,3..,sub.da                                  

################################################################################
# Compiled the CV results tested under different parts of parameter space
CompileCV=F
if(CompileCV){
  cv.ind <- NULL; for(alf in 1:nrow(afsp)){
    cv.ind <- rbind(cv.ind,read.csv(paste('./Output/',da.name,'_nin',nin,'_cvunit',cv.unit,'_alph',alpha.so[afsp[alf,1]]*100,
                                          '_cvout_Nmvx_Rallx_demo.csv',sep=''),header=T))
  }
  # Select the optimal parameter set with the minimal MSE
  paracv.demo <- secv.demo(cv.ind)
  write.csv(paracv.demo,paste('./Output/',da.name,'_nin',nin,'_cvunit',cv.unit,
                              '_OptimalCV_Nmvx_Rallx_NEW.csv',sep=''),row.names = F)
}


############################################################
# Fitting MDR S-map based on the parameters selected by CV
do.MDR <- F
cv.unit <- 0.025                           
ptype <- 'aenet'                           # enet:elastic-net or msaenet: adaptive elastic-net

# Select the optimal parameter set with the minimal MSE
paracv.demo <- read.csv(paste('./Output/',da.name,'_nin',nin,'_cvunit',cv.unit,'_OptimalCV_Nmvx_Rallx.csv',sep=''))

if(do.MDR){
  # Fitting the MDR S-map
  smap.demo <- MDRsmap.demo(paracv=paracv.demo,ptype=ptype,keep_intra=T,out.sample=T,
                            ds,ds_tp1,ds.test,dst_tp1,
                            dmatrix.list=dmatrix.train.mvx,
                            dmatrix.test.list=dmatrix.test.mvx)
  
  # Save forecast skills
  nr.out <- smap.demo[['nr.out']];
  # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
  if(SaveFile){
    write.csv(nr.out,paste(da.name,'_nin',nin,'_cvunit',cv.unit,'_',ptype,'_nrout_Nmvx_Rallx_demo_NEW.csv',sep=''),row.names=F)
    # Save interaction Jacobian matrices at all time points
    write.csv(smap.demo[['jcof']],paste(da.name,'_nin',nin,'_cvunit',cv.unit,'_',ptype,'_jcof_Nmvx_Rallx_demo_NEW.csv',sep=''),row.names=F)
  }
}
