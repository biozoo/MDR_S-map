# R code for comparing the reconstructed interaction Jacobian networks with their theoretical expectations 
# Article: Reconstructing large interaction networks from empirical time series data 
# Authors: Chun-Wei Chang, Takeshi Miki, Masayuki Ushio, Po-Ju Ke, Hsiao-Pei Lu, Fuh-Kwo Shiah, and Chih-hao Hsieh
# Published in: 

##############################################################################################
### Reconstructing interaction Jacobians for multi-species Ricker model 
### with time-varying interaction coefficients
rm(list = ls())
library(dplyr)
library(quantreg)

seed <- 49563
set.seed(seed)
# Set path
Root_Path <- 'D:\\data\\Meta_analysis\\AdpSmap\\Code'
Sub_Path <- '\\Demo_R_Code_MDR_S-MAP\\Additional models\\'
OPath <- paste(Root_Path,Sub_Path,sep='')
setwd(OPath)
source(paste(Root_Path,'\\Demo_R_Code_MDR_S-MAP\\Demo_MDR_function.R',sep=''))

# Dataset
(da.name <- 'RK_TV')
do <- read.csv(paste('./Data/',da.name,'_dynamics20210622.csv',sep=''),header=T,stringsAsFactors = F)
da.range <- 901:1000 # Subsample for data analysis
out.sample <- T # T/F for out-of-sample forecast
nout <- 2  # number of out-of-sample
do <- do[do[,1]%in%da.range,]
dot <- do[,1]
do <- do[,-1]
colnames(do) <- paste('X',1:ncol(do),sep='_')
ndo <- nrow(do)
nin <- ndo-nout # library sample size

# Excluding rare species with no clear temporal dynamics 
dto <- apply(do,1,sum) # total abundance
dpo <- do*repmat(dto^-1,1,ncol(do))
pcri <- 0;bcri <- 10^-3; # criteria for selecting species based on proportion of presence (pcri) and mean relative abundance (bri) 
doind2 <- (apply(dpo,2,mean,na.rm=T)>(bcri))&((apply(do>0,2,sum,na.rm=T)/nrow(do))>pcri) # index for selected species 
exsp2 <- setdiff(1:ncol(do),which(doind2))   # index for rare species 
do <- do[,-exsp2]                            # Dataset excluded rare species
(nsp <- ncol(do))                            # number of selected species

# In-sample
do.mean <- apply(do[1:nin,],2,mean,na.rm=T)  # mean abundance in in-sample
do.sd <- apply(do[1:nin,],2,sd,na.rm=T)      # SD of abundance in in-sample
dosdM <- repmat(c(do.sd)^-1,1,nsp)*repmat(c(do.sd),nsp,1)

# Reconstructed interaction Jacobian networks by MDR S-map
ints=read.csv(paste('./Output/',da.name,'_nin98_cvunit0.025_aenet_jcof_Nmvx_Rallx_demo.csv',sep=''),stringsAsFactors = F)

# Qualitative network
intsp <- ints[,-c(1:3)]
intsp[ints[,-c(1:3)]!=0] <- 1
intsp <- cbind(ints[,c(1:3)],intsp)

######################
# Theoretical networks exponential Jacobian matrix
jacobian_true=read.csv(paste('./Data/theoretical_DF_',da.name,'_demo.csv',sep=''),header = T)

######################
# Time-varying networks
########################
tsea <- unique(intsp[,'time'])
intsA <- JA <- JA2 <- NULL
cof <- NULL
TJ.ls <- EJ.ls <- list()
TJ.ols <- EJ.ols <- list()
cor.sp <- cor.sp.o <- NULL
for(i in tsea){
  # Theoretical interaction matrix
  JA.t <- as.matrix(dplyr::filter(jacobian_true,time==i)[,-1])
  # Scaled by the standard deviation of variables 
  JA.t2 <- JA.t*dosdM#*repmat(as.numeric(d[i,]),1,nsp)
  JA <- cbind(JA,c(JA.t))# time-varying true Jacobian 
  JA2 <- cbind(JA2,c(JA.t2))# time-varying true Jacobian (scaled)
  
  intsp.t <- intsp[intsp[,'time']==i&intsp[,'Insample']==1,-c(1:4)]
  ints.t <- as.matrix(ints[ints[,'time']==i&ints[,'Insample']==1,-c(1:4)])
  ints.t[intsp.t<1] <- 0
  rownames(intsp.t) <- rownames(intsp.t) <- colnames(intsp.t) <- colnames(ints.t)
  intsA <- cbind(intsA,c(ints.t))
  
  # Collect out-of sample networks
  to <- i+length(tsea)
  if(to<=length(tsea)+nout){
    JA.to <- as.matrix(dplyr::filter(jacobian_true,time==to)[,-1])*dosdM
    intsp.to <- intsp[intsp[,'time']==i&intsp[,'Insample']==0,-c(1:4)] # select network reconstructed by out-of-sample data
    ints.to <- as.matrix(ints[ints[,'time']==i&ints[,'Insample']==0,-c(1:4)])
    ints.to[intsp.t<1] <- 0
    TJ.ols[[i]] <- JA.to; EJ.ols[[i]] <- ints.to # Collecting out-of-sample networks
  }
  
  y.t <- ints.t; # Reconstructed network
  x.t <- JA.t2;  # Theoretical network
  TJ.ls[[i]] <- x.t; EJ.ls[[i]] <- y.t # Collecting in-sample networks
  
  # Node inference skill (min=0, max=1)
  cor.sp.t <- cor.sp.t.o <- NULL;
  for(j in 1:nsp){
    cor.sp.t <- c(cor.sp.t,cor(x.t[j,],y.t[j,]))
    cor.sp.t.o <- c(cor.sp.t.o,cor(x.t[,j],y.t[,j]))
  }
  cor.sp.t[cor.sp.t<0]=0
  cor.sp.t.o[cor.sp.t.o<0]=0
  cor.sp <- rbind(cor.sp,cor.sp.t);
  cor.sp.o <- rbind(cor.sp.o,cor.sp.t.o)
  
  # Statistical relationship between estimated interaction strength and theoretical expectation
  lm.t <- lm(unlist(c(y.t))~unlist(c(x.t)))
  lm.q <- rq(unlist(c(y.t))~unlist(c(x.t)))
  res.t <- unlist(c(x.t))-unlist(c(y.t))
  corr1 <- cor(x=unlist(c(x.t)), y=unlist(c(y.t)), method = 'pearson',use="pairwise.complete.obs")
  corr1[corr1<0]=0;
  
  cof <- rbind(cof,c(beta=as.numeric(lm.t$coefficients[2]),
                     beta.q=as.numeric(coefficients(lm.q)[2]),
                     Skill_Pearson=corr1,
                     rmse_t=sqrt(sum(res.t^2)/length(res.t)),
                     mae_t=sum(abs(res.t))/length(res.t)
  ))
  
} # end of i (time)

######################################################
## Fig S16a-c
win.graph(90,30)
par(mfrow=c(1,3),mar=c(4,4,1,1))
# A biplot demonstartes the empirical relationship between the estimated and theoretical interaction Jacobians 
# at a representative time point in which the overall inference skill equal to the median value  
tind <- which.min(abs(cof[,'Skill_Pearson']-median(cof[,'Skill_Pearson'])))
x.t <- TJ.ls[[tind]]; y.t <- EJ.ls[[tind]]
lm.t=lm(unlist(c(y.t))~unlist(c(x.t)))
plot(unlist(c(y.t))~unlist(c(x.t)),cex=1.5,cex.lab=1.5,
     ylab='Estimated interaction strength',xlab='Theoretical interaction strength')
abline(0,1,col='grey',lty=1);
abline(lm.t$coefficients,col=1,lty=1);
legend('bottomright',paste('r=',round(cof[tind,'Skill_Pearson'],3),
                           '\t t=',tind)
       ,bty='n')


# Summary statistics
cbind(mean=apply(cof,2,mean),sd=apply(cof,2,sd),median=apply(cof,2,median))

# Empirical distribution of overall inference skills derived from all time points
x <- cof[,'Skill_Pearson']

hist(cof[,'Skill_Pearson'],xlim=c(0.0,1.0),breaks=seq(0,1,0.05),border='white',
     main="",xlab="Overall inference skill (Pearson r)")
abline(v=mean(cof[,'Skill_Pearson']))
abline(v=median(cof[,'Skill_Pearson']),lty=2)
legend('topleft',paste(c(paste('mean =',round(mean(cof[,'Skill_Pearson']),3))),
                       '\nmedian =',round(median(cof[,'Skill_Pearson']),3)),
       bty='n'
)

# Node inference skills
y <- apply(cor.sp,2,mean)
x <- apply(cor.sp.o,2,mean)
boxplot(cbind(Outward=x,Inward=y),ylim=c(-0.1,1))
legend('bottomleft',paste(c(paste('mean =',round(mean(x),3))),
                          '\nmedian =',round(median(x),3)),
       bty='n'
)
legend('bottomright',paste(c(paste('mean =',round(mean(y),3))),
                           '\nmedian =',round(median(y),3)),
       bty='n'
)
