# R code for comparing the reconstructed interaction Jacobian networks with their theoretical expectations 
# Article: Reconstructing large interaction networks from empirical time series data 
# Authors: Chun-Wei Chang, Takeshi Miki, Masayuki Ushio, Po-Ju Ke, Hsiao-Pei Lu, Fuh-Kwo Shiah, and Chih-hao Hsieh
# Published in: 

##############################################################################################
### Reconstructing interaction Jacobians for multi-species Ricker model 
rm(list = ls())
library(igraph)
library(dplyr)
library(quantreg)
library(vegan)

seed <- 49563
set.seed(seed)
# Set path
Root_Path <- 'D:\\data\\Meta_analysis\\AdpSmap\\Code'
Sub_Path <- '\\Demo_R_Code_MDR_S-MAP\\MainResults_RickerModels\\'
OPath <- paste(Root_Path,Sub_Path,sep='')
setwd(OPath)
source(paste(Root_Path,'\\Demo_R_Code_MDR_S-MAP\\Demo_MDR_function.R',sep=''))

# Dataset
do <- read.csv('./result20191024_0_0_0_.csv',header=T,stringsAsFactors = F)
da.range <- 101:200 # Subsample for data analysis
cv.unit <- 0.025
ptype <- 'aenet'
out.sample <- T # T/F for out-of-sample forecast
nout <- 2  # number of out-of-sample

dot <- do[da.range,1]
do <- do[da.range,-1]
dto <- apply(do,1,sum)
dpo <- do*repmat(dto^-1,1,ncol(do))

# Exclusion of rare species 
pcri <- 0;bcri <- 10^-3;
doind2 <- (apply(dpo,2,mean,na.rm=T)>(bcri))&((apply(do>0,2,sum,na.rm=T)/nrow(do))>pcri)
exsp2 <- setdiff(1:ncol(do),which(doind2))

do <- do[,-exsp2]
nsp <- ncol(do)
ndo <- nrow(do)
nin <- ndo-nout

# Mean and SD of each node
do.mean <- apply(do[1:nin,],2,mean,na.rm=T)
do.sd <- apply(do[1:nin,],2,sd,na.rm=T)
dosdM <- repmat(c(do.sd)^-1,1,nsp)*repmat(c(do.sd),nsp,1)


# In-sample
d <- do[1:(nin-1),]                          # In-sample dataset at time t
d_tp1 <- do[2:(nin),]                        # In-sample dataset at time t+1
ds <- (d-repmat(do.mean,nrow(d),1))*repmat(do.sd,nrow(d),1)^-1 # Normalized in-sample dataset at time t
ds_tp1 <- (d_tp1-repmat(do.mean,nrow(d_tp1),1))*repmat(do.sd,nrow(d_tp1),1)^-1 # Normalized in-sample dataset at time t+1

# Out-sample
d.test <- do[nin:(ndo-1),]                 # Out-of-sample dataset at time t 
dt_tp1 <- do[(nin+1):ndo,]                 # Out-of-sample dataset at time t+1
ds.test <- (d.test-repmat(do.mean,nrow(d.test),1))*repmat(do.sd,nrow(d.test),1)^-1 # Normalized out-of-sample dataset at time t
dst_tp1 <- (dt_tp1-repmat(do.mean,nrow(dt_tp1),1))*repmat(do.sd,nrow(dt_tp1),1)^-1 # Normalized out-of-sample dataset at time t+1

# Compiled data at time t 
ds.all <- rbind(ds,ds.test)

# Reconstructed interaction Jacobian networks by MDR S-map
ints <- read.csv(paste('model1024_0_0_0_nin98_cvunit0.025_aenet_jcof_Nmvx_Rallx_demo.csv',sep=''),stringsAsFactors = F)
# Qualitative network
intsp <- ints[,-c(1:3)]
intsp[ints[,-c(1:3)]!=0] <- 1
intsp <- cbind(ints[,c(1:3)],intsp)

######################
# Theoretical networks derived from multi-species Ricker model
jacobian_true <- read.csv('theoretical_DF20191024_demo.csv',header=T)


######################
# Time-varying networks
########################
# Time-varying networks
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
  cor.sp.t <- cor.sp.t.o <- NULL;for(j in 1:nsp){
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
cof[cof[,3]<0,3]=0

######################################################
## Fig 2
win.graph(90,90)
par(mfrow=c(2,2),mar=c(4,4,1,1))
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
                           '\nslope=',round(cof[tind,'beta'],3),
                           '\t t=',tind)
       ,bty='n')


# Summary statistics
cbind(mean=apply(cof,2,mean),sd=apply(cof,2,sd),median=apply(cof,2,median))

# Empirical distribution of overall inference skills derived from all time points
x <- cof[,'Skill_Pearson']

hist(x,xlim=c(0.0,1.0),breaks=seq(0,1,0.05),border='white',
     main="",xlab="Overall inference skill (Pearson r)")
abline(v=mean(x))
abline(v=median(x),lty=2)
legend('topleft',paste(c(paste('mean =',round(mean(x),3))),
                       '\nmedian =',round(median(x),3)),
       bty='n'
)

# Empirical distribution of regression slope
x <- cof[,'beta']
hist(x,xlim=c(0.0,1.0),breaks=seq(0,1,0.05),border='white',
     main="",xlab="Regression slope ")
abline(v=mean(x))
abline(v=median(x),lty=2)
legend('topleft',paste(c(paste('mean =',round(mean(x),3))),
                       '\nmedian =',round(median(x),3)),
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

############################################################
# Out-of-sample forecast for interaction networks (Fig. S3)
win.graph(40,20);par(mfrow=c(1,2),mar=c(4,4,2,1))
cof.o <- erro.o <- NULL
for(i in 1:length(TJ.ols)){
  x.t <- TJ.ols[[i]] # Theoretical network
  y.t <- EJ.ols[[i]] # Reconstructed network
  #erro.o=rbind(erro.o,erof(X=x.t,Y=y.t)$erro)
  
  # Statistical relationship between estimated interaction strength and theoretical expectation
  lm.t <- lm(unlist(c(y.t))~unlist(c(x.t)))
  lm.q <- rq(unlist(c(y.t))~unlist(c(x.t)))
  res.t <- unlist(c(x.t))-unlist(c(y.t))
  corr1 <- cor(x=unlist(c(x.t)), y=unlist(c(y.t)), method = 'pearson',use="pairwise.complete.obs")

  # Plots
  plot(unlist(c(y.t))~unlist(c(x.t)),cex=1.5,cex.lab=1.5,
       ylab='Estimated interaction strength',xlab='Theoretical interaction strength')
  abline(0,1,col='grey',lty=1);
  abline(lm.t$coefficients,col=1,lty=1);
  legend('topleft',
         paste('r=',round(cor(unlist(c(y.t)),unlist(c(x.t))),3),
               '\nslope=',round(lm.t$coefficients[2],3),
               '\t t= out',i)
         ,bty='n')
  
  # Output
  cof.o <- rbind(cof.o,c(beta=as.numeric(lm.t$coefficients[2]),
                      beta.q=as.numeric(coefficients(lm.q)[2]),
                      Skill_Pearson=corr1,
                      rmse_t=sqrt(sum(res.t^2)/length(res.t)),
                      mae_t=sum(abs(res.t))/length(res.t)
  ))
}

cbind(mean=apply(cof.o,2,mean),sd=apply(cof.o,2,sd),median=apply(cof.o,2,median))

###################################################################################################
# Time-integrated interaction networks reconstructed by MDR S-map or correlation analysis (Fig. S4)
win.graph(400,700);par(mfrow=c(2,1),mar=c(4,4,4,1))
# Take median values for each interaction strength estimated from different time points
funn <- 'median'
intsA.i <- matrix(apply(intsA,1,match.fun(funn)),nsp,nsp) # time-integrated reconstrcuted network 
JA2.i <- matrix(apply(JA2,1,match.fun(funn)),nsp,nsp)     # time-integrated theoretical network 

Y <- intsA.i;
X <- JA2.i;

# Statistical relationship between estimated interaction strength and theoretical expectation
lm.t <- lm(unlist(c(Y))~unlist(c(X)))
lm.q <- rq(unlist(c(Y))~unlist(c(X)))
# Prediction skills
corr1 <- cor(unlist(c(Y)),unlist(c(X)), method = 'pearson',use="pairwise.complete.obs")
res.t <- unlist(c(X))-unlist(c(Y))
rmse_i <- sqrt(sum(res.t^2,na.rm=T)/length(res.t))
mae_i <- sum(abs(res.t))/length(res.t)


#Results of longterm median for each interaction Jacobian
plot(unlist(c(Y))~unlist(c(X)),ylab='Time-integrated estimated strength',
     xlab='Time-integrated theoretical strength',
     main='Reconstructed network by MDR S-map',cex=1)
abline(0,1,col='grey');
abline(lm.t$coefficients,col=1);
legend(0,-0.5,
       paste('r=',round(corr1,3),
             '\nslope=',round(lm.t$coefficients[2],3)),
       bty='n',text.col=2)


#Results of correlation coefficients
corM <- cortex(do) #  computation the correlation coefficients between each pair of node
y.t <- corM$cor_mat*(corM$p_value<0.05); # Let correlation coefficients=0 if it is not significant
#Y=corM$cor_mat;                      # Or still showing those correlation coefficients even if no significance 
x.t <- JA2.i;
diag(x.t) <-  NA # exclude intraspecific interactions

lm.t <- lm(unlist(c(y.t))~unlist(c(x.t)))
lm.q <- rq(unlist(c(y.t))~unlist(c(x.t)))
# Prediction skills
res.t <- unlist(c(x.t))-unlist(c(y.t))
corr1 <- cor(unlist(c(y.t)),unlist(c(x.t)), method = 'pearson',use="pairwise.complete.obs")
cor.test(unlist(c(y.t)),unlist(c(x.t)), method = 'pearson',use="pairwise.complete.obs")

plot(unlist(c(y.t))~unlist(c(x.t)),ylab='Correlation coefficients',xlab='Time-integrated theoretical strength',
     main='Reconstructed network by correlation',cex=1)
abline(0,1,col='grey');
abline(lm.t$coefficients,col=1);
legend(0,-0.2,
       paste('r=',round(corr1,3),
             '\nslope=',round(lm.t$coefficients[2],3)),
       bty='n',text.col=2)

