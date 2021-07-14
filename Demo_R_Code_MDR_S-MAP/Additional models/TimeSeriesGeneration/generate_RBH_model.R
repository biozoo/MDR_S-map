# Article: Reconstructing large interaction networks from empirical time series data 
# Authors: Chun-Wei Chang, Takeshi Miki, Masayuki Ushio, Po-Ju Ke, Hsiao-Pei Lu, Fuh-Kwo Shiah, and Chih-hao Hsieh
# Published in: 

#####Ricker-Beverton-Holt##############
rm(list = ls())
Root_Path <- 'D:\\data\\Meta_analysis\\AdpSmap\\Code'
Sub_Path <- '\\Demo_R_Code_MDR_S-MAP\\Additional models\\'
OPath <- paste(Root_Path,Sub_Path,sep='')
setwd(OPath)
SaveFile=F #T/F for saving files

#Initialize random seed
set.seed(101)

##############################################
## Truncated normal distribution
truncated_norm <- function(average, SD, thres)
{
  temp <- rnorm(1, mean=average, sd=SD)
  if (abs(temp) > thres) {
    temp <- 0
  }
  return(temp)
}

#Define parameters for Ricker-Beverton-Holt model
sp_no <- 1000
nv <- numeric(sp_no)  #true value
nv2 <- numeric(sp_no) #observation with error
initial_nv <- numeric(sp_no) #common initial values

r <- numeric(sp_no) #intrinsic growth rate
c <- numeric(sp_no) #BH parameter


####code for generating growth rate r####
max_var_r <- 0.7
ave_r <- 10.0
const_c <- 1.0
for(j in 1:sp_no) {
  r[j] = ave_r*(1.0 + runif(1, -max_var_r, max_var_r)) 
  c[j] = const_c*(1.0 + runif(1, -max_var_r, max_var_r))  
}

if(SaveFile){
  saveRDS(r, "./Parameters/growth_rRBH20210611_NEW.rds")
  saveRDS(c, "./Parameters/growth_cRBH20210611_NEW.rds")
  }

  
  
###############generate interaction matrix for Ricker-Beverton-Holt##################
interaction_matrix <- matrix(0:0,nrow=sp_no, ncol=sp_no)
max_int <- 0.5 
sdd<-2.0

for(k in 1:sp_no) {
  for(j in 1:sp_no) {
    interaction_matrix[k, j] = abs(truncated_norm(average=0, SD=sdd, thres=max_int))  #effect of sp. j on sp.k
  }
}

#For the diagonal elements
a_ii <- 1.0
for(k in 1:sp_no) {
  interaction_matrix[k, k] = a_ii*(1 + truncated_norm(average=0, SD=0.2, thres=0.5))
}
#save the result
if(SaveFile){saveRDS(interaction_matrix, "./Parameters/interactionRBH20210621_NEW.rds")}

######The DE map of Ricker-Beverton-Holt#######
RBH_map <- function(invec, t, int_matrix, rr, cc)
{
  next_vec <- numeric(sp_no)
  competition_vec <- int_matrix%*%invec
  next_vec <- rr*invec*exp(-competition_vec)/(1 + cc*(1.0 - exp(-competition_vec)))
  return(next_vec)
}

######initial condition (initial population size)#####
for(j in 1: sp_no) {
  initial_nv[j] <- 0.0
  initial_nv[j] <- 0.1*(1 + truncated_norm(average=0, SD=0.2, thres=0.9))
}
if(SaveFile){saveRDS(initial_nv, "./Parameters/initial_nvRBH20210611_NEW.rds")}

######Generate time series data#######################
time <- 0  # initial condition (initial time, 0)
end_time <- 1000
obs_error <- 0.0
nv <- initial_nv
time <- 0
result2<-t(rbind(time, as.data.frame(nv)))
for(time in 1:end_time)
{
  
  nv <- RBH_map(nv, time, interaction_matrix, r, c)
  cat(time,"\n")
  result2_2<-t(rbind(time, as.data.frame(nv)))
  #result1<-rbind(result1, result1_2)
  result2<-rbind(result2, result2_2)
}

sp_ID <- 5
plot(result2[900:1000,1], result2[900:1000,(1 + sp_ID)], type="l")
if(SaveFile){write.csv(result2,"./Data/RBH_dynamics20210611_NEW.csv", row.names=FALSE)}




#####################################################################################
###############The computations of interaction Jacobian matrices#####################
sp_no <- 1000
r <- readRDS("./Parameters/growth_rRBH20210611.rds") #intrinsic growth rate
c <- readRDS("./Parameters/growth_cRBH20210611.rds") 
interaction_matrix <- readRDS("./Parameters/interactionRBH20210611.rds") #the interaction matrix C_ij
class(interaction_matrix)
result2 <- read.csv("./Data/RBH_dynamics20210611.csv",header=T,stringsAsFactors = F)
class(result2)
#check Jacobian matrix


#Function to calculate theoretical Jacobian matrix
RBH_map_diff <- function(data_frame, t, int_matrix, rr,cc)
{
  invec <- numeric(sp_no)
  invec <- as.numeric(data_frame[t, -1]) #delete time information
  competition_vec <- int_matrix%*%invec #inner product A and x
  diff_matrix <- matrix(0:0,nrow=sp_no, ncol=sp_no)
  
  for(i in 1:sp_no) {
    for(j in 1:sp_no) {
      diff_matrix[i,j] <- -int_matrix[i,j]*rr[i]*invec[i]*(1.0 + cc[i])*exp(-competition_vec[i])/((1.0 + cc[i]*(1.0 - exp(-competition_vec[i])))*(1.0 + cc[i]*(1.0 - exp(-competition_vec[i]))))
        
    }
  }
  for(i in 1:sp_no) diff_matrix[i,i] <- diff_matrix[i,i] + rr[i]*exp(-competition_vec[i])/(1.0 + cc[i]*(1.0 - exp(-competition_vec[i])))
  
  return(diff_matrix)
}

# Computation of interaction Jacobian matrices for the last 100 time points
da.range=901:1000
Jaco_RBH=list()
for(step_t in (da.range+1)){
 jaco.t <-RBH_map_diff(result2, step_t, interaction_matrix, r,c)  
 Jaco_RBH[[step_t-901]]= jaco.t
}

if(SaveFile){saveRDS(Jaco_RBH, "./Data/RBH_Jaco.rds")}

###########################################################################
# Computation of interaction Jacobian matrices for the dominant species
(da.name='RBH')
source(paste(Root_Path,'\\Demo_R_Code_MDR_S-MAP\\Demo_MDR_function.R',sep=''))
da.range=901:1000 # Subsample for data analysis
do=read.csv(paste('./Data/',da.name,'_dynamics20210611.csv',sep=''),header=T,stringsAsFactors = F)
do=do[do[,1]%in%da.range,]
dot=do[,1]
do=do[,-1]
colnames(do)=paste('X',1:ncol(do),sep='_')

# Excluding rare species with no clear temporal dynamics 
dto=apply(do,1,sum) # total abundance
dpo=do*repmat(dto^-1,1,ncol(do))
pcri=0;bcri=10^-3; # criteria for selecting species based on proportion of presence (pcri) and mean relative abundance (bri) 
doind2=(apply(dpo,2,mean,na.rm=T)>(bcri))&((apply(do>0,2,sum,na.rm=T)/nrow(do))>pcri) # index for selected species 
exsp2=setdiff(1:ncol(do),which(doind2))   # index for rare species 
do=do[,-exsp2]                            # Dataset excluded rare species
(nsp=ncol(do))                            # number of selected species

jacobian_ls <- readRDS(paste('./Data/',da.name,'_Jaco.rds',sep=''))
jacobian_true=NULL
for(i in 1:length(da.range)){
  j.t=jacobian_ls[[i]][-exsp2,][,-exsp2]
  jacobian_true=rbind(jacobian_true,data.frame(time=i,j.t))
}

if(SaveFile){write.csv(jacobian_true,paste('./Data/theoretical_DF_',da.name,'_demo_NEW.csv',sep=''),row.names = F)}



