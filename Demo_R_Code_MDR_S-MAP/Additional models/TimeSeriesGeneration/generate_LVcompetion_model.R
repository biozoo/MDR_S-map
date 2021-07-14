# Article: Reconstructing large interaction networks from empirical time series data 
# Authors: Chun-Wei Chang, Takeshi Miki, Masayuki Ushio, Po-Ju Ke, Hsiao-Pei Lu, Fuh-Kwo Shiah, and Chih-hao Hsieh
# Published in: 

####Discrete LV competition model##############
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

#Define parameters for DE for LV competition
sp_no <- 1000  #total species number
nv <- numeric(sp_no)  #true value
initial_nv <- numeric(sp_no) #common initial values
r <- numeric(sp_no) #intrinsic

####code for generating growth rate r####
max_var_r <- 0.4
min_r <- 3.5

for(j in 1:sp_no) {
  r[j] = min_r +  runif(1, 0, max_var_r)
}

if(SaveFile){saveRDS(r, "./Parameters/growth_rLV20210611_NEW.rds")}

###############generate interaction matrix for competition model##################
interaction_matrix <- matrix(0:0,nrow=sp_no, ncol=sp_no)
const_aii <- min_r*0.5  #with this size, population never goes to negative w/o interactions
max_int<- 200.0/(sp_no*min_r)   #200 is oK 
sdd<-4.0

for(k in 1:sp_no) {
  for(j in 1:sp_no) {
    interaction_matrix[k, j] = abs(truncated_norm(average=0, SD=sdd, thres=max_int))  #effect of sp. j on sp.k
  }
}
hist(interaction_matrix)
sum(interaction_matrix > 0)/(sp_no*sp_no)

#For the intra-specific competition
for(k in 1:sp_no) {
  interaction_matrix[k, k] = const_aii*(1 + truncated_norm(average=0, SD=0.2, thres=0.3))
}

#save the result
if(SaveFile){saveRDS(interaction_matrix, "./Parameters/interactionLV20210611_NEW.rds")}


######The DE map of Discrete LV competition model#######
LV_map <- function(invec, t, int_matrix, rr)
{
  #for host i in 1:host_sp_no
  next_vec <- numeric(sp_no)
  competition_vec <- int_matrix%*%invec
  
  next_vec <- invec*(rr - competition_vec)
  return(next_vec)
}

######initial condition (initial population size)#####
for(j in 1: sp_no) initial_nv[j] <- 0.1*(1 + truncated_norm(average=0, SD=0.2, thres=0.99))
if(SaveFile){saveRDS(initial_nv, "./Parameters/initial_nvLV20210611_NEW.rds")}

######Generate time series data#####
end_time <- 1000
nv <- initial_nv
time <- 0
result2<-t(rbind(time, as.data.frame(nv)))
for(time in 1:end_time)
{
  nv <- LV_map(nv, time, interaction_matrix, r)
  cat(time,"\n")
  result2_2<-t(rbind(time, as.data.frame(nv)))
  result2<-rbind(result2, result2_2)
}

sp_ID <- 10
plot(result2[900:1000,1], result2[900:1000,(1 + sp_ID)], type="l")
if(SaveFile){write.csv(result2,"./Data/LV_dynamics20210611_NEW.csv", row.names=FALSE)}


#####################################################################################
###############The computations of interaction Jacobian matrices#####################
sp_no <- 1000
r <- readRDS("./Parameters/growth_rLV20210611.rds") #intrinsic growth rate
interaction_matrix <- readRDS("./Parameters/interactionLV20210611.rds") #the interaction matrix C_ij
class(interaction_matrix)
result2 <- read.csv("./Data/LV_dynamics20210611.csv",header=T,stringsAsFactors = F)
class(result2)

#Function to calculate theoretical Jacobian matrix
LV_map_diff <- function(data_frame, t, int_matrix, rr)
{
  invec <- numeric(sp_no)
  invec <- as.numeric(data_frame[t, -1]) #delete time information
  competition_vec <- int_matrix%*%invec #inner product A and x
  diff_matrix <- matrix(0:0,nrow=sp_no, ncol=sp_no)
  
  for(i in 1:sp_no) {
    for(j in 1:sp_no) {
      diff_matrix[i,j] <- -int_matrix[i,j]*invec[i]
    }
  }
  for(i in 1:sp_no) diff_matrix[i,i] <- diff_matrix[i,i] + rr[i] - competition_vec[i]
  
  return(diff_matrix)
}

# Computation of interaction Jacobian matrices for the last 100 time points
da.range <- 901:1000
Jaco_LV <- list()
for(step_t in (da.range+1)){
  jaco.t <-LV_map_diff(result2, step_t, interaction_matrix, r)  
  Jaco_LV[[step_t-901]] <- jaco.t
}

if(SaveFile){saveRDS(Jaco_LV, "./Data/LV_Jaco.rds")}


###########################################################################
# Computation of interaction Jacobian matrices for the dominant species
(da.name='LV')
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
pcri=0;bcri=1.2*10^-3; # criteria for selecting species based on proportion of presence (pcri) and mean relative abundance (bri) 
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

