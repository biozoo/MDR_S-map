# The R code for generating time series data from the models with time-varying interaction coefficients
# Article: Reconstructing large interaction networks from empirical time series data 
# Authors: Chun-Wei Chang, Takeshi Miki, Masayuki Ushio, Po-Ju Ke, Hsiao-Pei Lu, Fuh-Kwo Shiah, and Chih-hao Hsieh
# Published in: 

rm(list = ls())
Root_Path <- 'D:\\data\\Meta_analysis\\AdpSmap\\Code'
Sub_Path <- '\\Demo_R_Code_MDR_S-MAP\\Additional models\\'
OPath <- paste(Root_Path,Sub_Path,sep='')
setwd(OPath)
SaveFile <- F #T/F for saving files
JacoCompu <- F #T/F for computaing interaction Jacobian matrices

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

################# Multi-species host-parasitoid Nicholson-Bailey model #########
host_sp_no <- 200 #number of hosts = number of parasitoid
parasitoid_sp_no <- 400
sp_no <- host_sp_no + parasitoid_sp_no #total species number
nv <- numeric(sp_no)  #true value
nv2 <- numeric(sp_no) #observation with error
initial_nv <- numeric(sp_no) #common initial values

####Loading parameters and initial conditions
r <- readRDS("./Parameters/growth_rHP20210611.rds") #growth rate of hosts
aii <- readRDS("./Parameters/self_regHP20210611.rds") #self regulation of hosts
interaction_matrix <- readRDS("./Parameters/interactionHP20210611.rds") #the interaction matrix C_ij
initial_nv <- readRDS("./Parameters/initial_nvHP20210611.rds") #initial conditions
dim(interaction_matrix)

###seasonality of interactions
#initialize random seed
set.seed(101)
min_s <- 0.75
max_s <- 1.25
spri <- matrix(1.0, host_sp_no, parasitoid_sp_no)
summ <- matrix(runif(host_sp_no*parasitoid_sp_no, min_s, max_s), host_sp_no, parasitoid_sp_no)
autu <- matrix(runif(host_sp_no*parasitoid_sp_no, min_s, max_s), host_sp_no, parasitoid_sp_no)
wint <- matrix(runif(host_sp_no*parasitoid_sp_no, min_s, max_s), host_sp_no, parasitoid_sp_no)

######The case for host-parasitoid map model
HP_map_s <- function(invec, t, int_matrix, rr, aa, spri_m, summ_m, autu_m, wint_m)
{
  #seasonality of interaction
  if((t%%12)%in%c(3:5)){int_matrix.cyc <- spri_m*int_matrix}
  if((t%%12)%in%c(6:8)){int_matrix.cyc <- summ_m*int_matrix}
  if((t%%12)%in%c(9:11)){int_matrix.cyc <- autu_m*int_matrix}
  if((t%%12)%in%c(0,1,2)){int_matrix.cyc <- wint_m*int_matrix}
  
  #for host i in 1:host_sp_no
  next_vec <- numeric(sp_no)
  parasitism_vec <- int_matrix.cyc%*%invec[(host_sp_no + 1):sp_no]
  
  next_vec[1: host_sp_no] <- rr[1:host_sp_no]*invec[1:host_sp_no]*exp(-aa[1: host_sp_no]*invec[1: host_sp_no] - parasitism_vec[1: host_sp_no]) #F_i = r_i*H_i*exp(-a_i*H_i - [cP]_i)
  for(j in 1: parasitoid_sp_no) {
    next_vec[host_sp_no + j] <- 0.0
    for(k in 1:host_sp_no) {
      if(parasitism_vec[k] > 1.0e-10) next_vec[host_sp_no + j] <-  next_vec[host_sp_no + j] + int_matrix.cyc[k,j]*(invec[host_sp_no + j]/parasitism_vec[k])*invec[k]*exp(-aa[k]*invec[k])*(1.0 - exp(-parasitism_vec[k]))
    } 
  }
  return(next_vec)
}

###############################################
#Generate time series data
time <- 0  # initial condition (initial time, 0)
end_time <- 1000
obs_error <- 0.0
nv <- initial_nv
time <- 0
result2<-t(rbind(time, as.data.frame(nv)))
for(time in 1:end_time)
{
  nv <- HP_map_s(nv, time, interaction_matrix, r, aii, spri, summ, autu, wint)
  if(time%%100==0){cat(time,"\n")}
  result2_2<-t(rbind(time, as.data.frame(nv)))
  result2<-rbind(result2, result2_2)
}

host_ID <- 4
plot(result2[1:1000,1], result2[1:1000,(1 + host_ID)], type="l")
parasitoid_ID <- 13
plot(result2[1:1000,1], result2[1:1000,(1 + parasitoid_ID + host_sp_no)], type="l")

result2 <- as.data.frame(result2) #conversion from matrix to dataframe
if(SaveFile){write.csv(result2,"./Data/HP_TV_dynamics20210622_NEW.csv", row.names=FALSE)}


###############The computations of interaction Jacobian matrices#####################
# Warning: It requires long computation time for deriving interaction Jacobian matrices for host-parasitoid model

#Function to calculate theoretical Jacobian matrix
HP_map_s_diff <- function(data_frame, t, int_matrix, rr, aa, spri_m, summ_m, autu_m, wint_m)
{
  invec <- numeric(sp_no)
  invec <- as.numeric(data_frame[t, -1]) #delete time information
  
  #seasonality of interaction
  if((t%%12)%in%c(3:5)){int_matrix.cyc <- spri_m*int_matrix}
  if((t%%12)%in%c(6:8)){int_matrix.cyc <- summ_m*int_matrix}
  if((t%%12)%in%c(9:11)){int_matrix.cyc <- autu_m*int_matrix}
  if((t%%12)%in%c(0,1,2)){int_matrix.cyc <- wint_m*int_matrix}
  
  parasitism_vec <- int_matrix.cyc%*%invec[(host_sp_no + 1):sp_no] #inner product c and P
  diff_matrix <- matrix(0:0,nrow=sp_no, ncol=sp_no)
  
  #for dF_i/dH_i
  for(i in 1: host_sp_no) diff_matrix[i,i] <- rr[i]*exp(-aa[i]*invec[i] - parasitism_vec[i])*(1.0 - aa[i]*invec[i]) #dF_i/dH_i = r_i*exp(-a_i*H_i - [cP]_i)*(1 - a_i*H_i)
  #for dF_i/dP_j
  for(i in 1: host_sp_no) {
    for(j in 1:parasitoid_sp_no) diff_matrix[i, (host_sp_no+j)] <- -int_matrix.cyc[i,j]*rr[i]*invec[i]*exp(-aa[i]*invec[i] - parasitism_vec[i]) 
  } #dF_i/dP_j = -c_ij*r_i*H_i*exp(-a_i*H_i - [cP]_i)
  
  #for dG_j/dH_i
  for(j in 1:parasitoid_sp_no) {
    for(i in 1:host_sp_no)   diff_matrix[(host_sp_no+j), i] <- int_matrix.cyc[i,j]*(invec[host_sp_no + j]/parasitism_vec[i])*exp(-aa[i]*invec[i])*(1.0 - exp(-parasitism_vec[i]))*(1.0 - aa[i]*invec[i])
  } #for dG_j/dH_i = (c_ij*P_j/[cP]_i)*exp(-a_i*H_i)*(1 - exp(-[CP]_i))*(1.0 - a_i*H_i)
  
  #for dG_j/dP_y
  for(j in 1:parasitoid_sp_no) {
    for(k in 1:parasitoid_sp_no) {
      diff_matrix[(host_sp_no+j), (host_sp_no+k)] <- 0 #for making sure
      for(i in 1:host_sp_no){
        diff_matrix[(host_sp_no+j), (host_sp_no+k)]  <- diff_matrix[(host_sp_no+j), (host_sp_no+k)] + (int_matrix.cyc[i,j]*invec[(host_sp_no+j)]*invec[i]*exp(-aa[i]*invec[i])/(parasitism_vec[i]*parasitism_vec[i]))*(parasitism_vec[i]*int_matrix.cyc[i, k]*exp(-parasitism_vec[i]) - int_matrix.cyc[i,k]*(1.0 - exp(-parasitism_vec[i])))
      } #dG_j/dP_y(k) = (sum over i)c_ij*P_j*H_i*exp(-a_i*H_i)/([cP]_i^2) * ([cP]_i*c_ik*exp(-[cP]_i)-c_ik*(1- exp(-[cP]_i))) 
    }
  }
  #for dG_j/dP_j (overwriting of dGj/dPy)
  for(j in 1:parasitoid_sp_no) {
    for(i in 1:host_sp_no){
      diff_matrix[(host_sp_no+j), (host_sp_no+j)]  <- diff_matrix[(host_sp_no+j), (host_sp_no+j)]  + int_matrix.cyc[i,j]*invec[i]*exp(-aa[i]*invec[i])*(1.0 - exp(-parasitism_vec[i]))/parasitism_vec[i]
    } #dG_j/dP_j = (sum over i) c_ij*H_i*exp(-a_i*H_i)/([cP]_i^2) *([cP]_i*(1 + (c_ij*P_j*-1)*exp(-[cP]_i)) - c_ij*P_j*(1 - exp(-[cP]_i)))
  } 
  return(diff_matrix)
}

step_t <- 950
test950<-LV_map_s_diff(result2, step_t, interaction_matrix, r, spri, summ, autu, wint) #taking  < 1 second
hist(log10(abs(test950)))

###########################################################################
# Computation of interaction Jacobian matrices for the last 100 time points

if(JacoCompu){
  da.range <- 901:1000
  Jaco_HP <- list()
  for(step_t in (da.range+1)){
    jaco.t <- HP_map_s_diff(result2, step_t, interaction_matrix, r, aii, spri, summ, autu, wint)  
    Jaco_HP[[step_t-901]] <- jaco.t
  }
  if(SaveFile){saveRDS(Jaco_HP, "./Data/HP_TV_Jaco.rds")}
}

###########################################################################
# Computation of interaction Jacobian matrices for the dominant species
(da.name='HP_TV')
source(paste(Root_Path,'\\Demo_R_Code_MDR_S-MAP\\Demo_MDR_function.R',sep=''))
da.range=901:1000 # Subsample for data analysis
do=read.csv(paste('./Data/',da.name,'_dynamics20210622.csv',sep=''),header=T,stringsAsFactors = F)
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


######################################################################
################# Multi-species discrete LV competition model#########
######################################################################
#define parameters for DE for LV competition
sp_no <- 1000  #total species number
nv <- numeric(sp_no)  #true value
nv2 <- numeric(sp_no) #observation with error
initial_nv <- numeric(sp_no) #common initial values

r <- numeric(sp_no) #intrinsic

####Loading parameters and initial conditions
r <- readRDS("./Parameters/growth_rLV20210611.rds") #intrinsic growth rate
interaction_matrix <- readRDS("./Parameters/interactionLV20210611.rds") #the interaction matrix C_ij
dim(interaction_matrix)
initial_nv <- readRDS("./Parameters/initial_nvLV20210611.rds") #initial conditions

###seasonality of interactions
#initialize random seed
set.seed(101)
min_s <- 0.5
max_s <- 1.00
spri <- matrix(1.0, sp_no, sp_no)
summ <- matrix(runif(sp_no*sp_no, min_s, max_s), sp_no, sp_no)
autu <- matrix(runif(sp_no*sp_no, min_s, max_s), sp_no, sp_no)
wint <- matrix(runif(sp_no*sp_no, min_s, max_s), sp_no, sp_no)
diag(summ) <- rep(1,sp_no)
diag(autu) <- rep(1,sp_no)
diag(wint) <- rep(1,sp_no)

######The case for LV map model
LV_s_map <- function(invec, t, int_matrix, rr, spri_m, summ_m, autu_m, wint_m)
{

  next_vec <- numeric(sp_no)
  #seasonality of interaction
  if((t%%12)%in%c(3:5)){int_matrix.cyc <- spri_m*int_matrix}
  if((t%%12)%in%c(6:8)){int_matrix.cyc <- summ_m*int_matrix}
  if((t%%12)%in%c(9:11)){int_matrix.cyc <- autu_m*int_matrix}
  if((t%%12)%in%c(0,1,2)){int_matrix.cyc <- wint_m*int_matrix}
  
  competition_vec <- int_matrix.cyc%*%invec
  
  next_vec <- invec*(rr - competition_vec)
  return(next_vec)
}


#Generating Time-series
time <- 0  # initial condition (initial time, 0)
end_time <- 1000
obs_error <- 0.0
nv <- initial_nv
time <- 0
result2<-t(rbind(time, as.data.frame(nv)))
for(time in 1:end_time)
{
  nv <- LV_s_map(nv, time, interaction_matrix, r, spri, summ, autu, wint)
  if(time%%100==0){cat(time,"\n")}
  result2_2<-t(rbind(time, as.data.frame(nv)))
  result2<-rbind(result2, result2_2)
}
#View(result2)
sp_ID <- 5
plot(result2[900:1000,1], result2[900:1000,(1 + sp_ID)], type="l")
result2 <- as.data.frame(result2) #conversion from matrix to dataframe
if(SaveFile){write.csv(result2,"./Data/LV_TV_dynamics20210622_NEW.csv", row.names=FALSE)}    

###############To Calculate theoretical Jacobian
#Function to calculate theoretical Jacobian matrix
LV_map_s_diff <- function(data_frame, t, int_matrix, rr, spri_m, summ_m, autu_m, wint_m)
{
  invec <- numeric(sp_no)
  invec <- as.numeric(data_frame[t, -1]) #delete time information
  
  #seasonality of interaction
  if((t%%12)%in%c(3:5)){int_matrix.cyc <- spri_m*int_matrix}
  if((t%%12)%in%c(6:8)){int_matrix.cyc <- summ_m*int_matrix}
  if((t%%12)%in%c(9:11)){int_matrix.cyc <- autu_m*int_matrix}
  if((t%%12)%in%c(0,1,2)){int_matrix.cyc <- wint_m*int_matrix}
  
  competition_vec <- int_matrix.cyc%*%invec #inner product A and x
  diff_matrix <- matrix(0:0,nrow=sp_no, ncol=sp_no)
  
  for(i in 1:sp_no) {
    for(j in 1:sp_no) {
      diff_matrix[i,j] <- -int_matrix.cyc[i,j]*invec[i]
    }
  }
  for(i in 1:sp_no) diff_matrix[i,i] <- diff_matrix[i,i] + rr[i] - competition_vec[i]
  
  return(diff_matrix)
}

step_t <- 950
test950<-LV_map_s_diff(result2, step_t, interaction_matrix, r, spri, summ, autu, wint) #taking  < 1 second
hist(log10(abs(test950)))

###########################################################################
# Computation of interaction Jacobian matrices for the last 100 time points
if(JacoCompu){
  da.range=901:1000
  Jaco_LV=list()
  for(step_t in (da.range+1)){
    jaco.t <-LV_map_s_diff(result2, step_t, interaction_matrix, r, spri, summ, autu, wint)  
    Jaco_LV[[step_t-901]]= jaco.t
  }
  
  if(SaveFile){saveRDS(Jaco_LV, "./Data/LV_TV_Jaco_NEW.rds")}
}

###########################################################################
# Computation of interaction Jacobian matrices for the dominant species
(da.name='LV_TV')
source(paste(Root_Path,'\\Demo_R_Code_MDR_S-MAP\\Demo_MDR_function.R',sep=''))
da.range=901:1000 # Subsample for data analysis
do=read.csv(paste('./Data/',da.name,'_dynamics20210622.csv',sep=''),header=T,stringsAsFactors = F)
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


###################################################################
################# Multi-species Ricker-Beverton-Holt model#########
###################################################################
#define parameters
sp_no <- 1000
nv <- numeric(sp_no)  #true value
nv2 <- numeric(sp_no) #observation with error
initial_nv <- numeric(sp_no) #common initial values

#load parameter values and initial conditions
r <- readRDS("./Parameters/growth_rRBH20210611.rds") #intrinsic growth rate
c <- readRDS("./Parameters/growth_cRBH20210611.rds") 
interaction_matrix <- readRDS("./Parameters/interactionRBH20210611.rds") #the interaction matrix C_ij
initial_nv <- readRDS("./Parameters/initial_nvRBH20210611.rds") #initial conditions


###seasonality of interactions
#initialize random seed
set.seed(101)
min_s <- 0.75
max_s <- 1.25
spri <- matrix(1.0, sp_no, sp_no)
summ <- matrix(runif(sp_no*sp_no, min_s, max_s), sp_no, sp_no)
autu <- matrix(runif(sp_no*sp_no, min_s, max_s), sp_no, sp_no)
wint <- matrix(runif(sp_no*sp_no, min_s, max_s), sp_no, sp_no)
diag(summ) <- rep(1,sp_no)
diag(autu) <- rep(1,sp_no)
diag(wint) <- rep(1,sp_no)

######The case for RBH map model
RBH_s_map <- function(invec, t, int_matrix, rr, cc, spri_m, summ_m, autu_m, wint_m)
{
  next_vec <- numeric(sp_no)

  #seasonality of interaction
  if((t%%12)%in%c(3:5)){int_matrix.cyc <- spri_m*int_matrix}
  if((t%%12)%in%c(6:8)){int_matrix.cyc <- summ_m*int_matrix}
  if((t%%12)%in%c(9:11)){int_matrix.cyc <- autu_m*int_matrix}
  if((t%%12)%in%c(0,1,2)){int_matrix.cyc <- wint_m*int_matrix}

  competition_vec <- int_matrix.cyc%*%invec
  
  next_vec <- rr*invec*exp(-competition_vec)/(1 + cc*(1.0 - exp(-competition_vec)))
  #for(i in 1: sp_no) next_vec[i] <- rr[i]*invec[i]*exp(-competition_vec[i])/(1 + cc[i]*(1.0 - exp(-competition_vec[i])))
  return(next_vec)
}


#Generating Time series
time <- 0  # initial condition (initial time, 0)
end_time <- 1000
nv <- initial_nv
time <- 0
result2<-t(rbind(time, as.data.frame(nv)))
for(time in 1:end_time)
{
  nv <- RBH_s_map(nv, time, interaction_matrix, r, c, spri, summ, autu, wint)
  if(time%%100==0){cat(time,"\n")}
  result2_2<-t(rbind(time, as.data.frame(nv)))
  result2<-rbind(result2, result2_2)
}
#View(result2)
plot(result2[1000,2:(sp_no+1)])
sp_ID <- 6
plot(result2[1:1000,1], result2[1:1000,(1 + sp_ID)], type="l")
result2 <- as.data.frame(result2) #conversion from matrix to dataframe
if(SaveFile){write.csv(result2,"./Data/RBH_TV_dynamics20210622_NEW.csv", row.names=FALSE)}

###############To calculate the theoretical Jacobian
#Function to calculate theoretical Jacobian matrix
RBH_map_s_diff <- function(data_frame, t, int_matrix, rr, cc, spri_m, summ_m, autu_m, wint_m)
{
  invec <- numeric(sp_no)
  invec <- as.numeric(data_frame[t, -1]) #delete time information
  
  #seasonality of interaction
  if((t%%12)%in%c(3:5)){int_matrix.cyc <- spri_m*int_matrix}
  if((t%%12)%in%c(6:8)){int_matrix.cyc <- summ_m*int_matrix}
  if((t%%12)%in%c(9:11)){int_matrix.cyc <- autu_m*int_matrix}
  if((t%%12)%in%c(0,1,2)){int_matrix.cyc <- wint_m*int_matrix}
  
  competition_vec <- int_matrix.cyc%*%invec #inner product A and x
  diff_matrix <- matrix(0:0,nrow=sp_no, ncol=sp_no)
  
  for(i in 1:sp_no) {
    for(j in 1:sp_no) {
      diff_matrix[i,j] <- -int_matrix.cyc[i,j]*rr[i]*invec[i]*(1.0 + cc[i])*exp(-competition_vec[i])/((1.0 + cc[i]*(1.0 - exp(-competition_vec[i])))*(1.0 + cc[i]*(1.0 - exp(-competition_vec[i]))))
      
    }
  }
  for(i in 1:sp_no) diff_matrix[i,i] <- diff_matrix[i,i] + rr[i]*exp(-competition_vec[i])/(1.0 + cc[i]*(1.0 - exp(-competition_vec[i])))
  
  return(diff_matrix)
}


step_t <- 950
test950<-RBH_map_s_diff(result2, step_t, interaction_matrix, r, c, spri, summ, autu, wint) #taking  < 1 second
hist(log10(abs(test950)))

###########################################################################
# Computation of interaction Jacobian matrices for the last 100 time points
if(JacoCompu){
  da.range <- 901:1000
  Jaco_RBH <- list()
  for(step_t in (da.range+1)){
    jaco.t <- RBH_map_s_diff(result2, step_t, interaction_matrix, r,c, spri, summ, autu, wint)  
    Jaco_RBH[[step_t-901]]= jaco.t
  }
  
  if(SaveFile){saveRDS(Jaco_RBH, "./Data/RBH_TV_Jaco.rds")}
}

###########################################################################
# Computation of interaction Jacobian matrices for the dominant species
(da.name='RBH_TV')
source(paste(Root_Path,'\\Demo_R_Code_MDR_S-MAP\\Demo_MDR_function.R',sep=''))
da.range=901:1000 # Subsample for data analysis
do=read.csv(paste('./Data/',da.name,'_dynamics20210622.csv',sep=''),header=T,stringsAsFactors = F)
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



################################################################################
###  Time-varying multi-species Ricker model
################################################################################
#define parameters for DE for LV competition
Root_Path <- 'D:\\data\\Meta_analysis\\AdpSmap\\Code'
Sub_Path <- '\\Demo_R_Code_MDR_S-MAP\\Additional models\\'
OPath <- paste(Root_Path,Sub_Path,sep='')
setwd(OPath)

sp_no <- 1000  #total species number
nv <- numeric(sp_no)  #true value
nv2 <- numeric(sp_no) #observation with error
initial_nv <- numeric(sp_no) #common initial values

r <- numeric(sp_no) #intrinsic
SaveFile <- F #T/F for saving files
JacoCompu <- F #T/F for computaing interaction Jacobian matrices

####Loading parameters and initial conditions
r <- readRDS("./Parameters/RK_r_env_base20191024.rds")[,1] #intrinsic growth rate
interaction_matrix <- readRDS("./Parameters/RK_interaction.rds")
#the interaction matrix C_ij
dim(interaction_matrix)
initial_nv <- readRDS("./Parameters/RK_initial_nv1024.rds") #initial conditions

###seasonality of interactions
#initialize random seed
set.seed(101)
min_s <- 0.75
max_s <- 1.25
spri <- matrix(1.0, sp_no, sp_no)
summ <- matrix(runif(sp_no*sp_no, min_s, max_s), sp_no, sp_no)
autu <- matrix(runif(sp_no*sp_no, min_s, max_s), sp_no, sp_no)
wint <- matrix(runif(sp_no*sp_no, min_s, max_s), sp_no, sp_no)
diag(summ) <- rep(1,sp_no)
diag(autu) <- rep(1,sp_no)
diag(wint) <- rep(1,sp_no)

######The case for Ricker map model
RK_s_map <- function(invec, t, int_matrix, rr, spri_m, summ_m, autu_m, wint_m)
{
  next_vec <- numeric(sp_no)
  #seasonality of interaction
  if((t%%12)%in%c(3:5)){int_matrix.cyc <- spri_m*int_matrix}
  if((t%%12)%in%c(6:8)){int_matrix.cyc <- summ_m*int_matrix}
  if((t%%12)%in%c(9:11)){int_matrix.cyc <- autu_m*int_matrix}
  if((t%%12)%in%c(0,1,2)){int_matrix.cyc <- wint_m*int_matrix}
  
  growth_vec <- int_matrix.cyc%*%invec
  next_vec <- invec*exp(rr*(1.0 + growth_vec))
  return(next_vec)
}


#Generating Time-series
time <- 0  # initial condition (initial time, 0)
end_time <- 1000
obs_error <- 0.0
nv <- initial_nv
time <- 0
result2<-t(rbind(time, as.data.frame(nv)))
for(time in 1:end_time)
{
  nv <- RK_s_map(nv, time, interaction_matrix, r, spri, summ, autu, wint)
  if(time%%100==0){cat(time,"\n")}
  result2_2 <- t(rbind(time, as.data.frame(nv)))
  result2 <- rbind(result2, result2_2)
}
#View(result2)
sp_ID <- 5
plot(result2[900:1000,1], result2[900:1000,(1 + sp_ID)], type="l")
result2 <- as.data.frame(result2) #conversion from matrix to dataframe
if(SaveFile){write.csv(result2,"./Data/RK_TV_dynamics20210622_NEW.csv", row.names=FALSE)}    

#Function to calculate theoretical Jacobian matrix
RK_map_s_diff <- function(data_frame, t, int_matrix, rr, spri_m, summ_m, autu_m, wint_m)
{
  invec <- numeric(sp_no)
  invec <- as.numeric(data_frame[t, -1]) #delete time information
  #seasonality of interaction
  if((t%%12)%in%c(3:5)){int_matrix.cyc <- spri_m*int_matrix}
  if((t%%12)%in%c(6:8)){int_matrix.cyc <- summ_m*int_matrix}
  if((t%%12)%in%c(9:11)){int_matrix.cyc <- autu_m*int_matrix}
  if((t%%12)%in%c(0,1,2)){int_matrix.cyc <- wint_m*int_matrix}
  
  growth_vec <- int_matrix.cyc%*%invec
  temp_vec <- as.vector(invec*exp(rr*(1.0 + growth_vec)))
  
  #differential coefifcient component (1) (all element)
  DF <- rr*(temp_vec*int_matrix.cyc)
  
  #differential coefficient component (2) (diagnal only)
  for(k in 1: length(int_matrix.cyc[1,])) {
    DF[k, k] <- DF[k, k] + exp(rr*(1.0 + growth_vec))[k,]
  }
  return(DF)
}

step_t <- 950
test950 <- RK_map_s_diff(result2, step_t, interaction_matrix, r, spri, summ, autu, wint) #taking  < 1 second
hist(log10(abs(test950)))

###########################################################################
# Computation of interaction Jacobian matrices for the last 100 time points
if(JacoCompu){
  da.range <- 901:1000
  Jaco_RK <- list()
  for(step_t in (da.range+1)){
    jaco.t <-RK_map_s_diff(result2, step_t, interaction_matrix, r, spri, summ, autu, wint)  
    Jaco_RK[[step_t-901]]= jaco.t
  }
  
  if(SaveFile){saveRDS(Jaco_RK, "./Data/RK_TV_Jaco.rds")}
}


###########################################################################
# Computation of interaction Jacobian matrices for the dominant species
(da.name='RK_TV')
source(paste(Root_Path,'\\Demo_R_Code_MDR_S-MAP\\Demo_MDR_function.R',sep=''))
da.range=901:1000 # Subsample for data analysis
do=read.csv(paste('./Data/',da.name,'_dynamics20210622.csv',sep=''),header=T,stringsAsFactors = F)
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

