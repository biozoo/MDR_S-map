# Article: Reconstructing large interaction networks from empirical time series data 
# Authors: Chun-Wei Chang, Takeshi Miki, Masayuki Ushio, Po-Ju Ke, Hsiao-Pei Lu, Fuh-Kwo Shiah, and Chih-hao Hsieh
# Published in: 

#### Multi-species host-parasitoid Nicholson-Bailey model##############
rm(list = ls())
Root_Path <- 'D:\\data\\Meta_analysis\\AdpSmap\\Code'
Sub_Path <- '\\Demo_R_Code_MDR_S-MAP\\Additional models\\'
OPath <- paste(Root_Path,Sub_Path,sep='')
setwd(OPath)
SaveFile=F #T/F for saving files

#initialize random seed
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

#Define parameters for host parasitoid model, generation-non-overlapping
host_sp_no <- 200 #number of hosts = number of parasitoid
parasitoid_sp_no <- 400
sp_no <- host_sp_no + parasitoid_sp_no #total species number
nv <- numeric(sp_no)  #true value
initial_nv <- numeric(sp_no) #common initial values

r <- numeric(host_sp_no) #growth rate of hosts
aii <- numeric(host_sp_no) #self regulation of hosts

####Generate growth rate r####
max_var_r <- 0.8
ave_r <- 2.0
const_aii <- 0.5
for(j in 1:host_sp_no) {
  r[j] = ave_r*(1.0 + runif(1, -max_var_r, max_var_r)) 
  aii[j] = const_aii*(1.0 + runif(1, -max_var_r, max_var_r))  #half of the parasitoid impact??
}

if(SaveFile){
  saveRDS(r, "./Parameters/growth_rHP20210611_NEW.rds")
  saveRDS(aii, "./Parameters/self_regHP20210611_NEW.rds")
}

###############Generate interaction matrix ##################
interaction_matrix <- matrix(0:0,nrow=host_sp_no, ncol=parasitoid_sp_no)
max_int<-100.0/host_sp_no  #This is a typical bifurcation parameter, shifting from equlibrium, periodic to chaos
#max_int = 50/200: equilibrium, 80/200: periodic, 100/200 chaotic, 400/200, crash of parasitoid
sdd<-2.0

for(k in 1:host_sp_no) {
  for(j in 1:parasitoid_sp_no) {
    interaction_matrix[k, j] = abs(truncated_norm(average=0, SD=sdd, thres=max_int))  #effect of sp. j on sp.k
  }
}

#For the main host-parasitoid pairs
for(k in 1:host_sp_no) {
  interaction_matrix[k, k] = 1.0*(1 + truncated_norm(average=0, SD=0.2, thres=0.9))
}

#save the result
if(SaveFile){saveRDS(interaction_matrix, "./Parameters/interactionHP20210611_NEW.rds")}

######The DE map of host-parasitoid model#######
HP_map <- function(invec, t, int_matrix, rr,aa)
{
  #for host i in 1:host_sp_no
  next_vec <- numeric(sp_no)
  parasitism_vec <- int_matrix%*%invec[(host_sp_no + 1):sp_no]
  next_vec[1: host_sp_no] <- rr[1:host_sp_no]*invec[1:host_sp_no]*exp(-aa[1: host_sp_no]*invec[1: host_sp_no] - parasitism_vec[1: host_sp_no]) #F_i = r_i*H_i*exp(-a_i*H_i - [cP]_i)
  for(j in 1: parasitoid_sp_no) {
    next_vec[host_sp_no + j] <- 0.0
    for(k in 1:host_sp_no) {
      if(parasitism_vec[k] > 1.0e-10) next_vec[host_sp_no + j] <-  next_vec[host_sp_no + j] + int_matrix[k,j]*(invec[host_sp_no + j]/parasitism_vec[k])*invec[k]*exp(-aa[k]*invec[k])*(1.0 - exp(-parasitism_vec[k]))
    }
  }
  return(next_vec)
}

######Generate initial condition (initial population size)#####
for(j in 1: host_sp_no) initial_nv[j] <- 1.0*(1 + truncated_norm(average=0, SD=0.2, thres=0.9))
for(j in (host_sp_no+1): sp_no) initial_nv[j] <- 0.1*(1 + truncated_norm(average=0, SD=0.2, thres=0.9))
if(SaveFile){saveRDS(initial_nv, "./Parameters/initial_nvHP20210611_NEW.rds")}

######Generate time series data#####
time <- 0  # initial condition (initial time, 0)
end_time <- 1000
nv <- initial_nv
result2<-t(rbind(time, as.data.frame(nv)))
for(time in 1:end_time)
{
  #nv <- Ricker_map(nv, time, interaction_matrix[[1]], r_env[,time])
  nv <- HP_map(nv, time, interaction_matrix, r, aii)
  cat(time,"\n")
  result2_2<-t(rbind(time, as.data.frame(nv)))
  result2<-rbind(result2, result2_2)
}

host_ID <- 4
plot(result2[1:1000,1], result2[1:1000,(1 + host_ID)], type="l")
result2 <- as.data.frame(result2) #conversion from matrix to dataframe
if(SaveFile){
  write.csv(result2,"./Data/HP_dynamics20210611_NEW.csv", row.names=FALSE)
}

###############The computations of interaction Jacobian matrices#####################
# Warning: It requires long computation time for deriving interaction Jacobian matrices for host-parasitoid model
#Define parameters for DE for host-parasitoid model, generation-non-overlapping
host_sp_no <- 200 #number of hosts = number of parasitoid
parasitoid_sp_no <- 400
sp_no <- host_sp_no + parasitoid_sp_no #total species number

r <- readRDS("./Parameters/growth_rHP20210611.rds") #growth rate of hosts
aii <- readRDS("./Parameters/self_regHP20210611.rds") #self regulation of hosts
interaction_matrix <- readRDS("./Parameters/interactionHP20210611.rds") #the interaction matrix C_ij
class(interaction_matrix)
result2 <- read.csv("./Data/HP_dynamics20210611.csv",header=T,stringsAsFactors = F) #time series + time
class(result2)
#check Jacobian matrix


#Function to calculate theoretical Jacobian matrix
HP_map_diff <- function(data_frame, t, int_matrix, rr, aa)
{
  invec <- numeric(sp_no)
  invec <- as.numeric(data_frame[t, -1]) #delete time information
  parasitism_vec <- int_matrix%*%invec[(host_sp_no + 1):sp_no] #inner product c and P
  diff_matrix <- matrix(0:0,nrow=sp_no, ncol=sp_no)
  
  #for dF_i/dH_i
  for(i in 1: host_sp_no) diff_matrix[i,i] <- rr[i]*exp(-aa[i]*invec[i] - parasitism_vec[i])*(1.0 - aa[i]*invec[i]) #dF_i/dH_i = r_i*exp(-a_i*H_i - [cP]_i)*(1 - a_i*H_i)
  #for dF_i/dP_j
  for(i in 1: host_sp_no) {
    for(j in 1:parasitoid_sp_no) diff_matrix[i, (host_sp_no+j)] <- -int_matrix[i,j]*rr[i]*invec[i]*exp(-aa[i]*invec[i] - parasitism_vec[i]) 
  } #dF_i/dP_j = -c_ij*r_i*H_i*exp(-a_i*H_i - [cP]_i)
  
  #for dG_j/dH_i
  for(j in 1:parasitoid_sp_no) {
    for(i in 1:host_sp_no)   diff_matrix[(host_sp_no+j), i] <- int_matrix[i,j]*(invec[host_sp_no + j]/parasitism_vec[i])*exp(-aa[i]*invec[i])*(1.0 - exp(-parasitism_vec[i]))*(1.0 - aa[i]*invec[i])
  } #for dG_j/dH_i = (c_ij*P_j/[cP]_i)*exp(-a_i*H_i)*(1 - exp(-[CP]_i))*(1.0 - a_i*H_i)
  
  #for dG_j/dP_y
  for(j in 1:parasitoid_sp_no) {
    for(k in 1:parasitoid_sp_no) {
      diff_matrix[(host_sp_no+j), (host_sp_no+k)] <- 0 #for making sure
      for(i in 1:host_sp_no){
        diff_matrix[(host_sp_no+j), (host_sp_no+k)]  <- diff_matrix[(host_sp_no+j), (host_sp_no+k)] + (int_matrix[i,j]*invec[(host_sp_no+j)]*invec[i]*exp(-aa[i]*invec[i])/(parasitism_vec[i]*parasitism_vec[i]))*(parasitism_vec[i]*int_matrix[i, k]*exp(-parasitism_vec[i]) - int_matrix[i,k]*(1.0 - exp(-parasitism_vec[i])))
      } #dG_j/dP_y(k) = (sum over i)c_ij*P_j*H_i*exp(-a_i*H_i)/([cP]_i^2) * ([cP]_i*c_ik*exp(-[cP]_i)-c_ik*(1- exp(-[cP]_i))) 
    }
  }
  #for dG_j/dP_j (overwriting of dGj/dPy)
  for(j in 1:parasitoid_sp_no) {
    for(i in 1:host_sp_no){
      diff_matrix[(host_sp_no+j), (host_sp_no+j)]  <- diff_matrix[(host_sp_no+j), (host_sp_no+j)]  + int_matrix[i,j]*invec[i]*exp(-aa[i]*invec[i])*(1.0 - exp(-parasitism_vec[i]))/parasitism_vec[i]
    } #dG_j/dP_j = (sum over i) c_ij*H_i*exp(-a_i*H_i)/([cP]_i^2) *([cP]_i*(1 + (c_ij*P_j*-1)*exp(-[cP]_i)) - c_ij*P_j*(1 - exp(-[cP]_i)))
  } 
  return(diff_matrix)
}


# Computation of interaction Jacobian matrices for the last 100 time points
# Warning: It requires long computation time for deriving interaction Jacobian matrices for host-parasitoid model
da.range=901:1000
Jaco_HP=list()
for(step_t in (da.range+1)){
  jaco.t <-HP_map_diff(result2, step_t, interaction_matrix, r, aii)  
  Jaco_HP[[step_t-901]]= jaco.t
}

if(SaveFile){saveRDS(Jaco_HP, "./Data/HP_Jaco.rds")}

###########################################################################
# Computation of interaction Jacobian matrices for the dominant species
(da.name='HP')
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



