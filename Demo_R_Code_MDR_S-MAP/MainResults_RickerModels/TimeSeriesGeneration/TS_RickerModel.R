###Multi-species Ricker Model##############
rm(list = ls())
Root_Path <- 'D:\\data\\Meta_analysis\\AdpSmap\\Code'
Sub_Path <- '\\Demo_R_Code_MDR_S-MAP\\MainResults_RickerModels\\'
OPath <- paste(Root_Path,Sub_Path,sep='')
setwd(OPath)
SaveFile <- F #T/F for saving files
library(dplyr)

######Iterative map of multi-species Ricker model#######
Ricker_map <- function(invec, t, int_matrix, rr)
{
  growth_vec <- int_matrix%*%invec
  next_vec <- invec*exp(rr*(1.0 + growth_vec))
  return(next_vec)
}

#Function to calculate theoretical Jacobian matrix
Ricker_map_diff <- function(invec, t, int_matrix, rr)
{
  growth_vec <- int_matrix%*%invec
  temp_vec <- as.vector(invec*exp(rr*(1.0 + growth_vec)))
  
  #differential coefifcient component (1) (all element)
  DF <- rr*(temp_vec*int_matrix)
  
  #differential coefficient component (2) (diagnal only)
  for(k in 1: length(int_matrix[1,])) {
    DF[k, k] = DF[k, k] + exp(rr*(1.0 + growth_vec))[k,]
  }
  return(DF)
}

#############Vector of growth rate ##############
r <- readRDS("./TimeSeriesGeneration/RK_r_env_base20191024.rds")[,1]
####initial condition (initial population size)##
initial_nv <- readRDS("./TimeSeriesGeneration/RK_initial_nv1024.rds")
####Loading interaction matrix###################
interaction_matrix <- readRDS("./TimeSeriesGeneration/RK_interaction.rds")
class(interaction_matrix)

######################The first step analysis with the first interaction matrix and the first set of noise series##############
end_time <- 1000
start <- 800
end <- 1000
sp_no <- 1000

#Initial settings
time <- 0
nv <- initial_nv
class(nv)
result2<-t(rbind(time, as.data.frame(nv)))
for(time in 1:end_time)
{
  nv <- Ricker_map(nv, time, interaction_matrix, r)
  result2_t<-t(rbind(time, as.data.frame(nv)))
  result2<-rbind(result2, result2_t)
}

colnames(result2) <- c('time',paste('V',1:1000,sep=''))
rownames(result2) <- NULL
#give names
result3 <- result2[start:end,]

# Because the multi-species Ricker model exhibits chaotic dynamics, 
# the derived values is sensitive to the rounding error that differs among various R versions or OS 
# Our computation for generating time series data, 'result20191024_0_0_0_.csv', was based on R version 3.6.3 performed in Linux workstation (x86_64-pc-linux-gnu (64-bit))
# Therefore, it is very likely to obtain time series values different from that recorded in the original dataset even if using the same parameter sets
# Nevertheless, the dynamical behavior among species remains similar because they are in the same strange attractor.
plot(V1~time,result3,type='o')

# To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
if(SaveFile){write.csv(result3, 'result20191024_0_0_0_NEW.csv', quote=FALSE, row.names=FALSE)}


###########################################################################
# Computation of interaction Jacobian matrices for the last 100 time points
#load time series 
result_forJ <- read.csv('result20191024_0_0_0_.csv', header=T)
da.range <- 101:200 # Subsample for data analysis
# Dominant species ID
sp_ID <- which(apply(result_forJ[da.range,-1]/apply(result_forJ[da.range,-1],1,sum),2,mean,na.rm=T)>10^-3)

jacobian_true <- NULL
#loop to generate interaction Jacobian at the last 100 time points
for(step in 899:998) {
  nvv <- as.numeric(filter(result_forJ, time==step)[,-1])
  theoretical_DF_raw <- Ricker_map_diff(nvv, step, interaction_matrix , r)
  jacobian_true <- rbind(jacobian_true,data.frame(time=step-898, theoretical_DF_raw[sp_ID, sp_ID])) #subset of DF with dominant species
  rm(theoretical_DF_raw)
}
colnames(jacobian_true) <- c('time',paste('X',1:length(sp_ID),sep=''))
# To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
if(SaveFile){write.csv(jacobian_true, 'theoretical_DF20191024_demo_NEW.csv', quote=FALSE, row.names=FALSE)}


