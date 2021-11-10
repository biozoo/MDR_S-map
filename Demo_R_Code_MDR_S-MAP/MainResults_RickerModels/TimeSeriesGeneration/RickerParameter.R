# The R code for generating parameters in multi-species Ricker model
rm(list = ls())
Root_Path <- 'D:\\data\\Meta_analysis\\AdpSmap\\Code'
Sub_Path <- '\\Demo_R_Code_MDR_S-MAP\\MainResults_RickerModels\\'
OPath <- paste(Root_Path,Sub_Path,sep='')
setwd(OPath)
SaveFile <- F #T/F for saving files

truncated_norm <- function(average, SD, thres)
{
  temp <- rnorm(1, mean=average, sd=SD)
  if (abs(temp) > thres) {
    temp <- 0
  }
  return(temp)
}

#initialize random seed
set.seed(101)

#define parameters for ODE for prey predator model
sp_no <- 1000
nv <- numeric(sp_no)  #true value
nv2 <- numeric(sp_no) #observation with error
initial_nv <- numeric(sp_no) #common initial values

r <- numeric(sp_no)


####code for generating growth rate r####

max_var_r <- 0.5

for(j in 1:sp_no) {
  r[j] = 1.5*(1.0 + runif(1, -max_var_r, max_var_r)) 
}

if(SaveFile){saveRDS(r, "./result20190929/growth_r20190929.rds")}


###############generate auto-correlated time series of environmental parameter########
#initialize random seed for 20191024
set.seed(1024)

max_time <- 2000
env<-numeric(max_time)
cor <- 0.8
max_red <- 0.5
env[1]<-0.1

for(i in 2:max_time) env[i] = cor*env[i-1] + (1-cor)*runif(1, -max_red, max_red)
if(SaveFile){saveRDS(env, "./result20191024/env20191024.rds")}

#dependence of each species on environmental forcing
env_dep<-numeric(sp_no)
#runif(5, -max_dep, max_dep) 
env_dep <- runif(sp_no, -1.0, 1.0)
if(SaveFile){saveRDS(env_dep, "./result20191024/env_dep20191024.rds")}


###############generate interaction matrix##################

interaction_matrix <- matrix(0:0,nrow=sp_no, ncol=sp_no)
max_int<-0.5
sdd<-2.0

for(k in 1:(sp_no-1)) {
  for(j in (k+1):sp_no) {
    interaction_matrix[k, j] = truncated_norm(average=0, SD=sdd, thres=max_int)  #effect of sp. j on sp.k
  }
}
hist(interaction_matrix)
#hub species
#k = 5
#for(j in 1:k) {
#interaction_matrix[k, j] = truncated_norm(average=-0.2, SD=0.1, thres=0.5)

#}
for(k in 1:sp_no) {
  interaction_matrix[k, k] = -1.0
}
interaction_matrix[1,3]
for(k in 2:sp_no) {
  for(j in 1:(k-1)) {
    if(interaction_matrix[j, k] < 0) {
      interaction_matrix[k, j] = interaction_matrix[j, k]*(1.0 + 0.9*runif(1, -1, 1))  #competition
    }
    else {
      interaction_matrix[k, j] = -1*interaction_matrix[j, k]*(1.0 + 0.9*runif(1, -1, 1)) #trophic interaction, sp.j is prey while sp.k is predator
    }
  }
}
head(interaction_matrix)
#fraction of non_zero interactions
sum(abs(interaction_matrix) > 1.0e-2)/(2*sp_no*sp_no)
hist(interaction_matrix)
#save the result
saveRDS(interaction_matrix, "./result20190929/interaction20190929.rds")
readRDS("./result20190929/interaction20190929.rds")[1,2]
interaction_matrix[1,1:10]
interaction_matrix <- readRDS("./result20190929/interaction20190929.rds")
interaction_matrix2 <- readRDS("./result20191024/forClark20191024_revised_for_single_noise/interaction.rds")
interaction_matrix2[1,1:10]
#write.csv(result_filtered1, "result_01.csv", quote=FALSE, row.names=FALSE)
#interaction_matrix[1,]
#matrix_test <- readRDS("./result01_20190531/interaction.rds")
#matrix_test[1,]
exp(nv)
class(interaction_matrix)
interaction_matrix%*%nv


######The case for Ricker map model#######
Ricker_map <- function(invec, t, int_matrix, rr)
{
  growth_vec <- int_matrix%*%invec
  #This does not represent the correlation!!
  #Env_cor*Env[t] + (1 - Env_cor)*norm??
  next_vec <- invec*exp(rr*(1.0 + growth_vec))
  return(next_vec)
}


######initial condition (initial population size)#####
for(j in 1: sp_no) initial_nv[j] <- 0.1*(1 + truncated_norm(average=0, SD=0.2, thres=0.9))
saveRDS(initial_nv, "./result20190929/initial_nv.rds")
hist(initial_nv)
time <- 0  # initial condition (initial time, 0)
end_time <- 1000
#############generate growth rate matrix##############
r_env <- matrix(0:0,nrow=sp_no, ncol=end_time)

proc_error <- 0.0
forcing_size <- 0.1
for(k in 1: sp_no) {
  for(t in 1: end_time) {
    r_env[k, t] <- r[k]*(1.0 + forcing_size*env_dep[k]*env[t] + proc_error*runif(1, -1.0, 1.0))
  }
}
hist(r_env)
plot(c(1:1000), env[1:1000], type="l")
plot(c(1:1000), r_env[1,], type="l")
plot(env[1:1000], r_env[1,])
plot(env[1:1000], r_env[2,])
#save the result
#w proc_error <- 0, forcing_size <-0
saveRDS(r_env, "./result20190929/r_env_base20190929.rds")
saveRDS(r_env, "./result20191024/r_env_base20191024.rds")

#w proc_error <- 0.1, forcing_size <-0
saveRDS(r_env, "./result20190929/r_env_wProc20190929.rds")
saveRDS(r_env, "./result20191024/r_env_wProc20191024.rds")

#w proc_error <- 0.0, forcing_size <-0.1
saveRDS(r_env, "./result20190929/r_env_wForcing20190929.rds")
saveRDS(r_env, "./result20191024/r_env_wForcing20191024.rds")


#set common initial values
initial_nv <- readRDS("./result20190929/initial_nv.rds")
initial_nv <- as.numeric(readRDS("./result01_20190531/all_data_raw20190531.rds")[1,2:1001]) #for20191024
View(initial_nv)
saveRDS(initial_nv, "./result20191024/initial_nv1024.rds")

