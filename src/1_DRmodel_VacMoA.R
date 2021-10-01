############################################################
#Purpose: Dose-response model for unvaccinated and vaccinated
#Final edit: 1 Oct 2021
#Editor: Fumi Miura
############################################################
###Procedure
#0. install packages 
#1. Data table
#2. Dose-Response model
#3. Likelihood 
#4. MLE
#5. Summary of model fits
#6. Estimated parameters
############################################################

###0. install packages -----
install.packages("tidyverse")
library("tidyverse")

###1. Data table -----
#raw data
data <- read_csv("Data_HC_influenza_Oct2021.csv")
data$Dose <- as.numeric(data$Dose)
data$Total <- as.numeric(data$Total)
data$Infected <- as.numeric(data$Infected)
#all 
data_all <- data 
#all not-vaccinated
data_all_notvac <- data %>%
  filter(Vaccine == "No") %>%
  filter(Reassortant == "No") 
#all vaccinated
data_all_vac <- data %>%
  filter(Vaccine == "Yes") %>%
  filter(Reassortant == "No") 
#write rds 
write_rds(data_all, "data_all")

######2. DR model -----
###Unvaccinated dose response
#delta
DR_delta <- function(qd, para){
  a <- para[1]
  return(1-exp(-a*qd))
  } 
#gamma 
DR_gamma <- function(qd, para){
  theta <- para[1]
  k <- para[2]
  return(1-(1+theta*qd)^(-k)) 
  } 
#tweedie
DR_tweedie <- function(qd, para){  
  lambda <- para[1] 
  gamma <- para[2] 
  alpha <- para[3] 
  return(1 - exp(lambda*((1+gamma*qd)^(-alpha) - 1)))
  } 
#two-level 
DR_twolevel <- function(qd, para){
  a1 <- para[1]
  a2 <- para[2]
  p1 <- para[3]
  return(1- (p1*exp(-a1*qd) + (1-p1)*exp(-a2*qd))) 
  } 
#Gamma+pointmass
DR_gammaPointmass <- function(qd, para){
  a1 <- para[1]
  theta <- para[2]
  k <- para[3]
  p1 <- para[4]
  return(1- (p1*exp(-a1*qd) + (1-p1)*(1+theta*qd)^(-k)))
  }

###Vaccinated dose responses
#Gamma
AllNon_DR_gamma <- function(Veff, dose, para){
  return( (1-Veff)*DR_gamma(dose,para) )
}
Leaky_DR_gamma <- function(Veff, dose, para){
  return( DR_gamma((1-Veff)*dose,para) )
}
Mix_DR_gamma <- function(dose, para, Veff){
  alpha <- Veff[1]
  beta <- Veff[2]
  return((1-alpha)*DR_gamma((1-beta)*dose,para))
}

######3. Likelihood -----
###LogL for unvaccinated individuals
#delta
LogLikelihood_delta_e <- function(data){
  LL <- function(para){
    a <- exp(para[1])  
    sum( (-1)*(data[,2]-data[,3])*a*data[,1] + data[,3]*log( 1 - exp(-a*data[,1])))
    }
  return(LL)
} 
#gamma
LogLikelihood_gamma_e <- function(data){ 
  LL <- function(para){
    theta <- exp(para[1])
    k <- exp(para[2])
    sum( (-1)*(data[,2]-data[,3])*k*log(1 + theta*data[,1]) + data[,3]*log(1 - (1 + theta*data[,1])^(-k)) )
  }
  return(LL)
}
#tweedie
LogLikelihood_tweedie_e <- function(data){
  LL <- function(para){
    lambda <- exp(para[1]) 
    gamma <- exp(para[2]) 
    alpha <- exp(para[3]) 
    sum( (data[,2]-data[,3])*(  lambda*((1+gamma*data[,1])^(-alpha) - 1) ) 
         + data[,3]*log(1 - exp(lambda*((1+gamma*data[,1])^(-alpha) - 1))) )
  }
  return(LL)
}
#two level
LogLikelihood_twolevel_logistic <- function(data){
  LL <- function(para){
    a1 <- exp(para[1])
    a2 <- exp(para[2])
    p1 <- 1/(1+exp(-para[3]))
    sum( (data[,2]-data[,3])*log((p1*exp(-a1*data[,1]) + (1-p1)*exp(-a2*data[,1])) ) + data[,3]*log(1- ( p1*exp(-a1*data[,1]) + (1-p1)*exp(-a2*data[,1])) ) )
  }
  return(LL)
}
#Gamma+pointmass
LogLikelihood_gammaPointmass_logistic_e <- function(data){
  LL <- function(para){
    a1 <- exp(para[1])
    theta <- exp(para[2])
    k <- exp(para[3])
    p1 <- 1/(1+exp(-para[4]))
    sum( (data[,2]-data[,3])*log((p1*exp(-a1*data[,1]) + (1-p1)*(1+theta*data[,1])^(-k))) + data[,3]*log(1- (p1*exp(-a1*data[,1]) + (1-p1)*(1+theta*data[,1])^(-k))) )
  }
  return(LL)
}
#Gamma+pointmass_0immu
LogLikelihood_gammaPointmass_logistic_0immu_e <- function(data){
  LL <- function(para){
    theta <- exp(para[1])
    k <- exp(para[2])
    p1 <- 1/(1+exp(-para[3]))
    sum( (data[,2]-data[,3])*log((p1 + (1-p1)*(1+theta*data[,1])^(-k))) + data[,3]*log(1- (p1 + (1-p1)*(1+theta*data[,1])^(-k))) )
  }
  return(LL)
}
###LogL for vaccinated individuals
#Mix, where background is gamma distributed
LogLikelihood_mix_e <- function(data,para){ 
  LL <- function(Veff){ 
    theta <- para[1]
    k <- para[2]
    alpha <- 1/(1+exp(-Veff[1]))
    beta <- 1/(1+exp(-Veff[2]))
    sum( (data[,2]-data[,3])*log(1-(1-alpha)*(1-(1 + theta*(1-beta)*data[,1])^(-k))) + data[,3]*log((1-alpha)*(1-(1 + theta*(1-beta)*data[,1])^(-k))) )
  }
  return(LL)
} 
#All-or-nothing, where background is gamma distributed
LogLikelihood_AoN_e <- function(data,para){ 
  LL <- function(Veff){ 
    theta <- para[1]
    k <- para[2]
    alpha <- 1/(1+exp(-Veff[1]))
    beta <- 0
    sum( (data[,2]-data[,3])*log(1-(1-alpha)*(1-(1 + theta*(1-beta)*data[,1])^(-k))) + data[,3]*log((1-alpha)*(1-(1 + theta*(1-beta)*data[,1])^(-k))) )
  }
  return(LL)
} 
#Leaky, where background is gamma distributed
LogLikelihood_Leak_e <- function(data,para){ 
  LL <- function(Veff){ 
    theta <- para[1]
    k <- para[2]
    alpha <- 0
    beta <- 1/(1+exp(-Veff[1]))
    sum( (data[,2]-data[,3])*log(1-(1-alpha)*(1-(1 + theta*(1-beta)*data[,1])^(-k))) + data[,3]*log((1-alpha)*(1-(1 + theta*(1-beta)*data[,1])^(-k))) )
  }
  return(LL)
}

######4. MLE -----
#delta
ini_para_delta <- c(-8)
opt.delta_e <- list(
  optim( fn= LogLikelihood_delta_e(data_all_notvac), par=ini_para_delta, control = list(fnscale = -1), hessian = T),
  optim( fn= LogLikelihood_delta_e(data_all_vac),    par=ini_para_delta, control = list(fnscale = -1), hessian = T)
  )
#gamma
ini_para_gamma <- c(-1.5,-1.2)
opt.gamma_e <- list(
  optim( fn= LogLikelihood_gamma_e(data_all_notvac), par=ini_para_gamma, control = list(fnscale = -1), hessian = T),
  optim( fn= LogLikelihood_gamma_e(data_all_vac),    par=ini_para_gamma, control = list(fnscale = -1), hessian = T)
)
#tweedie
ini_para_twee <- c(log(10), log(5*10^(-5)), log(10))
opt.tweedie_e <- list(
  optim( fn= LogLikelihood_tweedie_e(data_all_notvac), par=ini_para_twee, control = list(fnscale = -1), hessian = T),
  optim( fn= LogLikelihood_tweedie_e(data_all_vac),    par=ini_para_twee, control = list(fnscale = -1), hessian = T)
  )
#two level
ini_para_twolevel <- c(log(10^(-8)),log(10^(-6)),-1)
opt.twolevel_logistic <- list(
  optim( fn= LogLikelihood_twolevel_logistic(data_all_notvac), par=ini_para_twolevel, control = list(fnscale = -1), hessian = T),
  optim( fn= LogLikelihood_twolevel_logistic(data_all_vac),    par=ini_para_twolevel, control = list(fnscale = -1), hessian = T)
  )
#Gamma+pointmass
ini_para_gamPM <- c(log(8*10^(-20)),log(5*10^(-6)), log(0.8), 0)
opt.gammaPointmass_logistic <- list(
  optim( fn= LogLikelihood_gammaPointmass_logistic_e(data_all_notvac), par=ini_para_gamPM, control = list(fnscale = -1), hessian = T),
  optim( fn= LogLikelihood_gammaPointmass_logistic_e(data_all_vac),    par=ini_para_gamPM, control = list(fnscale = -1), hessian = T)
  )
#Gamma+pointmass_0immu
ini_para_gamPM_0 <- c(log(5*10^(-6)), log(0.8), 0)
opt.gammaPointmass_logistic_0immu <- list(
  optim( fn= LogLikelihood_gammaPointmass_logistic_0immu_e(data_all_notvac), par=ini_para_gamPM_0, control = list(fnscale = -1), hessian = T),
  optim( fn= LogLikelihood_gammaPointmass_logistic_0immu_e(data_all_vac),    par=ini_para_gamPM_0, control = list(fnscale = -1), hessian = T)
  )
###optim() for proposed Leaky, AoN, and Mix models
nonvac_gamma_para <- exp(opt.gamma_e[[1]]$par) #set the parameters for susceptibility dist of the control group
opt.mix_gamma <- optim( fn= LogLikelihood_mix_e(data_all_vac,nonvac_gamma_para),    par=c(-5,-10), control = list(fnscale = -1), hessian = T)
opt.AoN_gamma <- optim( fn= LogLikelihood_AoN_e(data_all_vac,nonvac_gamma_para),    par=c(0), control = list(fnscale = -1), hessian = T)
opt.Leak_gamma <- optim( fn= LogLikelihood_Leak_e(data_all_vac,nonvac_gamma_para),    par=c(-5), control = list(fnscale = -1), hessian = T)

######5. Summary of model fits -----
###log-likelihood and AICs
#non vaccinated, general PDFs
loglike_est_nonvac <- c(opt.delta_e[[1]]$value, opt.gamma_e[[1]]$value, opt.tweedie_e[[1]]$value, opt.twolevel_logistic[[1]]$value, opt.gammaPointmass_logistic[[1]]$value, opt.gammaPointmass_logistic_0immu[[1]]$value) 
num_para_nonvac <- c(1,2,3,3,4,3)
AIC_est_nonvac <- 2*num_para_nonvac - 2*loglike_est_nonvac
Sum_modelfit_nonvac <- as.data.frame(rbind(loglike_est_nonvac, num_para_nonvac, AIC_est_nonvac))
names(Sum_modelfit_nonvac) <- c("delta", "gamma", "tweedie", "two-level", "gamma+point-mass", "gamma+point-mass_imm0") 

#vaccinated, general PDFs
loglike_est_vac <- c(opt.delta_e[[2]]$value, opt.gamma_e[[2]]$value, opt.tweedie_e[[2]]$value, opt.twolevel_logistic[[2]]$value, opt.gammaPointmass_logistic[[2]]$value, opt.gammaPointmass_logistic_0immu[[2]]$value) 
num_para_vac <- c(1,2,3,3,4,3)
AIC_est_vac <- 2*num_para_vac - 2*loglike_est_vac
Sum_modelfit_vac <- as.data.frame(rbind(loglike_est_vac, num_para_vac, AIC_est_vac))
names(Sum_modelfit_vac) <- c("delta_v", "gamma_v", "tweedie_v", "two-level_v", "gamma+point-mass_v", "gamma+point-mass_imm0_v")

Model_comp <- cbind(Sum_modelfit_nonvac, Sum_modelfit_vac)
rownames(Model_comp) <- c("Log-Likelihood", "No. of parameter", "AIC")

#Add results of proposed mixture models, combine all
AIC_Veff_mix_gamma  <- 2*2 - 2*opt.mix_gamma$value
AIC_Veff_AoN_gamma  <- 2*1 - 2*opt.AoN_gamma$value
AIC_Veff_Leak_gamma <- 2*1 - 2*opt.Leak_gamma$value

Model_comp_all <- cbind(Model_comp,
                        c(opt.AoN_gamma$value, 1, AIC_Veff_AoN_gamma),
                        c(opt.Leak_gamma$value, 1, AIC_Veff_Leak_gamma),
                        c(opt.mix_gamma$value, 2, AIC_Veff_mix_gamma))
rownames(Model_comp_all) <- c("Log-Likelihood", "No. of parameter", "AIC")
colnames(Model_comp_all) <- c("Delta", "Gamma", "Tweedie", "Two-level", "Gamma+point-mass", "Gamma+point-mass at 0", "Delta_v", "Gamma_v", "Tweedie_v", "Two-level_v", "Gamma+point-mass_v", "Gamma+point-mass at 0_v", "All-or-Nothing", "Leaky", "Mixture")

#Table S1, Table S2
write.csv(Model_comp_all, "TableS1S2_VacMoA.csv")

######6. Estimated parameters -----
est_para_list <- list(
  #non vaccinated
  delta_nonvac   = list(est=exp(opt.delta_e[[1]]$par)),
  gamma_nonvac   = list(est=exp(opt.gamma_e[[1]]$par)),
  tweedie_nonvac = list(est=exp(opt.tweedie_e[[1]]$par)),
  twolev_nonvac  = list(est=c(exp(opt.twolevel_logistic[[1]]$par[1:2]),1/(1+exp(-opt.twolevel_logistic[[1]]$par[3])))),
  gamPM_nonvac   = list(est=c(exp(opt.gammaPointmass_logistic[[1]]$par[1:3]), 1/(1+exp(-opt.gammaPointmass_logistic[[1]]$par[4])))),
  gamPM0_nonvac  = list(est=c(exp(opt.gammaPointmass_logistic_0immu[[1]]$par[1:2]), 1/(1+exp(-opt.gammaPointmass_logistic_0immu[[1]]$par[3])))),
  #vaccinated
  tweedie_vac = list(est=exp(opt.tweedie_e[[2]]$par)),
  VacMix_gamma = list(est=1/(1+exp(-opt.mix_gamma$par))),
  VacAoN_gamma = list(est=1/(1+exp(-opt.AoN_gamma$par))),
  VacLeak_gamma = list(est=1/(1+exp(-opt.Leak_gamma$par)))
)

#write rds 
write_rds(est_para_list, "est_para_list")

#table 1 
#see "2_Bootstrap_VacMoA.R"