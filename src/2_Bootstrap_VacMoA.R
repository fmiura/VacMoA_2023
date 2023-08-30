############################################################
#Purpose: Dose-response model for unvaccinated and vaccinated
#Final edit: 1 Oct 2021
#Editor: Fumi Miura
############################################################
###Procedure
#1. Run "1_DRmodel_VacMoA.R"
#2. Dose-Response model and bootstrap function
#3. Refitting 
#4. Draw trajectories of dose-response curves
#5. Draw dose-dependent VE
#6. Estimated parameters with 95%CI (Table-1)
############################################################

###1. Run "1_DRmodel_VacMoA.R" -----
#estimated parameters, likelihoods, dose-response models are required

###2. Dose-Response model and bootstrap function -----
##Dose-response model for vaccinated (Mixture model)
DR_mixPara_gamma <- function(qd, para){
  theta <- para[1]
  k <- para[2]
  alpha <- para[3]
  beta <- para[4]
  return((1-alpha)*(1-(1+theta*(1-beta)*qd)^(-k))) 
}
##Bootstrap function
Boot_func <- function(data, DRmodel, DRpara, N_boot=1000){
  Prob_est <- DRmodel(data$Dose, DRpara)
  Trial_num <- data$Total
  Boot_sample <- matrix(0,length(Trial_num),N_boot)
  for(i in 1:length(Trial_num)){
    Boot_sample[i,] <- rbinom(N_boot,Trial_num[i],Prob_est[i])
  }
  Boot_data <- matrix(0,length(Boot_sample[,1]),3)
  Boot_data[,1] <- data$Dose
  Boot_data[,2] <- data$Total
  
  return(list(Boot_sample = Boot_sample,
              Boot_data   = Boot_data))
}

###3. Refitting  -----
##Generate bootstrapped samples
Boot_gamma_nonvac   <- Boot_func(data_all_notvac, DR_gamma, est_para_list$gamma_nonvac$est, 1000)
Boot_mix_gamma_vac  <- Boot_func(data_all_vac,    DR_mixPara_gamma, c(est_para_list$gamma_nonvac$est, est_para_list$VacMix_gamma$est), 1000)
Boot_AoN_gamma_vac  <- Boot_func(data_all_vac,    DR_mixPara_gamma, c(est_para_list$gamma_nonvac$est, est_para_list$VacAoN_gamma$est, 0), 1000)
Boot_Leak_gamma_vac <- Boot_func(data_all_vac,    DR_mixPara_gamma, c(est_para_list$gamma_nonvac$est, 0, est_para_list$VacLeak_gamma$est), 1000)
Boot_tweedie_vac    <- Boot_func(data_all_vac,    DR_tweedie, est_para_list$tweedie_vac$est, 1000)

##gamma_nonvac -----
gamma_nonvac_refit <- list()
Boot_gamma_nonvac_para <- matrix(0,1000,2)
for(i in 1:1000){
  Boot_gamma_nonvac$Boot_data[,3] <- Boot_gamma_nonvac$Boot_sample[,i]
  gamma_nonvac_refit[i] <- optim( fn= LogLikelihood_gamma_e(Boot_gamma_nonvac$Boot_data), 
                               par=c(-6, 0), control = list(fnscale = -1), hessian = T)
  Boot_gamma_nonvac_para[i,] <- exp(gamma_nonvac_refit[[i]])
}

##mix_gamma_vac -----
mix_gamma_vac_refit <- list()
Boot_mix_gamma_vac_para <- matrix(0,1000,4)
for(i in 1:1000){
  Boot_mix_gamma_vac$Boot_data[,3] <- Boot_mix_gamma_vac$Boot_sample[,i]
  mix_gamma_vac_refit[i] <- optim( fn= LogLikelihood_mix_e(Boot_mix_gamma_vac$Boot_data, est_para_list$gamma_nonvac$est),
                                   par=c(-5,-10), control = list(fnscale = -1), hessian = T)
  Boot_mix_gamma_vac_para[i,] <- c(est_para_list$gamma_nonvac$est, 1/(1+exp(-mix_gamma_vac_refit[[i]])))
}

##AoN_gamma_vac -----
AoN_gamma_vac_refit <- list()
Boot_AoN_gamma_vac_para <- matrix(0,1000,4)
for(i in 1:1000){
  Boot_AoN_gamma_vac$Boot_data[,3] <- Boot_AoN_gamma_vac$Boot_sample[,i]
  AoN_gamma_vac_refit[i] <- optim( fn= LogLikelihood_AoN_e(Boot_AoN_gamma_vac$Boot_data, est_para_list$gamma_nonvac$est),
                                   par=c(0), control = list(fnscale = -1), hessian = T)
  Boot_AoN_gamma_vac_para[i,] <- c(est_para_list$gamma_nonvac$est, 1/(1+exp(-AoN_gamma_vac_refit[[i]])),0)
}

##Leak_gamma_vac -----
Leak_gamma_vac_refit <- list()
Boot_Leak_gamma_vac_para <- matrix(0,1000,4)
for(i in 1:1000){
  Boot_Leak_gamma_vac$Boot_data[,3] <- Boot_Leak_gamma_vac$Boot_sample[,i]
  Leak_gamma_vac_refit[i] <- optim( fn= LogLikelihood_Leak_e(Boot_Leak_gamma_vac$Boot_data, est_para_list$gamma_nonvac$est),
                                    par=c(-5), control = list(fnscale = -1), hessian = T)
  Boot_Leak_gamma_vac_para[i,] <- c(est_para_list$gamma_nonvac$est, 0, 1/(1+exp(-Leak_gamma_vac_refit[[i]])))
}

##tweedie_vac -----
tweedie_vac_refit <- list()
Boot_tweedie_vac_para <- matrix(0,1000,3)
for(i in 1:1000){
  Boot_tweedie_vac$Boot_data[,3] <- Boot_tweedie_vac$Boot_sample[,i]
  tweedie_vac_refit[i] <- optim( fn= LogLikelihood_tweedie_e(Boot_tweedie_vac$Boot_data),
                                 par=c(log(10), log(5*10^(-5)), log(10)), control = list(fnscale = -1), hessian = T)
  Boot_tweedie_vac_para[i,] <- exp(tweedie_vac_refit[[i]])
}

###4. Draw trajectories of dos-response curves -----
traject_func <- function(BootPara, DRmodel, x_plot){
  inter_bt <- matrix(0, length(BootPara[,1]), length(x_plot))
  for(i in 1:length(BootPara[,1])){
    inter_bt[i,] <- DRmodel(x_plot, BootPara[i,])
  }
  percentile_bt <- matrix(0, 3, length(x_plot))
  for(i in 1:length(x_plot)){
    percentile_bt[,i] <- quantile(inter_bt[,i], c(0.025, 0.5, 0.975))
  }
  return(list(inter_bt = inter_bt,
              percentile_bt = percentile_bt))
}

gamma_nonvac_traject <- traject_func(Boot_gamma_nonvac_para, DR_gamma, 10^seq(-3,12,by=0.05)) #note: x_plot should be matched with x3 in "3_ggplot_VacMoA.R". 
mix_gamma_vac_traject <- traject_func(Boot_mix_gamma_vac_para, DR_mixPara_gamma, 10^seq(-3,12,by=0.05))
AoN_gamma_vac_traject <- traject_func(Boot_AoN_gamma_vac_para, DR_mixPara_gamma, 10^seq(-3,12,by=0.05))
Leak_gamma_vac_traject <- traject_func(Boot_Leak_gamma_vac_para, DR_mixPara_gamma, 10^seq(-3,12,by=0.05))
tweedie_vac_traject <- traject_func(Boot_tweedie_vac_para, DR_tweedie, 10^seq(-3,12,by=0.05))

###5. Draw dose-dependent VE -----
VE_dose_func <- function(BootPara_nonvac, DRmodel_nonvac, BootPara_vac, DRmodel_vac, x_plot){
  VE_bt <- matrix(0, length(BootPara_nonvac[,1]), length(x_plot))
  inter_bt <- matrix(0, length(BootPara_nonvac[,1]), length(x_plot))
  for(i in 1:length(BootPara_nonvac[,1])){
    inter_bt[i,] <- 1 - (DRmodel_vac(x_plot, BootPara_vac[i,]) / DRmodel_nonvac(x_plot, BootPara_nonvac[i,]))
  }
  percentile_bt <- matrix(0, 3, length(x_plot))
  for(i in 1:length(x_plot)){
    percentile_bt[,i] <- quantile(inter_bt[,i], c(0.025, 0.5, 0.975))
  }
  return(list(inter_bt = inter_bt,
              percentile_bt = percentile_bt))
}

VE_mix_gamma_vac <- VE_dose_func(Boot_gamma_nonvac_para, DR_gamma, Boot_mix_gamma_vac_para, DR_mixPara_gamma, 10^seq(-3,12,by=0.05))
VE_AoN_gamma_vac <- VE_dose_func(Boot_gamma_nonvac_para, DR_gamma, Boot_AoN_gamma_vac_para, DR_mixPara_gamma, 10^seq(-3,12,by=0.05))
VE_Leak_gamma_vac <- VE_dose_func(Boot_gamma_nonvac_para, DR_gamma, Boot_Leak_gamma_vac_para, DR_mixPara_gamma, 10^seq(-3,12,by=0.05))
VE_tweedie_vac <-  VE_dose_func(Boot_gamma_nonvac_para, DR_gamma, Boot_tweedie_vac_para, DR_tweedie, 10^seq(-3,12,by=0.05))

###6. Estimated parameters with 95%CI (Table-1) -----
#gamma_nonvac
quantile(Boot_gamma_nonvac_para[,1], c(0.025, 0.975))
quantile(Boot_gamma_nonvac_para[,2], c(0.025, 0.975))
#mix_gamma_vac
quantile(Boot_mix_gamma_vac_para[,1], c(0.025, 0.975)) #fixed
quantile(Boot_mix_gamma_vac_para[,2], c(0.025, 0.975)) #fixed
quantile(Boot_mix_gamma_vac_para[,3], c(0.025, 0.975))
quantile(Boot_mix_gamma_vac_para[,4], c(0.025, 0.975))
#AoN_gamma_vac
quantile(Boot_AoN_gamma_vac_para[,1], c(0.025, 0.975)) #fixed
quantile(Boot_AoN_gamma_vac_para[,2], c(0.025, 0.975)) #fixed
quantile(Boot_AoN_gamma_vac_para[,3], c(0.025, 0.975))
quantile(Boot_AoN_gamma_vac_para[,4], c(0.025, 0.975)) #fixed as 0 
##Leak_gamma_vac
quantile(Boot_Leak_gamma_vac_para[,1], c(0.025, 0.975)) #fixed
quantile(Boot_Leak_gamma_vac_para[,2], c(0.025, 0.975)) #fixed
quantile(Boot_Leak_gamma_vac_para[,3], c(0.025, 0.975)) #fixed as 0
quantile(Boot_Leak_gamma_vac_para[,4], c(0.025, 0.975))
#tweedie_vac
quantile(Boot_tweedie_vac_para[,1], c(0.025, 0.975))
quantile(Boot_tweedie_vac_para[,2], c(0.025, 0.975))
quantile(Boot_tweedie_vac_para[,3], c(0.025, 0.975))
