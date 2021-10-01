############################################################
#Purpose: ggplot codes for Vaccine MoA paper
#Final edit: 1 Oct 2021
#Editor: Fumi Miura
############################################################
###Procedure
#0. package 
#1. human challenge data and estimated parameters
#2. plots
############################################################

###0. package -----
library(ggplot2)
library(patchwork)
library(tweedie)
#Run "1_DRmodel_VacMoA.R"
#Run "2_Bootstrap_VacMoA.R"

###1.Human challenge data and estimated parameters -----
##human challenge data
data_plot <- dplyr::mutate(data_all, Prob = Infected/Total, LogDose = log10(Dose))
data_plot <- as.data.frame(data_plot)
data_fig1 <- as.data.frame(list(
  dose = c(c(10^0, 10^2, 10^5, 10^8), c(10^0, 10^2, 10^4, 10^8)),
  log10dose = c(c(0, 2, 5, 8), c(0, 2, 4, 8)),
  Total = c(c(16, 40, 30, 10), c(40, 40, 10, 30)),
  Infected = c(c(4, 30, 27, 10), c(0, 2, 4, 18)),
  Prob = c(c(4/16,30/40,27/30, 10/10), c(0/40, 2/40, 4/10, 20/30)),
  Vaccine = c(rep("No",4), rep("Yes",4))
))

###2. Plots -----
#####Fig1 -----
x1 <- 10^seq(-3,12,by=0.05)
gam_example    <- as.data.frame(list(dose=x1, log10dose=log10(x1), pred=DR_gamma(x1,c(5,0.2))))
delta_example  <- as.data.frame(list(dose=x1, log10dose=log10(x1), pred=DR_delta(x1,c(0.01))))
twolev_example <- as.data.frame(list(dose=x1, log10dose=log10(x1), pred=DR_twolevel(x1,c(10^-6,1.5,0.3))))
AoN_example    <- as.data.frame(list(dose=x1, log10dose=log10(x1), pred=Mix_DR_gamma(x1,c(5,0.2), c(0.481,0))))
Leak_example   <- as.data.frame(list(dose=x1, log10dose=log10(x1), pred=Mix_DR_gamma(x1,c(5,0.2), c(0,0.9998))))
Mix_example   <- as.data.frame(list(dose=x1, log10dose=log10(x1), pred=Mix_DR_gamma(x1,c(5,0.2), c(0.3,0.998))))

fig_1A_DRmodel <- ggplot() + 
  xlim(-3,12) +
  ylim(0,1) +
  geom_point(data=data_fig1 %>% filter(Vaccine=="No"), aes(x=log10dose, y=Prob, size=Total, color=Vaccine),alpha=0.4) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(11, "RdBu")[11],RColorBrewer::brewer.pal(11, "RdBu")[1])) +
  geom_line(data=gam_example, aes(x=log10dose,y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1) + 
  geom_line(data=delta_example, aes(x=log10dose,y=pred), color="black", size=1) + 
  geom_line(data=twolev_example, aes(x=log10dose,y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[10], size=1) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none") +
  labs(x="log10(Dose)", y = "Probability of infection")

z_f1 <- seq(0,10,length=10^4)
gamma_ex_d <- as.data.frame(list(z=z_f1, dens=dgamma(z_f1, shape = 5, scale = 0.2))) 
delta_ex_d <- as.data.frame(list(z=z_f1, dens=c(rep(0,989),rep(1, length(z_f1[990:1100])), rep(0,length(z_f1)-1100))))
twolev_ex_d <- as.data.frame(list(z=z_f1, 
                                  dens=c(rep(0.3, 100), rep(0,1400),
                                         rep(0.7, 100), rep(0,length(z_f1)-1600))))

fig_1A_gamma <- ggplot() +
  xlim(0,5) +
  ylim(0,1.5) + 
  geom_area(data = gamma_ex_d, aes(x=z,y=dens), color=RColorBrewer::brewer.pal(11, "RdBu")[11], fill=RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.5) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 24))+
  labs(x="", y = "Density")
  
fig_1A_delta <- ggplot() +
  xlim(0,5) +
  ylim(0,1.5) + 
  geom_area(data = delta_ex_d, aes(x=z,y=dens), color="black", fill="black", alpha=0.5) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 24))+
  labs(x="", y = "")

fig_1A_twolev <- ggplot() +
  xlim(0,5) +
  ylim(0,1.5) + 
  geom_area(data = twolev_ex_d, aes(x=z,y=dens), color=RColorBrewer::brewer.pal(11, "RdBu")[10], fill=RColorBrewer::brewer.pal(11, "RdBu")[10], alpha=0.5) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 24))+
  labs(x="Susceptibility", y = "")
fig_1A_panel <- (fig_1A_delta / fig_1A_gamma / fig_1A_twolev)
fig_1A <- (fig_1A_DRmodel | fig_1A_panel)

fig_1B <- ggplot() + 
  xlim(-3,12) +
  ylim(0,1) +
  geom_point(data=data_fig1 %>% filter(Vaccine=="Yes"), aes(x=log10dose, y=Prob, size=Total, color=Vaccine),alpha=0.4) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(11, "RdBu")[1],RColorBrewer::brewer.pal(11, "RdBu")[11])) +
  geom_line(data=gam_example, aes(x=log10dose,y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1, alpha=0.4) + 
  geom_line(data=AoN_example, aes(x=log10dose,y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[4], size=1) +
  geom_line(data=Leak_example, aes(x=log10dose,y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[3], size=1) +
  geom_line(data=Mix_example, aes(x=log10dose,y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[2], size=1) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none") +
  labs(x="log10(Dose)", y = "Probability of infection")

(fig_1A | fig_1B) + plot_layout(widths = c(3, 1, 3)) + plot_annotation(tag_levels = c('A', rep('',3), 'B'))


#####Fig2 -----
x2 <- 10^seq(-3,12,by=0.05)
Nonvac_plot <- as.data.frame(list(dose=x2, log10dose=log10(x2), pred= DR_gamma(x2,est_para_list$gamma_nonvac$est), low = gamma_nonvac_traject$percentile_bt[1,], upp = gamma_nonvac_traject$percentile_bt[3,]))
AoN_vac_plot <- as.data.frame(list(dose=x2, log10dose=log10(x2), pred= Mix_DR_gamma(x2,est_para_list$gamma_nonvac$est, c(est_para_list$VacAoN_gamma$est,0)), low = AoN_gamma_vac_traject$percentile_bt[1,], upp = AoN_gamma_vac_traject$percentile_bt[3,]))
Leak_vac_plot <- as.data.frame(list(dose=x2, log10dose=log10(x2), pred= Mix_DR_gamma(x2,est_para_list$gamma_nonvac$est, c(0,est_para_list$VacLeak_gamma$est)), low = Leak_gamma_vac_traject$percentile_bt[1,], upp = Leak_gamma_vac_traject$percentile_bt[3,]))
Mix_vac_plot <- as.data.frame(list(dose=x2, log10dose=log10(x2), pred= Mix_DR_gamma(x2,est_para_list$gamma_nonvac$est, est_para_list$VacMix_gamma$est), low = mix_gamma_vac_traject$percentile_bt[1,], upp = mix_gamma_vac_traject$percentile_bt[3,]))
Gen_vac_plot <- as.data.frame(list(dose=x2, log10dose=log10(x2), pred=  DR_tweedie(x2,est_para_list$tweedie_vac$est), low = tweedie_vac_traject$percentile_bt[1,], upp = tweedie_vac_traject$percentile_bt[3,]))

fig_2A <- ggplot() + 
  xlim(-3,12) +
  ylim(0,1) +
  geom_point(data=data_plot, aes(x=LogDose,y=Prob, color=Vaccine, size=Total),alpha=0.4) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(11, "RdBu")[11],RColorBrewer::brewer.pal(11, "RdBu")[1])) +
  geom_line(data=Nonvac_plot, aes(x=log10dose, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1) +
  geom_ribbon(data=Nonvac_plot, aes(x=log10dose,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  geom_line(data=AoN_vac_plot, aes(x=log10dose, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[4], size=1) +
  geom_ribbon(data=AoN_vac_plot, aes(x=log10dose,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[4], alpha=0.1) +
  geom_line(data=Leak_vac_plot, aes(x=log10dose, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[3], size=1) +
  geom_ribbon(data=Leak_vac_plot, aes(x=log10dose,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[3], alpha=0.1) +
  geom_line(data=Mix_vac_plot, aes(x=log10dose, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[2], size=1,linetype="dashed") +
  geom_ribbon(data=Mix_vac_plot, aes(x=log10dose,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[2], alpha=0.1) +
  theme_classic(base_size = 24) +
  labs(x="log10(Dose)", y = "Probability of infection")

fig_2B <- ggplot() + 
  xlim(-3,12) +
  ylim(0,1) +
  geom_point(data=data_plot, aes(x=LogDose,y=Prob, color=Vaccine, size=Total),alpha=0.4) + 
  scale_color_manual(values = c(RColorBrewer::brewer.pal(11, "RdBu")[11],RColorBrewer::brewer.pal(11, "RdBu")[1])) +
  geom_line(data=Nonvac_plot, aes(x=log10dose, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=1) +
  geom_ribbon(data=Nonvac_plot, aes(x=log10dose,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  geom_line(data=Gen_vac_plot, aes(x=log10dose, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[1], size=1) +
  geom_ribbon(data=Gen_vac_plot, aes(x=log10dose,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[1], alpha=0.1) +
  theme_classic(base_size = 24) +
  labs(x="log10(Dose)", y = "Probability of infection")

###Fig-2(C), density
##Tweedie (vaccinated)
power <- (est_para_list$tweedie_vac$est[3]+2)/(est_para_list$tweedie_vac$est[3]+1)
mu <- est_para_list$tweedie_vac$est[1] * est_para_list$tweedie_vac$est[2] * est_para_list$tweedie_vac$est[3]
phi <-  (est_para_list$tweedie_vac$est[1]^(1-power)) * ((est_para_list$tweedie_vac$est[2] * est_para_list$tweedie_vac$est[3])^(2-power)) / (2-power)

##random variable 
z <- seq(0,10,length=10^4)
z_l <- seq(0,10*(1-est_para_list$VacLeak_gamma$est),length=10^4) #scaled by Veff_leak
z_twe <- seq(0,10*(0.001332391/0.9162384),length=10^4) #scaled by mean(Tweedie_nonvac) / mean(Gamma_nonvac)

##density
nonvac_gam_dens <- as.data.frame(list(
  susceptibility = z,
  density = dgamma(z,shape = est_para_list$gamma_nonvac$est[1], scale = est_para_list$gamma_nonvac$est[2])
  ))
AoN_gam_dens <- as.data.frame(list(
  susceptibility = z,
  density = c(rep(est_para_list$VacAoN_gamma$est,50), rep(0,(length(z)-50))) + (1-est_para_list$VacAoN_gamma$est) * dgamma(z,shape =  est_para_list$gamma_nonvac$est[1], scale =  est_para_list$gamma_nonvac$est[2])
))
Leak_gam_dens <- as.data.frame(list(
  susceptibility = z_l,
  density = (1/(1-est_para_list$VacLeak_gamma$est)) * dgamma(z, shape =  est_para_list$gamma_nonvac$est[1], scale =  est_para_list$gamma_nonvac$est[2])
))
Mix_gam_dens <- as.data.frame(list(
  susceptibility = z_l,
  density = c(rep(est_para_list$VacAoN_gamma$est,50), rep(0,(length(z)-50))) + (1-est_para_list$VacAoN_gamma$est) * (1/(1-est_para_list$VacLeak_gamma$est)) * dgamma(z, shape =  est_para_list$gamma_nonvac$est[1], scale =  est_para_list$gamma_nonvac$est[2])
))
GenVac_Twee_dens <- as.data.frame(list(
  susceptibility = z_twe,
  density = dtweedie(z_twe, power=power, mu=mu, phi=phi)
))

##ggplot
#nonvac_gam_dens, unvaccinated
fig_2C1 <- ggplot() + 
  xlim(0,max(nonvac_gam_dens$susceptibility)) + 
  ylim(0,max(nonvac_gam_dens$density)) + 
  geom_area(data = nonvac_gam_dens, aes(x=susceptibility,y=density), color=RColorBrewer::brewer.pal(11, "RdBu")[11], fill=RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.5) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 24))+
  labs(x="Susceptibility", y = "Density") +
  annotate("text",  x=Inf, y = Inf, label = "Non-vaccinated", vjust=1, hjust=1, color=RColorBrewer::brewer.pal(11, "RdBu")[11], size=6)
#AoN_gam_dens, All-or-nothing
fig_2C2 <- ggplot() +
  xlim(0,max(AoN_gam_dens$susceptibility)) + 
  ylim(0,max(AoN_gam_dens$density)) + 
  geom_area(data = AoN_gam_dens, aes(x=susceptibility,y=density), color=RColorBrewer::brewer.pal(11, "RdBu")[4], fill=RColorBrewer::brewer.pal(11, "RdBu")[4], alpha=0.5) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 24))+
  labs(x="Susceptibility", y = "Density") +
  annotate("text",  x=Inf, y = Inf, label = "All-or-Nothing", vjust=1, hjust=1, color=RColorBrewer::brewer.pal(11, "RdBu")[4], size=6)
#Leak_gam_dens, Leaky
fig_2C3 <- ggplot() + 
  xlim(0,max(Leak_gam_dens$susceptibility)) + 
  ylim(0,max(Leak_gam_dens$density)) + 
  geom_area(data = Leak_gam_dens, aes(x=susceptibility,y=density), color=RColorBrewer::brewer.pal(11, "RdBu")[3], fill=RColorBrewer::brewer.pal(11, "RdBu")[3], alpha=0.5) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 24))+
  labs(x="Susceptibility", y = "Density") +
  annotate("text",  x=Inf, y = Inf, label = "Leaky", vjust=1, hjust=1, color=RColorBrewer::brewer.pal(11, "RdBu")[3], size=6)
#Mix_gam_dens, Mixture 
fig_2C4 <- ggplot() +
  xlim(0,max(Mix_gam_dens$susceptibility)) + 
  ylim(0,max(Mix_gam_dens$density)) + 
  geom_area(data = Mix_gam_dens, aes(x=susceptibility,y=density), color=RColorBrewer::brewer.pal(11, "RdBu")[2], fill=RColorBrewer::brewer.pal(11, "RdBu")[2], alpha=0.5) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 24))+
  labs(x="Susceptibility", y = "Density") +
  annotate("text",  x=Inf, y = Inf, label = "Mixture", vjust=1, hjust=1, color=RColorBrewer::brewer.pal(11, "RdBu")[2], size=6)
#GenVac_Twee_dens, Generalized vaccinated model with Tweedie
fig_2C5 <- ggplot() +
  xlim(0,max(GenVac_Twee_dens$susceptibility)) + 
  ylim(0,max(GenVac_Twee_dens$density)) + 
  geom_area(data = GenVac_Twee_dens, aes(x=susceptibility,y=density), color=RColorBrewer::brewer.pal(11, "RdBu")[1], fill=RColorBrewer::brewer.pal(11, "RdBu")[1], alpha=0.5) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 24))+
  labs(x="Susceptibility", y = "Density") +
  annotate("text",  x=Inf, y = Inf, label = "Generalized", vjust=1, hjust=1, color=RColorBrewer::brewer.pal(11, "RdBu")[1], size=6)

fig_2C <- (fig_2C1 / fig_2C2 / fig_2C3 / fig_2C4 / fig_2C5) 
fig_2 <- (fig_2A / fig_2B) | fig_2C

fig_2 + plot_annotation(tag_levels = 'A')

#####Fig3 -----
#x values for plot (=dose)
x3 <- 10^seq(-3,12,by=0.05)
#estimated VE(d)
VE_AoN  <- 1-Mix_DR_gamma(x3,est_para_list$gamma_nonvac$est, c(est_para_list$VacAoN_gamma$est,0))/DR_gamma(x3,est_para_list$gamma_nonvac$est)
VE_Leak <- 1-Mix_DR_gamma(x3,est_para_list$gamma_nonvac$est, c(0,est_para_list$VacLeak_gamma$est))/DR_gamma(x3,est_para_list$gamma_nonvac$est)
VE_mix  <- 1-Mix_DR_gamma(x3,est_para_list$gamma_nonvac$est, est_para_list$VacMix_gamma$est)/DR_gamma(x3,est_para_list$gamma_nonvac$est)
VE_Gen  <- 1-DR_tweedie(x3,est_para_list$tweedie_vac$est)/DR_gamma(x3,est_para_list$gamma_nonvac$est)
VE_result <- as.data.frame(list(dose=x3,log10dose=log10(x3), VE_mix=VE_mix, VE_AoN=VE_AoN, VE_Leak=VE_Leak, VE_Gen=VE_Gen))
#add 95% boostrap CI
VE_result_Mix <- as.data.frame(list(dose=x3,log10dose=log10(x3), 
                                    est=VE_mix, 
                                    low=VE_mix_gamma_vac$percentile_bt[1,], 
                                    upp=VE_mix_gamma_vac$percentile_bt[3,]))
VE_result_AoN <- as.data.frame(list(dose=x3,log10dose=log10(x3), 
                                    est=VE_AoN, 
                                    low=VE_AoN_gamma_vac$percentile_bt[1,], 
                                    upp=VE_AoN_gamma_vac$percentile_bt[3,]))
VE_result_Leak <- as.data.frame(list(dose=x3,log10dose=log10(x3), 
                                    est=VE_Leak, 
                                    low=VE_Leak_gamma_vac$percentile_bt[1,], 
                                    upp=VE_Leak_gamma_vac$percentile_bt[3,]))
VE_result_Gen <- as.data.frame(list(dose=x3,log10dose=log10(x3), 
                                     est=VE_Gen, 
                                     low=VE_tweedie_vac$percentile_bt[1,], 
                                     upp=VE_tweedie_vac$percentile_bt[3,]))
#replace 0> and 1< values with 0 and 1 for plot
VE_result_Mix$low[VE_result_Mix$low<0] <- 0
VE_result_Mix$upp[VE_result_Mix$upp>1] <- 1
VE_result_AoN$low[VE_result_AoN$low<0] <- 0
VE_result_AoN$upp[VE_result_AoN$upp>1] <- 1
VE_result_Leak$low[VE_result_Leak$low<0] <- 0
VE_result_Leak$upp[VE_result_Leak$upp>1] <- 1
VE_result_Gen$low[VE_result_Gen$low<0] <- 0
VE_result_Gen$upp[VE_result_Gen$upp>1] <- 1

#plot
fig_3A <- ggplot() +
  xlim(-2,max(VE_result_AoN$log10dose)) + 
  ylim(0,1) + 
  geom_ribbon(data = VE_result_AoN, aes(x=log10dose,ymin=low,ymax=upp), fill = "grey30", alpha=0.1) +
  geom_line(data = VE_result_AoN, aes(x=log10dose,y=VE_AoN), color=RColorBrewer::brewer.pal(11, "RdBu")[4], alpha=0.75, size=1.1) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50") , text = element_text(size = 24))+
  labs(x="log10(Dose)", y = "Vaccine efficacy") +
  annotate("text",  x=Inf, y = Inf, label = "All-or-Nothing", vjust=1, hjust=1, color=RColorBrewer::brewer.pal(11, "RdBu")[4], size=6)

fig_3B <- ggplot() +
  xlim(-2,max(VE_result_Leak$log10dose)) + 
  ylim(0,1) + 
  geom_ribbon(data = VE_result_Leak, aes(x=log10dose,ymin=low,ymax=upp), fill = "grey30", alpha=0.1) +
  geom_line(data = VE_result_Leak, aes(x=log10dose,y=VE_Leak), color=RColorBrewer::brewer.pal(11, "RdBu")[3], alpha=0.75, size=1.1) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 24))+
  labs(x="log10(Dose)", y = "Vaccine efficacy") +
  annotate("text",  x=Inf, y = Inf, label = "Leaky", vjust=1, hjust=1, color=RColorBrewer::brewer.pal(11, "RdBu")[3], size=6)

fig_3C1 <- ggplot() +
  xlim(-2,max(VE_result_Mix$log10dose)) + 
  ylim(0,1) + 
  geom_ribbon(data = VE_result_Mix, aes(x=log10dose,ymin=low,ymax=upp), fill = "grey30", alpha=0.1) +
  geom_line(data = VE_result_Mix, aes(x=log10dose,y=VE_mix), color=RColorBrewer::brewer.pal(11, "RdBu")[2], alpha=0.75, size=1.1) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 24))+
  labs(x="log10(Dose)", y = "Vaccine efficacy") +
  annotate("text",  x=Inf, y = Inf, label = "Mixture", vjust=1, hjust=1, color=RColorBrewer::brewer.pal(11, "RdBu")[2], size=6) 

fig_3C2 <- ggplot() +
  xlim(-2,max(VE_result_Mix$log10dose)) + 
  geom_line(data = VE_result_Mix, aes(x=log10dose,y=VE_mix), color=RColorBrewer::brewer.pal(11, "RdBu")[2], alpha=0.75, size=0.8) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), 
        axis.title.y = element_text(size = rel(0.75)), 
        axis.title.x = element_text(size = rel(0.75)),
        axis.text = element_text(size = rel(0.75)), 
        text = element_text(size = 18) )+
  scale_y_continuous(breaks=seq(0.48110,0.48118,by=0.00004), limits = c(0.48110,0.48118))+
  labs(x="log10(Dose)", y = "Efficacy") 

fig_3C <- fig_3C1 + inset_element(fig_3C2, left = 0.6, bottom = 0.01, right = 0.99, top = 0.4)

fig_3D <- ggplot() +
  xlim(-2,max(VE_result_Gen$log10dose)) + 
  ylim(0,1) + 
  geom_ribbon(data = VE_result_Gen, aes(x=log10dose,ymin=low,ymax=upp), fill = "grey30", alpha=0.1) +
  geom_line(data = VE_result_Gen, aes(x=log10dose,y=VE_Gen), color=RColorBrewer::brewer.pal(11, "RdBu")[1], alpha=0.75, size=1.1) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"), text = element_text(size = 24))+
  labs(x="log10(Dose)", y = "Vaccine efficacy") +
  annotate("text",  x=Inf, y = Inf, label = "Generalized", vjust=1, hjust=1, color=RColorBrewer::brewer.pal(11, "RdBu")[1], size=6)

fig_3 <- (fig_3A | fig_3C) / (fig_3B | fig_3D)

fig_3 + plot_annotation(tag_levels = 'A')