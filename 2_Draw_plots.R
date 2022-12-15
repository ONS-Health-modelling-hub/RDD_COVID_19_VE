library(dplyr)
library(tidyr)
library(readr)
library(ggpubr)
library(ggplot2)


#\\nsdata2/MG
#setwd("T:/Branch folders/A&R/HALE/Vaccination/Effectiveness")

user ="Charlotte"
if (user== "Charlotte"){
setwd("C:/Users/bermic/Office for National Statistics/COVID-19 - Deaths_and_vaccination/RDD VE/plotting scripts updated")
}
if (user== "Vahé"){
  setwd("C:/Users/nafilv/Office for National Statistics(1)/COVID-19 - Deaths_and_vaccination/RDD VE/plotting scripts updated")
}



#-----------------------------------#
# Fig 1 - Cumulative incidence plot #
#-----------------------------------#
fit_df_age_group <- read_csv( "data\\cum_inc_vaccination_age_group.csv")


ggplot(data = fit_df_age_group, aes(x=t, prop_vaccinated, colour=age_group))+
  geom_step() + 
  labs(x="Days after 8 December 2020", 
       y="Proportion Received First Vaccination")+
  scale_colour_grey(name=expression(underline("Age group")),
                    labels=c("75-79", "80-84"))+
  theme_bw() +
  theme(legend.position=c(0.8, 0.7),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(),
        text = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  ylim(0,1)


ggsave(paste0(drive, "Plots2\\AJE-01109-2021 Bermingham Figure 1.pdf"), height=5, width=7, dpi=300)
ggsave(paste0(drive, "Plots2\\AJE-01109-2021 Bermingham Figure 1.eps"), height=5, width=7, dpi=300)
ggsave(paste0(drive, "Plots2\\AJE-01109-2021 Bermingham Figure 1.png"), height=5, width=7, dpi=300)


#-----------------------------------#
# Fig 2 - Regression discontinuity design plots 
# for treatment and outcome variables by age in months 
#-----------------------------------#

rdd  <- read_csv("data/rdd_all_39_follow_up_15.csv") 

rdd_covid=rdrobust::rdrobust(y=rdd$log_odds_covid , x=rdd$age, c=960,
                   fuzzy = rdd$p_treat, p=1, h=60, kernel="uniform", weights = rdd$num_people)

disc <- read_csv('data/rdd_discontinuity_values_39_follow_up_15_vacc_type.csv')



fontsize <- 12

# formated discontinuity estimates
d_anti = paste0(formatC(disc$antibody_me, 3), " [", formatC(disc$antibody_LCL	, 3), 
                " to ",formatC(disc$antibody_UCL	, 3), "]" )

d_covid = paste0(formatC(disc$covid_me*100000, 1,format = "f"), " [", formatC(disc$covid_LCL*100000	, 1,format = "f"), 
                 " to ",formatC(disc$covid_UCL*100000	,1, format = "f"), "]" )
d_other = paste0(formatC(disc$other_me*100000, 1,format = "f"), " [", formatC(disc$other_LCL*100000, 1,format = "f"),
                 " to ",formatC(disc$other_UCL	*100000, 1,format = "f"), "]" )
## graph - antibodies
disc_treated <- ggplot(data= rdd, aes(x = age/12, y = p_treat, group=dummy, weight=num_people))+
  geom_point(colour="black", alpha=0.3)+
  geom_smooth(method = lm, colour='black')+
  geom_vline(xintercept=80, linetype="dotted")+
  labs(title = "A)",
       x="Age, Years", y="Weighted Proportion \nProtected")+
  theme_bw() +
  ylim(0, 0.3)+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = fontsize, colour = "black"),
        axis.text = element_text(size = fontsize, colour = "black"),
        plot.title = element_text(size = fontsize, hjust = -0.11))
disc_treated
disc_treated + ggsave(paste0("Plots2/AJE-01109-2021 Bermingham Figure 2A.eps"), height=5, width=7, dpi=300)

disc_covid <- ggplot(data= rdd, aes(x = age/12, y = p_covid*100000, group=dummy, weight=num_people))+
  geom_point(colour="black", alpha=0.3)+
  geom_smooth(method = lm, colour='black')+
  geom_vline(xintercept=80, linetype="dotted")+
  labs(title = "B)",
       x="Age, Years", y="15-day Mortality Rate, \n100,000 People")+
  theme_bw() +
  ylim(0, 400) +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = fontsize, colour = "black"),
        axis.text = element_text(size = fontsize, colour = "black"),
        plot.title = element_text(size = fontsize, hjust = -0.12))
disc_covid
disc_covid + ggsave(paste0("Plots2/AJE-01109-2021 Bermingham Figure 2B.eps"), height=5, width=7, dpi=300)


disc_non_covid <- ggplot(data= rdd, aes(x = age/12, y = p_other*100000, group=dummy, weight=num_people))+
  geom_point(colour="black", alpha=0.3)+
  geom_smooth(method = lm, colour='black')+
  geom_vline(xintercept=80, linetype="dotted")+
  labs(title = "C)",
       x="Age in Years", y="15-day Mortality Rate, \n100,000 People")+
  theme_bw() +
  ylim(0, 400) +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = fontsize, colour = "black"),
        axis.text = element_text(size = fontsize, colour = "black"),
        plot.title = element_text(size = fontsize, hjust = -0.12))
disc_non_covid 
disc_non_covid + ggsave(paste0("Plots2/AJE-01109-2021 Bermingham Figure 2C.eps"), height=5, width=7, dpi=300)

# option 1
ggpubr::ggarrange(disc_treated, disc_covid, disc_non_covid, nrow=3)+
  ggsave("Plots2/AJE-01109-2021 Bermingham Figure 2.png",width = 7, height = 9) +
  ggsave("Plots2/AJE-01109-2021 Bermingham Figure 2.pdf",width = 7, height = 9, dpi=300)
# option 2
#ggpubr::ggarrange(treated, covid,  nrow=1)+
#  ggsave("Plots2/rdd_graph_wide.png",width = 10, height = 6)


####-------------------------------------------####
#### graphs for discontinuities by time period ####
####-------------------------------------------####

save_location <- "final_plots/"


main_est_colour= "#a9a9a9"
fontsize <- 10
hjust <- -0.25

# main dataset
data <- read.csv("data/results_all_15_all_20221111.csv")

# combine 
data <- dplyr::filter(data, day.number >= 30, day.number <=80) %>%
  mutate(date = as.Date('2020-12-08')+day.number,
         VE = VE*100,
         VE_UCL = VE_UCL*100,
         VE_LCL = VE_LCL*100)



# Fig 2 - A Discontinuity in the weighted proportion of people who have been vaccinated and reached threshold antibody level over the 15-day period




treat_discontinuity <- ggplot(data=filter(data, day.number!=39), aes(x=date, y=treat_from_treat_rdd))+
  geom_point(colour=main_est_colour, size=1)+
  geom_errorbar(data=filter(data, day.number!=39), aes(ymin=treat_LCL, ymax=treat_UCL), colour=main_est_colour)+
  geom_point(data=filter(data, day.number==39), aes(x=date, y=treat_from_treat_rdd), shape=17, size=1)+
  geom_errorbar(data=filter(data, day.number==39), aes(ymin=treat_LCL, ymax=treat_UCL)) +
  labs(title = "A)",
       y = "Difference in Protection",
       x = "Index Date") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(),
        text = element_text(size = fontsize, colour = "black"),
        axis.text = element_text(size = fontsize, colour = "black"),
        legend.text = element_text(size = fontsize),
        legend.title = element_text(size = fontsize),
        plot.title = element_text(size = fontsize, hjust = -0.15)) +
  ylim(0, 0.4) +
  scale_x_date(date_labels = '%b %d', breaks = seq(min(data$date), as.Date("2021-03-12"), by="17 days"), limits = as.Date(c('2021-01-06', '2021-03-01')))
treat_discontinuity

treat_discontinuity + ggsave(paste0(save_location, "AJE-01109-2021 Bermingham Figure 3A.eps"), width=7/2, height=9/3, dpi=300)
treat_discontinuity + ggsave(paste0(save_location, "AJE-01109-2021 Bermingham Figure 3A.pdf"), width=7/2, height=9/3, dpi=300)     

# Fig 2 - B discontinuity in the 15 day COVID-19 mortality rate 

COVID_discontinuity <- ggplot(data=filter(data, day.number!=39), aes(x=date, y=covid_discontinuity * 100000, ymin=covid_LCL* 100000, ymax=covid_UCL* 100000))+
  geom_point(colour=main_est_colour, size=1)+
  geom_errorbar(colour=main_est_colour) +
  geom_point(data=filter(data, day.number==39), shape=17, size=1)+
  geom_errorbar(data=filter(data, day.number==39)) +
  labs(title = "B)",
       y = "Effect on 15-day COVID-19 \nMortality Rate per 100,000 People",
       x = "Index Date") +
  scale_x_date(date_labels = '%b %d', breaks = seq(min(data$date), as.Date("2021-03-12"), by="17 days"), limits = as.Date(c('2021-01-06', '2021-03-01')))+
  theme_bw()+
  ylim(-46, 46)+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(),
        text = element_text(size = fontsize, colour = "black"),
        axis.text = element_text(size = fontsize, colour = "black"),
        legend.text = element_text(size = fontsize),
        legend.title = element_text(size = fontsize),
        plot.title = element_text(size = fontsize, hjust = hjust))

COVID_discontinuity
ggsave(paste0(save_location, "AJE-01109-2021 Bermingham Figure 3B.eps"), width=7/2, height=9/3, dpi=300)
ggsave(paste0(save_location, "AJE-01109-2021 Bermingham Figure 3B.pdf"), width=7/2, height=9/3, dpi=300)



# Fig 2 C) discontinuity in the 15 day non-COVID-19 mortality rate

non_COVID_discontinuity <- ggplot(data=filter(data, day.number!=39), aes(x=date, y=other_discontinuity* 100000, ymin=other_LCL* 100000, ymax=other_UCL* 100000))+
  geom_point(colour=main_est_colour, size=1)+
  geom_errorbar(colour=main_est_colour) +
  geom_point(data=filter(data, day.number==39), shape=17, size=1)+
  geom_errorbar(data=filter(data, day.number==39)) +
  labs(title = "C)",
       y = "Effect on 15-day Non-COVID-19 \nMortality Rate per 100,000 People",
       x = "Index Date") +
  scale_x_date(date_labels = '%b %d', breaks = seq(min(data$date), as.Date("2021-03-12"), by="17 days"), limits = as.Date(c('2021-01-06', '2021-03-01')))+
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(),
        text = element_text(size = fontsize, colour = "black"),
        axis.text = element_text(size = fontsize, colour = "black"),
        legend.text = element_text(size = fontsize),
        legend.title = element_text(size = fontsize),
        plot.title = element_text(size = fontsize, hjust = -0.25)) +
  ylim(-41, 40)
non_COVID_discontinuity         
ggsave(paste0(save_location, "AJE-01109-2021 Bermingham Figure 3C.eps"), width=7/2, height=9/3, dpi=300)
ggsave(paste0(save_location, "AJE-01109-2021 Bermingham Figure 3C.pdf"), width=7/2, height=9/3, dpi=300)

# Fig 2 ; D) effect of vaccination on #COVID-19 mortality in people aged 80 (= (B)/(A))

LATE <- ggplot(data=filter(data, day.number!=39), aes(x = date, y = LATE_covid * 100000, 
                                                      ymin = LATE_covid_LCL * 100000, 
                                                      ymax = LATE_covid_UCL * 100000))+
  geom_point(colour=main_est_colour, size=1)+
  geom_errorbar(colour=main_est_colour) +
  geom_point(data=filter(data, day.number==39), shape=17, size=1)+
  geom_errorbar(data=filter(data, day.number==39)) +
  labs(title = "D)",
       y = "Effect on 15-day COVID-19 \nMortality Rate  per 100,000 People",
       x = "Index Date") +
  scale_x_date(date_labels = '%b %d', breaks = seq(min(data$date), as.Date("2021-03-12"), by="17 days"), limits = as.Date(c('2021-01-06', '2021-03-01')))+
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(),
        text = element_text(size = fontsize, colour = "black"),
        axis.text = element_text(size = fontsize, colour = "black"),
        legend.text = element_text(size = fontsize),
        legend.title = element_text(size = fontsize),
        plot.title = element_text(size = fontsize, hjust = -0.28)) +
  ylim(-400, 200)
LATE
LATE+ ggsave(paste0(save_location, "AJE-01109-2021 Bermingham Figure 3D.eps"), width=7/2, height=9/3, dpi=300)
LATE+ ggsave(paste0(save_location, "AJE-01109-2021 Bermingham Figure 3D.pdf"), width=7/2, height=9/3, dpi=300)



# Fig 2 E) vaccine effectiveness (= 1 - odds ratio for being protected).


VE <- ggplot(data=filter(data, day.number!=39), aes(x=date, y= VE,
                                                               ymin = VE_LCL,
                                                               ymax = VE_UCL)) +
  geom_point(colour=main_est_colour, size=1)+
  geom_errorbar(colour=main_est_colour) +
  geom_point(data=filter(data, day.number==39), shape=17, size=1)+
  geom_errorbar(data=filter(data, day.number==39)) +
  labs(title = "E)",
       y = "Relative Reduction in Risk of \nCOVID-19 Death, %",
       x = "Index Date") +
  scale_x_date(date_labels = '%b %d', breaks = seq(min(data$date), as.Date("2021-03-12"), by="17 days"), limits = as.Date(c('2021-01-06', '2021-03-01')))+
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(),
        text = element_text(size = fontsize, colour = "black"),
        axis.text = element_text(size = fontsize, colour = "black"),
        legend.text = element_text(size = fontsize),
        legend.title = element_text(size = fontsize),
        plot.title = element_text(size = fontsize, hjust = -0.28)) +
  ylim(-200, 100)
VE
ggsave(paste0(save_location, "AJE-01109-2021 Bermingham Figure 3E.eps"), width=7/2, height=9/3, dpi=300)
ggsave(paste0(save_location, "AJE-01109-2021 Bermingham Figure 3E.pdf"), width=7/2, height=9/3, dpi=300)


# for paper
ggpubr::ggarrange(treat_discontinuity, COVID_discontinuity,non_COVID_discontinuity, LATE , VE, nrow=3, ncol=2 )+
  ggsave(paste0(save_location, "AJE-01109-2021 Bermingham Figure 3.pdf"), width=7, height=9, dpi=300)
