

#---------------------------------------------------#
# Code to run the RDD
# input data:
# analytical_linked_vaccination_dataset_20210524: 
# linked Census to NIMS, GDPPR and mortality
#
# dependencies: functions.R
# Outputs: 
# * rdd_all_39_follow_up_15.csv: Protection and mortality rates by age in month;
#   used for the estimation of the main RDD results and plots
# * results_all_15_all.csv: RDD results for 15 day follow-up, 
#   with different starting time for the follow-up period; 
# * cut_offs_39_follow_up_15.csv: placebo tests using different cut-offs for treatment
#   (testing that there is no discontinuity in the outcome or treatment at other ages)
# * continuity_data.csv: check for continuity of covariates at the cut-off
# * bandwidths_39_follow_up_15.csv: RDD estimates with different bandwidths
# * doughnut_rdd_39_follow_up_15.csv: RDD estimates , removing people very close
#   to the cut-off
#---------------------------------------------------#




library(plyr)
library(dplyr)
library(sparklyr)
library(survival)
library(ggplot2)
library(tidyverse)
library(ggfortify)
library(lubridate)
library(broom)
library(rdrobust)


#----------------------------------------#
# Set up the spark connection
#----------------------------------------#


config <- spark_config() 
config$spark.dynamicAllocation.maxExecutors <- 30
config$spark.executor.cores <- 5
config$spark.executor.memory <- "20g"
config$spark.driver.maxResultSize <- "10g"
sc <- spark_connect(master = "yarn-client",
                    app_name = "R_Example",
                    config = config,
                    version = "2.3.0")

#----------------------------------------#
# Set file locations and load functions
#----------------------------------------#


# location to store the results
dir ="cen_dth_gps/Vaccination analysis/Vaccine effectiveness/RDD"
out_dir= paste0(dir, "/Results")
out_dir_paper_results= paste0(dir, "/Results/paper_update_results")

# source functions
fun_dir = paste0(dir, "/_Functions")
fun_list <- paste(fun_dir, list.files(fun_dir), sep="/")
for(i in 1:length(fun_list)) {source(fun_list[i])}

#--------------------------
# Create sparklyr dataframe
#--------------------------

# first day of vaccinations
day_zero <- "2020-12-08"

# Define spark dataframe
df <- sdf_sql(sc, "SELECT * FROM cen_dth_gps.analytical_linked_vaccination_dataset_20210524_fixed") 

# determine end of study date (last dod minus 7 days)

min_max <- df %>%
dplyr::summarise(min(first_vacc_date), max(first_vacc_date), min(antibody_date),  max(antibody_date), min(dod), max_dod = max(dod))%>%
collect()
print(min_max)

eos_date = as.character(min_max$max_dod -7) 
print(eos_date)

## count N people vaccinated before eos_date, not linked to GPES/Census
unlinked <- df  %>%
  filter(is.na(b_asthma) & (!is.na(first_vacc_date) & first_vacc_date<=eos_date ))  %>%
  filter(dod>=day_zero | is.na(dod), age_when_eligible_months >= 75*12 & age_when_eligible_months <85*12)

n_unlinked = count(unlinked)%>%
             collect()

# filter to only include those in gpes (hence alive in 2019), alive on 7th Dec
# resident in england and in census and gpes

df <- df %>%
  filter(present_in_gpes == 1 & !is.na(b_asthma) & present_in_census == 1)  %>%
  # remove people who died before the first day of the vaccination data
  filter(dod>=day_zero | is.na(dod)) %>%
  # Remove those not living in England
  mutate(country = substr(rgn_derived, 1, 1)) %>%
  filter(country == "E")


## percent of populaiton we cover
test <- df %>%
  filter(is.na(dod) | dod>"2019-06-01") %>%
  filter(datediff("2019-06-01", dob)>=365.25*18) %>%
  count() %>%
  collect()

test/44263393

# some vaccine types are blank for vaccinated people - record as 'unknown' instead
df <- df %>%
  mutate(first_vacc_manufacturer = ifelse((vaccinated == 1) & is.na(first_vacc_manufacturer), 'unknown', first_vacc_manufacturer))

df <- df %>%
# remove Moderna vaccines and change unknown to be Pfizer
  filter(is.na(first_vacc_manufacturer) | first_vacc_manufacturer!= 'Moderna' & first_vacc_manufacturer!= 'Moderna, Inc') %>%
  mutate(first_vacc_manufacturer = ifelse(first_vacc_manufacturer=='unknown', 'Pfizer', first_vacc_manufacturer))

# filter to include only people age 75-84 when eligible and create age group groups
# create death variable for whether a death has occurred
df <- df %>% 
  filter(age_when_eligible_months >= 75*12 & age_when_eligible_months <85*12) %>%
  dplyr::mutate(age_group_when_eligible = ifelse(age_when_eligible_years >= 80, "80-84", "75-79"),
               death = ifelse(!is.na(covid_death_status),1,0))

# if was vaccinated on or after the eos_date, change the vaccination
# date to NULL and change the vaccinated flag to 0 
# if someone died after the eos_date, change the covid_death_status and dod to NULL
# if the antibody date is after the eos_date, change it to NULL
df <- df %>%
  dplyr::mutate(first_vacc_date = ifelse(first_vacc_date >=eos_date, NULL, first_vacc_date),
        vaccinated = ifelse(first_vacc_date >=eos_date|is.na(first_vacc_date), 0, 1),
        antibody_date = ifelse(antibody_date >=eos_date, NULL, antibody_date))

##  # count observations 
n_all <- df %>% 
         group_by(vaccinated)%>%
         count()%>%
         collect()%>%
         spread(vaccinated, n)

# remove time from date-time entries in vaccination dates
# Round age in months to nearest number of whole months
df <- df%>% 
  dplyr::mutate(first_vacc_date = date(timestamp(first_vacc_date)),
                second_vacc_date = date(timestamp(second_vacc_date)),
                age_when_eligible_months = floor(age_when_eligible_months)) 

# filter to include only those who are not vaccinated or were vaccinated
# after 8.12.2020 (this removes the few eroneous vaccination entries
# before this date)
df <- df %>%
  filter(first_vacc_date >= "2020-12-08" | is.na(first_vacc_date)) 

colnames(df)

#----------------------------------------#
#        Sample flow
#----------------------------------------#

# count observations 
n_2 <- df %>%
       group_by(vaccinated)%>%
         count()%>%
         collect()%>%
         spread(vaccinated, n)

sample_flow <- rbind(n_all, n_2)%>%
                rename(unvaccinated = `0`,
                      vaccinated = `1`)%>%
               mutate(n= vaccinated+unvaccinated)
sample_flow                      

write.csv(sample_flow, paste0(out_dir , "/sample_flow.csv"))


#----------------------------------------#
# Table 1 - summary statistics
#----------------------------------------#

## Identifying deaths involving Covid-19

df <- df %>%
  mutate(u071 = ifelse(fic10und == "U071" | fic10men1 == "U071" | fic10men2 == "U071" |
                       fic10men3 == "U071" | fic10men4 == "U071" | fic10men5 == "U071" |
                      fic10men6 == "U071" | fic10men7 == "U071" | fic10men8 == "U071" |
                      fic10men9 == "U071" | fic10men10 == "U071" | fic10men11 == "U071" |
                      fic10men12 == "U071" | fic10men13 == "U071" | fic10men14 == "U071" |
                      fic10men15 == "U071", 1, 0),
        due_to = ifelse(fic10und == "U071" | fic10und == "U072", 1, 0))

df %>% group_by(covid_death_status, u071) %>% count() %>% collect()
df %>% group_by(covid_death_status, due_to) %>% count() %>% collect()

## recode vaccine manufacturer 
df <- df %>%
      mutate(manufacturer = case_when( vaccinated == 0 ~ "Not vaccinated",
                                      vaccinated == 1 ~ first_vacc_manufacturer))

## variables of interest 
columns <- c('covid_death_status','vaccinated', 'sex', 'age_when_eligible_years',
             'clin_ext_vuln_flag', 'ethnicity','manufacturer')

continuous_vars <- c('age_when_eligible_years')
binary_vars <- c('vaccinated')

sum_stats <- get_tables_per_var(df, columns)
sum_stats

saveRDS(sum_stats, paste0(out_dir, '/summary_stats_list.rds'))

## stratified by age 

sum_stats_80 <- get_tables_per_var(filter(df,age_when_eligible_years >= 80), columns)
sum_stats_75 <- get_tables_per_var(filter(df,age_when_eligible_years < 80), columns)

saveRDS(sum_stats_80, paste0(out_dir, '/summary_stats_list_80.rds'))
saveRDS(sum_stats_75, paste0(out_dir, '/summary_stats_list_75.rds'))


## format tables ##
var_lookup <- read.csv(paste0('cen_dth_hes/COVID19_occupations', '/varnames_lookup.csv'))%>%
select(variable=names, labels, type)

   
sum_stats <- readRDS(paste0(out_dir,'/summary_stats_list.rds'))
sum_stats_80 <- readRDS(paste0(out_dir,'/summary_stats_list_80.rds'))
sum_stats_75 <- readRDS(paste0(out_dir,'/summary_stats_list_75.rds'))


all <- format_table(sum_stats, name = "all" )                  
a75 <- format_table(sum_stats_75, name= "75-79")                  
a80 <- format_table(sum_stats_80, name= "80-84")                  


table1 <- select(all, variable, value, all ) %>%
left_join(select(a75, variable, value, `75-79` ) )    %>%
left_join(select(a80,  variable, value , `80-84`))   

## Save final table ##
write.csv(table1, paste0(out_dir, '/table1.csv'))




# -----------------------------------------------
# Proportion of vaccinated people linked to GPES
# -----------------------------------------------

linkage <- sdf_sql(sc, "SELECT * FROM cen_dth_gps.analytical_linked_vaccination_dataset_20210524") %>%
 filter(age_when_eligible_months >= 75*12 & age_when_eligible_months <85*12,
          (!is.na(first_vacc_date) & first_vacc_date<=eos_date & first_vacc_date >= day_zero),
        dod>=day_zero | is.na(dod))%>%
 mutate(linked = ifelse(!is.na(b_asthma),1,0))%>%
 group_by(linked)%>%
count()%>%
collect()

write.csv(linkage, paste0(out_dir, '/linkage.csv'))

#-------------------------------------------------
# Plot vaccinations by age group 
# accounting for competing risk (death)
# (Web Figure 3) 
#-------------------------------------------------


# dataframe to use for survival analysis with censure (event=0) for death or end 
# of study and event=1 for vaccinated
df_vaccination <- df %>%
       dplyr::mutate(event = case_when(!is.na(first_vacc_date) & (first_vacc_date<dod | is.na(dod))  ~ 1,
                                   death==1 ~ 0)) %>%
       dplyr::mutate(t = case_when(event == 1 ~ datediff(first_vacc_date, day_zero),
                     event == 0 ~ datediff(dod, day_zero),
                     is.na(event) ~ datediff(eos_date, day_zero)))

# add a 0 instead of NA for people who are alive and unvaccinated
# at the end of the time period so that they are censured
df_vaccination <- df_vaccination %>%  
  na.replace(event = 0)

# create R dataframe for survival analysis
data_df_vaccination <- df_vaccination %>%
           select(t, event, age_group_when_eligible, age_when_eligible_months, age_when_eligible_years)%>%
           collect()

# fit survival curve for vaccination probability by age group
mfit <- survfit(Surv(t, as.factor(event)) ~ age_group_when_eligible, 
                 data = data_df_vaccination)


# extract results
fit_df_age_group <- data.frame(t = summary(mfit)$time,
                     n.risk= summary(mfit)$n.risk[,1],
                    prop_vaccinated = summary(mfit)$pstate[,2],
                    age_group = summary(mfit)$strata)
fit_df_age_group$age_group <- revalue(fit_df_age_group$age_group, 
                                      c("age_group_when_eligible=75-79"="75-79", 
                                        "age_group_when_eligible=80-84"="80-84"))

# draw graph 
ggplot(data = fit_df_age_group, aes(x=t, prop_vaccinated, colour = age_group) )+
  geom_step() +
  labs(title= "Cumulative incidence of vaccination by age group",
       x="Days after Dec 08 2020", y="Proportion vaccinated") +
  scale_colour_discrete(name = "",  
                      labels = c("75-79", "80-84"))+
  theme_bw() +
  theme(legend.position = "top") +
  ylim(0, 1)

# save results
ggsave(paste0(out_dir_paper_results, "/plot_data/cum_inc_vaccination_age_group.png"), height=6, width=4/3*6)
write.csv(fit_df_age_group, paste0(out_dir_paper_results, "/plot_data/cum_inc_vaccination_age_group.csv"))



#----------------------------------------------------------#
# Run model for a particular start date and follow up time #
#----------------------------------------------------------#

t_start <- 39
follow_up <- 15

#--------------------------
# CIS antibody fits
#--------------------------

prob_antibody_pos <- read.csv(file = paste0(dir, '/logistic_80y_Figure1.csv'))

IgG_levels <- read.csv(file = paste0(dir, '/linear_80y_Figure3.csv'))


#--------------------------
# Main results - Run model 
#--------------------------

## run once
results <- full_rd_model(t_start = 39, follow_up = 14, df=df)
results

data.frame(as.list(results))
write.csv(data.frame(as.list(results)),
          paste0(out_dir_paper_results, "/results_39_followup_14.csv"), row.names = FALSE)

out_dir_paper_results


#--------------------------
# Different start times    #
#--------------------------
follow_up <- 15

## loop through different start times

start_t <- 20
end_t <- 80
t_start_values <- start_t:end_t

values <- sapply(t_start_values, full_rd_model, follow_up=follow_up, df = df)


results_df <- data.frame("day number"=rbind(values)[1,], 
                      "follow up time" = rbind(values)[2,], 
                       "treatment_discontinuity"=rbind(values)[3,], 
                        "treat_upper"=rbind(values)[4,],
                        "treat_lower"=rbind(values)[5,],
                        "covid_discontinuity"=rbind(values)[6,], 
                       "other_discontinuity"=rbind(values)[7,], 
                        "LATE_covid"=rbind(values)[8,], 
                        "LATE_covid_LCL"=rbind(values)[9,], 
                        "LATE_covid_UCL"=rbind(values)[10,], 
                        "LATE_covid_p_value"=rbind(values)[11,],
                        "LATE_other"=rbind(values)[12,], 
                        "LATE_other_LCL"=rbind(values)[13,], 
                        "LATE_other_UCL"=rbind(values)[14,], 
                        "LATE_other_p_value"=rbind(values)[15,],
                      "treat_from_treat_rdd"=rbind(values)[16,],
                      "treat_LCL"=rbind(values)[17,],
                      "treat_UCL"=rbind(values)[18,],
                      "treat_se"=rbind(values)[19,],
                        "LATE_covid_se"=rbind(values)[20,],
                         "VE" = rbind(values)[21,],
                      "VE_LCL" = rbind(values)[22,],
                      "VE_UCL" = rbind(values)[23,])

print(results_df)
write.csv(results_df,
          paste0(out_dir_paper_results, "/results_all_",
            toString(follow_up),"_days_", toString(start_t),
                "_", toString(end_t), ".csv"), row.names = FALSE)

results_df

# VN: not sure I understand the file_numbers
# combine files when running separately
file_numbers <- list('20_30', '31_40', '41_50', '51_60', '61_70','71_80')

all_t_data<- read.csv(file = paste0(out_dir_paper_results, 
                                       '/results_all_15_days_', 
                                    file_numbers[1], '.csv'))

for (numbers in file_numbers[-1]) {
  t_segment <- read.csv(file = paste0(out_dir_paper_results, 
                                      '/results_all_15_days_', numbers, '.csv'))
  all_t_data <- rbind(all_t_data, t_segment)
}

write.csv(all_t_data,
          paste0(out_dir_paper_results, "/results_all_15_all.csv"), row.names = FALSE)


#####################################
# Different analysis period lengths #
#####################################

values_analysis_period <- sapply(seq(5, 25, 5), full_rd_model, t_start = 39, df = df)
values <- values_analysis_period

results_df <- data.frame("day number"=rbind(values)[1,], 
                      "follow up time" = rbind(values)[2,], 
                       "treatment_discontinuity"=rbind(values)[3,], 
                        "treat_upper"=rbind(values)[4,],
                        "treat_lower"=rbind(values)[5,],
                        "covid_discontinuity"=rbind(values)[6,], 
                       "other_discontinuity"=rbind(values)[7,], 
                        "LATE_covid"=rbind(values)[8,], 
                        "LATE_covid_LCL"=rbind(values)[9,], 
                        "LATE_covid_UCL"=rbind(values)[10,], 
                        "LATE_covid_p_value"=rbind(values)[11,],
                        "LATE_other"=rbind(values)[12,], 
                        "LATE_other_LCL"=rbind(values)[13,], 
                        "LATE_other_UCL"=rbind(values)[14,], 
                        "LATE_other_p_value"=rbind(values)[15,],
                      "treat_from_treat_rdd"=rbind(values)[16,],
                      "treat_LCL"=rbind(values)[17,],
                      "treat_UCL"=rbind(values)[18,],
                      "treat_se"=rbind(values)[19,],
                        "LATE_covid_se"=rbind(values)[20,],
                         "VE" = rbind(values)[21,],
                      "VE_LCL" = rbind(values)[22,],
                      "VE_UCL" = rbind(values)[23,])

print(results_df)
write.csv(results_df,
          paste0(out_dir_paper_results, "/results_different_follow_up_t_start_",
                  t_start, ".csv"), row.names = FALSE)





#######################
# Create plots #
#######################
t_start <- 78
follow_up <- 15

plot_data <- read.csv(file = paste0(out_dir_paper_results, 
                              "/plot_data/plot_data_t_start_in_gpes_", toString(t_start),
                 "_follow_up_", toString(follow_up), ".csv"))

plot_data



#--------------------------
# Fit and plot rdd for percentage treated by age in months
#--------------------------


antibody_rdd <- fit_rdd(plot_data, "age", "p_treat", 960)
antibody_rdd$main_effect
antibody_rdd$main_effect_error
antibody_rdd$p_value
antibody_rdd$LCL
antibody_rdd$UCL

cut_off_date = as.character(as.Date("20201208", format = "%Y%m%d") +t_start) 
title = paste0("Weighted proportion treated over ", toString(follow_up), " days from ", toString(cut_off_date))
x_lab = "Age on Dec 8 2020 (years)"
y_lab = "Weighted proportion treated"
save_name = paste0(out_dir_paper_results,"/plot_data/", "weighted_prop_treat_", 
                   toString(follow_up), "_days_from_",toString(t_start), 
                   "_months_vacc_type_new_follow_up.png")

antibody_plot <- plot_rdd(antibody_rdd$rdd_dataframe, "age", "p_treat", 960, antibody_rdd$fit, "D",
        TRUE, title, x_lab, y_lab, save_name)
antibody_plot

write.csv(antibody_rdd$rdd_dataframe, paste0(out_dir_paper_results,
                                             "/plot_data/rdd_treatment_", 
                                             toString(t_start), 
                                             "_follow_up_", 
                                             toString(follow_up), 
                                             "_vacc_type_new_follow_up.csv"))

antibody_data <- plot_rdd_data(antibody_rdd$rdd_dataframe, "age", "p_treat", 960, antibody_rdd$fit, "D")
antibody_data
#### Fit and plot rdd for Other deaths

mortality_type <- "p_other"
mortality_rdd_other <- fit_rdd(plot_data, "age", mortality_type, 960)
mortality_rdd_other                 

title = paste0("Probability of non Covid-19 mortality by age over ", toString(follow_up), " days after ", toString(cut_off_date))
x_lab = "Age on Dec 8 2020 (years)"
y_lab = paste0("Non Covid-19 mortality probability within ", toString(follow_up), " days")
save_name = paste0(out_dir_paper_results,"/plot_data/other-mortality_", 
                   toString(follow_up), "days_from_", toString(t_start), 
                   "_months_vacc_type_new_follow_up.png")

other_plot <- plot_rdd(mortality_rdd_other$rdd_dataframe, "age", mortality_type, 960, mortality_rdd_other$fit, "D",
        TRUE, title, x_lab, y_lab, save_name)
other_plot 

write.csv(mortality_rdd_other$rdd_dataframe, paste0(out_dir_paper_results,
                                                    "/plot_data/rdd_other_", 
                                                    toString(t_start), 
                                                    "_follow_up_", 
                                                    toString(follow_up), 
                                                    "_vacc_type_new_follow_up.csv"))

other_data <- plot_rdd_data(mortality_rdd_other$rdd_dataframe, "age", "p_treat", 960, mortality_rdd_other$fit, "D")
other_data

#### Fit and plot rdd for Covid deaths

mortality_type <- "p_covid"

## calculate log odds ratios
#  plot_data <- plot_data %>% 
#    mutate(log_odds_covid = log( p_covid / (1-p_covid)),
#         log_odds_treat = log( p_treat / (1-p_treat)),
#         log_odds_other = log( p_other / (1-p_other)))
#
#plot_data[plot_data$p_covid != 0,]

mortality_rdd_covid <- fit_rdd(plot_data, "age", mortality_type, 960)
mortality_rdd_covid        

title = paste0("Probability of COVID-19 mortality by age over ", toString(follow_up), " days after ", toString(cut_off_date))
x_lab = "Age on Dec 8 2020 (years)"
y_lab = paste0("Covid mortality probability within ", toString(follow_up), " days")
save_name = paste0(out_dir_paper_results,"/plot_data/covid-mortality_", 
                   toString(follow_up), "days_from_", toString(t_start), 
                   "_months_vacc_type_new_follow_up.png")

covid_plot <- plot_rdd(mortality_rdd_covid$rdd_dataframe, "age", mortality_type, 960, mortality_rdd_covid$fit, "D",
        TRUE, title, x_lab, y_lab, save_name)
covid_plot

write.csv(mortality_rdd_covid$rdd_dataframe, paste0(out_dir_paper_results,
                                                    "/plot_data/rdd_covid_", 
                                                    toString(t_start), 
                                                    "_follow_up_", 
                                                    toString(follow_up), 
                                                    "_vacc_type_new_follow_up.csv"))

covid_data <- plot_rdd_data(mortality_rdd_covid$rdd_dataframe, "age", "p_treat", 960, mortality_rdd_covid$fit, "D")
covid_data

# save values for text in figure caption in paper
discontinuity_values <- data.frame(antibody_me = antibody_rdd$main_effect,
                                   antibody_err = antibody_rdd$main_effect_error,
                                  antibody_p = antibody_rdd$p_value,
                                   antibody_LCL = antibody_rdd$LCL,
                                   antibody_UCL = antibody_rdd$UCL,
                                  covid_me = mortality_rdd_covid$main_effect,
                                   covid_err = mortality_rdd_covid$main_effect_error,
                                  covid_p = mortality_rdd_covid$p_value,
                                    covid_LCL = mortality_rdd_covid$LCL,
                                   covid_UCL = mortality_rdd_covid$UCL,
                                  other_me = mortality_rdd_other$main_effect,
                                   other_err = mortality_rdd_other$main_effect_error,
                                  other_p = mortality_rdd_other$p_value,
                                  other_LCL = mortality_rdd_other$LCL,
                                   other_UCL = mortality_rdd_other$UCL)

write.csv(discontinuity_values, paste0(out_dir_paper_results,
                                       "/plot_data/rdd_discontinuity_values_", 
                                       toString(t_start), "_follow_up_", 
                                       toString(follow_up), "_vacc_type_new_follow_up.csv"),
                                      row.names=FALSE)


# combine data for plotting outside of TRE

all_data <- merge(plot_data, antibody_data, on="age", all=TRUE)
all_data <- all_data %>% 
  rename(prediction_antibody = prediction,
        upper_antibody=upper,
        lower_antibody=lower)

all_data <- merge(all_data, covid_data, on="age", all=TRUE)
all_data <- all_data %>% 
  rename(prediction_covid = prediction,
        upper_covid=upper,
        lower_covid=lower)

all_data <- merge(all_data, other_data, on="age", all=TRUE)
all_data <- all_data %>% 
  rename(prediction_other = prediction,
        upper_other=upper,
        lower_other=lower)


write.csv(all_data, 
          paste0(out_dir_paper_results,"/plot_data/rdd_all_", toString(t_start), 
                 "_follow_up_", toString(follow_up), "_new_follow_up.csv"),
         row.names = FALSE)

#----------------------------------
#----------------------------------
# Sensitivity analyses 
#----------------------------------
#----------------------------------


#----------------------------------
#Continuity in covariates 
#----------------------------------

total <- df %>%
  group_by(age_when_eligible_months) %>%
  count() %>%
  rename(total=n)%>%
  collect()

female <- df %>%
  filter(sex==2) %>%
  group_by(age_when_eligible_months) %>%
  count() %>%
  rename(female_n=n) %>%
  collect()

imd <- df %>%
  filter(imd_quintile==1) %>%
  group_by(age_when_eligible_months) %>%
  count() %>%
  rename(imd_n=n) %>%
  collect()

clin_vuln <- df %>%
  filter(clin_ext_vuln_flag==1) %>%
  group_by(age_when_eligible_months) %>%
  count() %>%
  rename(vuln_n=n) %>%
  collect()

data_all <- merge(total, female, all=TRUE)
data_all <- merge(data_all, imd, all=TRUE)
data_all <- merge(data_all, clin_vuln, all=TRUE)
data_all <- data_all %>%
  mutate(perc_female = female_n/total,
        perc_imd = imd_n/total,
        perc_clin_vuln = vuln_n/total,
        D = ifelse(age_when_eligible_months >= 960, 1, 0))


write.csv(data_all, 
          paste0(out_dir_paper_results,"/sensitivity_analysis/continuity_data.csv"),
         row.names = FALSE)

#----------------------------------------------
# Discontinuity at other places and bandwidth #
#----------------------------------------------


#read data in

t_start <- 78
follow_up <- 15

plot_data <- read.csv(paste0(out_dir_paper_results, "/plot_data/plot_data_t_start_in_gpes_", toString(t_start),
                 "_follow_up_", toString(follow_up), ".csv"))



rdd_with_params <- function(cut_off, bandwidth, dataframe){
  rd_covid <- rdrobust(y=dataframe$p_covid, x=dataframe$age, c=cut_off, 
           fuzzy=dataframe$p_treat, p=1, h=bandwidth, kernel='uniform', 
                      weights = dataframe$num_people)

  rd_other <- rdrobust(y=dataframe$p_other, x=dataframe$age, c=cut_off, 
           fuzzy=dataframe$p_treat, p=1, h=bandwidth, kernel='uniform', 
                      weights = dataframe$num_people)

  rd_treat <- rdrobust(y=dataframe$p_treat, x=dataframe$age, c=cut_off, 
           p=1, h=bandwidth, kernel='uniform', 
                      weights = dataframe$num_people)

  results <- c(
  "LATE_other" = rd_other$coef[1],
  "LATE_other_lcl" = rd_other$ci[1,1],
  "LATE_other_ucl" = rd_other$ci[1,2],
  "LATE_other_p" = rd_other$pv[1],
  "other_discontinuity" = rd_other$tau_cl[2]-rd_other$tau_cl[1],

  "LATE_covid" = rd_covid$coef[1],
  "LATE_covid_lcl" = rd_covid$ci[1,1],
  "LATE_covid_ucl" = rd_covid$ci[1,2],
  "LATE_covid_p" = rd_covid$pv[1],
  "LATE_covid_se" = rd_covid$se[1],

  "covid_discontinuity" = rd_covid$tau_cl[2]-rd_covid$tau_cl[1],
  "treat_discontinuity" = rd_covid$beta_p_r[1,2]-rd_covid$beta_p_l[1,2],
  "treat_upper" = rd_covid$beta_p_l[1,2],
  "treat_lower" = rd_covid$beta_p_r[1,2],
  "treat_from_treat_rdd" = rd_treat$coef[1],
  "treat_LCL" = rd_treat$ci[1,1],
  "treat_UCL" = rd_treat$ci[1,2],
  "treat_se" = rd_treat$se[1])
  return(results)
}

#----------------------------------
# different cut offs
#----------------------------------
cut_off=80*12

bandwidth = 60
results <- rdd_with_params(cut_off, bandwidth, plot_data)
results

cut_offs <- seq(910, 1010, 10)
results <- sapply(cut_offs, rdd_with_params, bandwidth=bandwidth, 
                  dataframe=plot_data)
results <- rbind(results)

results_df <- data.frame("cut_off" = cut_offs,
                         "LATE_other"=results[1,], 
                        "LATE_other_lcl" = results[2,],
                        "LATE_other_ucl" = results[3,],
                        "LATE_other_p" = results[4,],
                        "other_discontinuity" = results[5,],
                         
                        "LATE_covid" = results[6,],
                        "LATE_covid_lcl" = results[7,],
                        "LATE_covid_ucl" = results[8,],
                        "LATE_covid_p" = results[9,],
                        "LATE_covid_se" = results[10,],

                        "covid_discontinuity" = results[11,],
                        "treat_discontinuity" = results[12,],
                        "treat_upper" = results[13,],
                        "treat_lower" = results[14,],
                        "treat_from_treat_rdd" = results[15,],
                        "treat_LCL" = results[16,],
                        "treat_UCL" = results[17,],
                        "treat_se" = results[18,])

results_df

write.csv(results_df, 
          paste0(out_dir_paper_results,"/sensitivity_analysis/cut_offs_", 
                 toString(t_start), "_follow_up_", toString(follow_up), "_new_follow_up.csv"),
         row.names = FALSE)
results_df

#----------------------------------
# different bandwidths
#----------------------------------

cut_off=80*12
bandwidths = seq(10, 60, 10)

results <- sapply(X=bandwidths, rdd_with_params, cut_off=cut_off, 
                  dataframe=plot_data)
results <- rbind(results)

results_df <- data.frame("bandwidths" = bandwidths,
                         "LATE_other"=results[1,], 
                        "LATE_other_lcl" = results[2,],
                        "LATE_other_ucl" = results[3,],
                        "LATE_other_p" = results[4,],
                        "other_discontinuity" = results[5,],
                         
                        "LATE_covid" = results[6,],
                        "LATE_covid_lcl" = results[7,],
                        "LATE_covid_ucl" = results[8,],
                        "LATE_covid_p" = results[9,],
                        "LATE_covid_se" = results[10,],

                        "covid_discontinuity" = results[11,],
                        "treat_discontinuity" = results[12,],
                        "treat_upper" = results[13,],
                        "treat_lower" = results[14,],
                        "treat_from_treat_rdd" = results[15,],
                        "treat_LCL" = results[16,],
                        "treat_UCL" = results[17,],
                        "treat_se" = results[18,])

write.csv(results_df, 
          paste0(out_dir_paper_results,"/sensitivity_analysis/bandwidths_", 
                 toString(t_start), "_follow_up_", toString(follow_up), 
                 "_new_follow_up.csv"),
         row.names = FALSE)
results_df

#------------------------------
# doughnut rdd
#------------------------------



doughnut_rdd <- function(hole_half_width, cut_off, bandwidth, full_data){
  lower <- 960-hole_half_width
  higher <- 960+ hole_half_width
  print(lower)
  print(higher)

  doughnut_data <- full_data[(full_data$age <= lower | full_data$age >= higher),]
  
  results <- rdd_with_params(cut_off, bandwidth, doughnut_data)
  
  return(results)
}

test <- plot_data[(plot_data$age <= 959 | plot_data$age >= 961),]
rdd_with_params(cut_off, bandwidth, test)

cut_off <- 80*12
bandwidth <- 60
hole_half_width <- 0:50

results <- sapply(hole_half_width, doughnut_rdd, bandwidth = bandwidth, 
                  cut_off=cut_off, full_data=plot_data)
results <- rbind(results)

results_df <- data.frame("hole_half_widths" = hole_half_width,
                         "LATE_other"=results[1,], 
                        "LATE_other_lcl" = results[2,],
                        "LATE_other_ucl" = results[3,],
                        "LATE_other_p" = results[4,],
                        "other_discontinuity" = results[5,],
                         
                        "LATE_covid" = results[6,],
                        "LATE_covid_lcl" = results[7,],
                        "LATE_covid_ucl" = results[8,],
                        "LATE_covid_p" = results[9,],
                        "LATE_covid_se" = results[10,],

                        "covid_discontinuity" = results[11,],
                        "treat_discontinuity" = results[12,],
                        "treat_upper" = results[13,],
                        "treat_lower" = results[14,],
                        "treat_from_treat_rdd" = results[15,],
                        "treat_LCL" = results[16,],
                        "treat_UCL" = results[17,],
                        "treat_se" = results[18,])

write.csv(results_df, 
          paste0(out_dir_paper_results,"/sensitivity_analysis/doughnut_rdd_", 
                 toString(t_start), "_follow_up_", toString(follow_up), "_new_follow_up.csv"),
         row.names = FALSE)


#######
# END #
#######











