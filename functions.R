#--------------------------
# Utility functions
#--------------------------

df_length <- function(df){
  # calculates the length (number of rows) of a sparklyr dataframe
  length <- df %>%
  dplyr::count() %>%
  collect()
  
  return(length)
}

#------------------------------
# Descriptive stats functions
#------------------------------



get_tables_per_var <- function(df, vars){
  
  res <- list()
  
  for (v in vars){
    print(v)
    
    v <- as.symbol(v)
    v_str <- as.character(v)
    
    df_var <- df %>%
      select(v)
    
     if (!(v_str %in% continuous_vars)){
      df_summary <- df_var %>%
        select(v) %>%
        group_by(v) %>%
        tally() %>%
        collect() %>%
        as.data.frame()%>%
        arrange()

      df_summary$percentage <- df_summary$n / sum(df_summary$n) * 100
  
      df_summary$percentage <- round(df_summary$percentage, 2)

      df_summary$n <- format(df_summary$n, big.mark=",", scientific=FALSE)
      df_summary[[v_str]] <- as.character(df_summary[[v_str]])
      
#      new_col <- rep('temp', nrow(df_summary))
#      n <- df_summary$n
#      p <- df_summary$p
#      print( df_summary$p)
#      for (i in 1:length(new_col)){
#        new_col[i] <- paste0(n[i], ' (', format(p[i], nsmall = 2), ')')
#      }
#      
      df_summary[['Count (n)']] <- paste0(df_summary$n, ' (', format(df_summary$percentage, nsmall = 2), ')')
      print(df_summary )
      df_summary <- select(df_summary, -n, -percentage)
      
      if (v_str %in% binary_vars){
        df_summary <-df_summary %>%
        filter()
      }
      
    }
    
    else if (v_str %in% continuous_vars){
      df_summary <- df_var %>%
      sdf_describe(cols=v_str) %>%
      filter(summary=='stddev' | summary=='mean') %>%
      collect() %>%
      arrange() %>%
      t() %>%
      as.data.frame()
      
      colnames(df_summary) <- c('mean', 'sd')
      df_summary <- df_summary[-1,]
      df_summary$sd <- as.numeric(as.character(df_summary$sd))
      df_summary$sd <- round(df_summary$sd, 2)
      df_summary$sd <- format(df_summary$sd, big.mark=",", scientific=FALSE)
      df_summary$mean <- as.numeric(as.character(df_summary$mean))
      df_summary$mean <- round(df_summary$mean, 2)
      df_summary$mean <- format(df_summary$mean, big.mark=",", scientific=FALSE)
      
      
      new_col <- rep('temp', nrow(df_summary))
      mean <- df_summary$mean
      sd <- df_summary$sd
      
      for (i in 1:length(new_col)){
        new_col[i] <- paste0(mean[i], ' (', sd[i], ')')
      }
      
      df_summary[['Mean (SD)']] <- new_col
      df_summary <- select(df_summary, -mean, -sd)
      
    }
    
    df_summary$variable <- rep(v_str, nrow(df_summary))
    
    res <- append(res, list(df_summary))
    
    rm(df_var, df_summary)
    gc()
    
  }
  
  return(res)
}


# function to format the tables
format_table<- function(data, name='Count (n)'){
  ind <- c(sapply(continuous_vars, function(x) which(columns == x)[[1]]))
  cat <- data[-ind]
  con <- data[ind]


  for (i in 1:length(cat)){
    print(colnames(cat[[i]]))
    colnames(cat[[i]]) <- c('value', name, 'variable')
    cat[[i]]<- arrange(cat[[i]],value)
  }

  cat <- do.call('rbind', cat) 
  cat <- cat %>%
    left_join(var_lookup)%>%
    group_by(variable)%>%
    filter(is.na(type)|type == ""|type=="binary"&as.numeric(value) ==1)%>%
    as.data.frame()
  
    con <- do.call('rbind', con)%>%
  mutate(value="", labels="", type="")
names(con) <-c(name, 'value', "variable" , "labels" ,   "type" )

rbind(cat, con)
  
}

#--------------------------
# Weighting functions
#--------------------------



prop_vacc_by_type <- function(delay, dataframe_analysis_period, start_date, end_date){
  print(delay)
  if ((as.Date("2020-12-08") + delay) > end_date){
    print("no vaccinations withthis delay before end of analysis period")
    fit_df_age <- data.frame()}
  else{
  df_vacc <- dataframe_analysis_period %>%
    mutate(vacc_plus = date_add(first_vacc_date,delay)) %>%
    mutate(vacc_plus = ifelse(vacc_plus > end_date, NA, vacc_plus)) %>%
#    mutate(vacc_plus = ifelse(vacc_plus < start_date, NA, vacc_plus)) %>%
    mutate(event = case_when(!is.na(vacc_plus) & (vacc_plus<dod | is.na(dod)) & (first_vacc_manufacturer == 'Pfizer')  ~ 'Pfizer',
                          !is.na(vacc_plus) & (vacc_plus<dod | is.na(dod)) & (first_vacc_manufacturer == 'AstraZeneca')  ~ 'AstraZeneca',
                                   death==1 ~ 'Censure')) %>%
    mutate(t = case_when((event == 'Pfizer' | event == 'AstraZeneca') ~ datediff(vacc_plus, start_date),
           event == 'Censure' ~ datediff(dod, start_date),
            is.na(event) ~ datediff(end_date, start_date))) %>%
    mutate(t = ifelse(t<0, 0, t))

  # add a 0 instead of NA for people who are alive and unvaccinated
  # at the end of the time period so that they are censured
  df_vacc <- df_vacc %>%  
    na.replace(event = 'Censure')
  
  # create R dataframe for survival analysis
  data_df_vacc <- df_vacc %>%
    select(t, event, age_when_eligible_months) %>%
          collect()

  data_df_vacc$event <- factor(data_df_vacc$event, levels = c("Censure", "Pfizer", "AstraZeneca"))

  # fit dataframe using A-J survival model
  mfit<- survfit(Surv(t, as.factor(event), type="mstate") ~ age_when_eligible_months, 
              data=data_df_vacc)

  #extract data from fitted model
  fit_df_age <- data.frame(t = summary(mfit)$time,
                     n.risk= summary(mfit)$n.risk[,1],
                    Pfizer = summary(mfit)$pstate[,2],
                    AstraZeneca = summary(mfit)$pstate[,3],
                    age = summary(mfit)$strata)

  # tranform age to an integer number of months
  fit_df_age$age <- as.integer(gsub('age_when_eligible_months=', '', fit_df_age$age))}
  
  # add column for delay
  fit_df_age <- fit_df_age %>%
    mutate(min_vacc_days = paste0("day_", as.integer(delay)))
  
  return(fit_df_age) 
}



wide_form_prop_vaccinated <- function(dataframe, vaccine, vaccine_drop){

  # split into AZ and PF dataframes and convert to wideform
  cum_inc_all_days <- dataframe %>%
    select(-vaccine_drop) %>%
    rename(prop_vaccinated = vaccine)

  df_all_delay <- cum_inc_all_days %>%
    select(t, prop_vaccinated, min_vacc_days, age) %>%
    pivot_wider(names_from = min_vacc_days, values_from = prop_vaccinated) 

  # order by t in age in month groups
  df_all_delay <- df_all_delay[order(df_all_delay$age, df_all_delay$t),]

  # for missing values, fill down (as the prop vacicnated does not change on these days)
  # add zeros at the beginning if no values
  df_all_delay <- df_all_delay %>%
    group_by(age) %>%
    fill(-t, -age) %>%
    ungroup()
  
  df_all_delay[is.na(df_all_delay)] <- 0
  
  return (df_all_delay)
  
}

weighted_prop_treated_df_modelled <- function(df_all_delay, weight_data){
  days <- substring(names(select(df_all_delay, -age, -t)), first = 5)

  max_day <- tail(days, 1)

  # get the prop treated the max day times the value of sigmoid for this number of days
  weight <- df_all_delay[[paste0("day_", as.character(max_day))]]*weight_data$prediction[weight_data$time==max_day]

  if (length(days)>1){
    # add the additional proportion vaccinated at least max-1 days * sigmoid for this day,
    # then continue to the minimum
    for (i in 1:(length(days)-1)){
      day = days[length(days)-i]
      weight <- weight + (df_all_delay[[paste0("day_", as.character(day))]] - df_all_delay[[paste0("day_", as.character(days[length(days)-i+1]))]])*weight_data$prediction[weight_data$time==day]

  }}

  # add the weighted values as a column and select columns needed
  df_all_delay$weighted_prop_treated <- weight  
  df_all_delay <- df_all_delay %>%
    select(t, age, weighted_prop_treated)
 
  return(df_all_delay)
}

weighted_average_vacc_both <- function(age_months, prop_treated_dataframe, follow_up) {

  ## for each age in months, if it is in the dataframe, filter dataframe to just this age
  # add a column for the difference between the subsequent proportions vaccinated (with the
  # first being the first proportion). This gives the additional proportion of people
  # in the treatment state on each subsequent day.
  ## add a column for the weighting which is the proportion of the follow up period spent
  # in that treatment state so the number of days left in the follow up after t divided
  # by the follow up
  ## add a new column that multiplies the weight by the additional proportion of people
  # that enter the treatment status on that day, and adds this to the previous day
  # to give a measure of the weighted proportion treated by each day

  if (age_months %in% prop_treated_dataframe$age) {
     fit_df_age_specific <- prop_treated_dataframe %>%
      filter(prop_treated_dataframe$age==age_months) 
    
    fit_df_age_specific <- fit_df_age_specific %>%
      mutate(p_non_cum = c(fit_df_age_specific$weighted_prop_treated[1], diff(fit_df_age_specific$weighted_prop_treated)),
        weight = ((follow_up-1)-t)/(follow_up-1)) %>%
      mutate(weighted_p_non_cum = cumsum(weighted_p_non_cum = p_non_cum*weight))

  # the weighted proportion treated over the follow period is the cumulative sum
  # on the last day of the follow up period
#  max_t <- max(fit_df_age_specific$t) 
  cum_sum_age <- fit_df_age_specific %>%
    filter(t==follow_up-1) 
  weighted_p <- cum_sum_age$weighted_p_non_cum
    
    # if there are no entries for that age, return 0 as no people of that age have been
    # treated during the follow up period
  } else {weighted_p <- 0}

  # return the weighted proportion pf people treated over the follow up period
  return(weighted_p)
}







#--------------------------
# RDD functions
#--------------------------


fit_rdd <- function(dataframe, assignment_col, outcome_col, c){
  # function to fit a linear rdd and calculate the value of the discontinuity
  
  # add a column for a dummy variable that is 1 above the cutoff and 0 below
  dataframe <- dataframe %>% 
    dplyr::mutate(D = ifelse(dataframe[[assignment_col]] >= c, 1, 0))
  
  # formula for linear rdd
  formula = paste0(outcome_col, " ~ D * I(", assignment_col, " - ", toString(c), ")")

  # fit rdd
  fit <- lm(formula=formula, data=dataframe, weights=num_people)
  print(summary(fit))
  confints = confint(fit, 'D', level=0.95)
  
  # get the main effect size and its error
  main_effect <- summary(fit)$coef[2,1]
  main_effect_error <- summary(fit)$coef[2,2]
  p <- summary(fit)$coefficients[2,4]
  return(list("summary" = summary(fit), "main_effect" = main_effect, 
              "main_effect_error" = main_effect_error, "fit"=fit,
              "p_value"=p,
             "rdd_dataframe" = dataframe,
             "LCL"=confints[1],
             "UCL"=confints[2]))
}




plot_rdd <- function(dataframe, assignment_col, outcome_col, c, fit, group_col, save,
                    title, x_lab, y_lab, save_name){
  # plot data used for rdd and the relevant rdd fit

  # calculate predicted lines using the rdd fit either side of the cut off
  line_higher <- predict(fit, newdata = filter(dataframe, dataframe[[group_col]]==1), interval="confidence",
               level=0.95)
  
  line_lower <- predict(fit, newdata = filter(dataframe, dataframe[[group_col]]==0), interval="confidence",
               level=0.95)

  # create dataframes for data above and below cut off together with the predicted
  # lines and errors
  data_lower = data.frame(age=filter(dataframe,dataframe[[group_col]]==0)[[assignment_col]], 
                          prediction = line_lower[,"fit"], 
                        upper = line_lower[,"upr"], lower = line_lower[,"lwr"])
  data_higher = data.frame(age=filter(dataframe,dataframe[[group_col]]==1)[[assignment_col]], 
                           prediction = line_higher[,"fit"], 
                        upper = line_higher[,"upr"], lower = line_higher[,"lwr"])
  
  # set the dummy variable column as a factor and rename
  dataframe[[group_col]] <- as.factor(dataframe[[group_col]])
  dataframe[[group_col]] <- revalue(dataframe[[group_col]], c("0"="Control group", "1"="Treatment group"))
  

  # plot data and fit for rdd
  plot <- ggplot() +
  geom_point(data=dataframe, aes(x = .data[[assignment_col]]/12, y = .data[[outcome_col]], color=.data[[group_col]]))+ 
  geom_line(data = data_lower, aes(x=age/12, y=prediction)) +
  geom_ribbon(data = data_lower, aes(x=age/12, ymin=lower, ymax=upper), linetype=2, alpha=0.3)+
  geom_line(data = data_higher, aes(x=age/12, y=prediction)) +
  geom_ribbon(data = data_higher, aes(x=age/12, ymin=lower, ymax=upper), linetype=2, alpha=0.3) +
  geom_vline(xintercept=c/12) +
#  geom_segment(aes(x = 80, y = 0.003, xend = 80, yend = 0.002, colour=3),
#                  arrow = arrow(length = unit(0.5, "cm")))+
  labs(title= title, x=x_lab, y=y_lab)# + theme(text = element_text(size = 20))

  # save plot if needed 
  if (save==TRUE){
    ggsave(save_name, height=6, width=4/3*6)
  
  return(plot)
  }

}

plot_rdd_data <- function(dataframe, assignment_col, outcome_col, 
                          c, fit, group_col, suffix){
  # plot data used for rdd and the relevant rdd fit

  # calculate predicted lines using the rdd fit either side of the cut off
  line_higher <- predict(fit, newdata = filter(dataframe, dataframe[[group_col]]==1), interval="confidence",
               level=0.95)
  
  line_lower <- predict(fit, newdata = filter(dataframe, dataframe[[group_col]]==0), interval="confidence",
               level=0.95)

  # create dataframes for data above and below cut off together with the predicted
  # lines and errors
  data_lower = data.frame(age=filter(dataframe,dataframe[[group_col]]==0)[[assignment_col]], 
                          prediction = line_lower[,"fit"], 
                        upper = line_lower[,"upr"], lower = line_lower[,"lwr"])
  data_higher = data.frame(age=filter(dataframe,dataframe[[group_col]]==1)[[assignment_col]], 
                           prediction = line_higher[,"fit"], 
                        upper = line_higher[,"upr"], lower = line_higher[,"lwr"])
  
  # set the dummy variable column as a factor and rename
  dataframe[[group_col]] <- as.factor(dataframe[[group_col]])
  dataframe[[group_col]] <- revalue(dataframe[[group_col]], c("0"="Control group", "1"="Treatment group"))
  
  all_prediction <- rbind(data_lower, data_higher)
  
  return(all_prediction)
}

survival_df <- function(df_follow_up, start_date, end_date){
  # make a dataframe for survival analysis of covid and non covid deaths
  # with a particular start day (so start with everyone alive on this day)
#  cut_off_date = as.character(as.Date("20201208", format = "%Y%m%d") +t_start) 

#  # remove people who died before cutoff t_start
#  df_follow_up <- df %>%
#  filter(dod>=cut_off_date | is.na(dod))

  # dataframe for multistate survival analysis
  # with events censure/Other death/Covid deaths and the time to the event
  # from the cut off date
  df_follow_up <- df_follow_up %>%
    dplyr::mutate(event_post_cut_off = case_when(death == 0 ~ "Censure",
                                    covid_death_status == 0 ~ "Other death",
                                    covid_death_status == 1 ~ "COVID-19 death")) %>%                               
    dplyr::mutate(t_post_cut_off = case_when(death == 1 ~ datediff(dod, start_date),
                                     event_post_cut_off=="Censure" ~ datediff(end_date, start_date)))
  return(df_follow_up)
}



#--------------------------
# Model
#--------------------------


full_rd_model <- function(t_start, follow_up, df){
  
  start_time_sys <- Sys.time()
  
  #--------------------------
  # Analysis period dataset
  #--------------------------
  start_date = as.character(as.Date("20201208", format = "%Y%m%d") +t_start) 
  print(start_date)
  end_date = as.Date(start_date) + follow_up - 1
  print(end_date)
  
  print(t_start)
  print(follow_up)
  
  # remove people who died before the analysis period
  # set deaths variables for deaths that occur after the analysis period to NULL/0 as required
  # set 'unknown' and Moderna vaccination to 'Pfizer'
  df_follow_up <- df %>%
    filter((dod>=start_date) | is.na(dod)) %>%
    select(age_when_eligible_months, dod, first_vacc_manufacturer, first_vacc_date,
           death, covid_death_status, second_vacc_date)
  
  # #  number of second vaccinations
  #  test <- df_follow_up %>%
  #    filter(!is.na(second_vacc_date) & second_vacc_date <= end_date) %>%
  #    count() %>%
  #    collect()
  #  test
  #  
  #  n_75_80 <- df_follow_up %>%
  #    filter(age_when_eligible_months >= 900 & age_when_eligible_months < 960) %>%
  #    count() %>%
  #    collect()
  #  n_75_80
  
  # number of people in the dataset for each age in months
  n_by_age <- df_follow_up %>%
    group_by(age_when_eligible_months) %>%
    count() %>%
    rename(num_people = n) %>%
    collect()
  
  
  #--------------------------
  # Proportion vaccinated on each day of analysis period
  #--------------------------
  # people vacicnated before the analysis period are included in the total for day 0
  
  ### construct a dataframe of the cumulative incidence of each treatement definition
  # for 1,2,3,4 weeks using multistate survival curves for each vaccine type
  ### returns a list of dataframes, 
  # each dataframe is for a different delay since vaccination
  df_list <- lapply(seq(0, 28, 7), prop_vacc_by_type, 
                    dataframe_analysis_period=df_follow_up, start_date = start_date,
                    end_date = end_date)
  
  
  # combine the dataframes into one 
  cum_inc_all_days <- bind_rows(df_list)
  
  
  # convert to wide form (delays as columns) for AZ and Pfizer
  prop_vaccinated_Pfizer <- wide_form_prop_vaccinated(cum_inc_all_days, "Pfizer", "AstraZeneca")
  
  prop_vaccinated_AZ <- wide_form_prop_vaccinated(cum_inc_all_days, "AstraZeneca", "Pfizer")
  
  weighted_prop_treated_pfizer <- weighted_prop_treated_df_modelled(prop_vaccinated_Pfizer, 
                                                                    filter(prob_antibody_pos, 
                                                                           model=='pf_one_negative') )
  
  weighted_prop_treated_az <- weighted_prop_treated_df_modelled(prop_vaccinated_AZ, 
                                                                filter(prob_antibody_pos, 
                                                                       model=='az_one_negative') )
  
  total_prop_treated <- merge(weighted_prop_treated_pfizer,
                              weighted_prop_treated_az,
                              by=c('t', 'age')) %>%
    mutate(weighted_prop_treated = weighted_prop_treated.x + weighted_prop_treated.y) %>%
    select(-weighted_prop_treated.x, -weighted_prop_treated.y)
  
  print(unique(total_prop_treated$t))
  print(total_prop_treated)
  # calculate the weighted proportion of people treated during the follow up period
  weighted_cum_vacc_df <- data.frame(age =unique(total_prop_treated$age),
                                     weighted_cum_vacc = sapply(unique(total_prop_treated$age), 
                                                                weighted_average_vacc_both, 
                                                                prop_treated_dataframe=total_prop_treated,
                                                                follow_up=follow_up)) 
  
  
  weighted_cum_vacc_df <- weighted_cum_vacc_df %>%
    select(age, weighted_cum_vacc)
  
  
  #--------------------------
  # Dataframe for survival after cutoff for multistate models
  #--------------------------
  # Probability of covid/non covid mortality from a particular cut off
  # date as start date. Plot multistate survival curves.
  
  # t, event dataframe for covid and other deaths with t_start as day 0
  
  df_follow_up_mort <- df_follow_up %>%
    dplyr::mutate(event_post_cut_off = case_when(death == 0 ~ "Censure",
                                                 covid_death_status == 0 ~ "Other death",
                                                 covid_death_status == 1 ~ "COVID-19 death")) %>%                               
    dplyr::mutate(t_post_cut_off = case_when(death == 1 ~ datediff(dod, start_date),
                                             event_post_cut_off=="Censure" ~ datediff(eos_date, start_date)))
  
  
  
  # extract r dataframe
  data_df_follow_up <- df_follow_up_mort %>%
    select(t_post_cut_off, event_post_cut_off, age_when_eligible_months) %>%
    collect()
  
  
  #--------------------------
  # Fit survival for age in months
  #--------------------------
  
  mfit<- survfit(Surv(t_post_cut_off, as.factor(event_post_cut_off),type="mstate") ~ age_when_eligible_months, 
                 data=data_df_follow_up)
  
  fit_df <- data.frame(t = summary(mfit)$time,
                       n.risk= summary(mfit)$n.risk[,1],
                       covid_death = summary(mfit)$pstate[,2],
                       other_death = summary(mfit)$pstate[,3],
                       age = summary(mfit)$strata)
  
  fit_df$age <- as.integer(gsub('age_when_eligible_months=', '', fit_df$age))
  
  to_merge <- data.frame(t=rep(0:follow_up, each=120), age=rep(900:1019, follow_up+1))
  fit_df <- merge(to_merge, fit_df, all.x=TRUE)
  
  # order by t in age in month groups
  fit_df <- fit_df[order(fit_df$age, fit_df$t),]
  
  # for missing values, fill down
  # add zeros at the beginning if no values
  fit_df <- fit_df %>%
    group_by(age) %>%
    fill(-t, -age) %>%
    ungroup()
  
  fit_df[is.na(fit_df)] <- 0
  
  
  # set the day number to follow up to
  fit_df_time <- fit_df %>%
    filter(t==follow_up-1)
  
  #--------------------------
  # Fit both treatment and outcome in one to get LATE
  #--------------------------
  
  c=80*12
  
  # put treatemtn and outcome data in one dataframe
  data_treat_outcome <- merge(weighted_cum_vacc_df, 
                              select(fit_df_time, age, covid_death, other_death),
                              by.x='age', by.y='age', all=TRUE)
  
  # add number of people column for weighting
  data_treat_outcome <- merge(data_treat_outcome, 
                              n_by_age,
                              by.x='age', by.y='age_when_eligible_months', all.x=TRUE)
  
  
  
  # rename to simplify names and add dummy column
  data_treat_outcome <- data_treat_outcome %>%
    rename(p_treat = weighted_cum_vacc, p_covid = covid_death, p_other = other_death) %>% 
    dplyr::mutate(dummy = ifelse(age >= c, 1, 0))
  
  #  write.csv(data_treat_outcome,
  #          paste0(out_dir_paper_results, "/plot_data/plot_data_t_start_in_gpes_", toString(t_start),
  #                 "_follow_up_", toString(follow_up), ".csv"), row.names = FALSE)
  
  
  # calculate log odds ratios
  rdd <- data_treat_outcome %>% 
    mutate(log_odds_covid = log( p_covid / (1-p_covid)),
           log_odds_treat = log( p_treat / (1-p_treat)),
           log_odds_other = log( p_other / (1-p_other)))
  
  rdd <- rdd[rdd$p_covid != 0,]
  
  ### RDD models
  ## covid
  # log_odds (for VE)
  rdd_covid=rdrobust(y=rdd$log_odds_covid , x=rdd$age, c=960,
                     fuzzy = rdd$p_treat, p=1, h=60, kernel="uniform", weights = rdd$num_people)
  # rate (for LATE)
  rdd_covid_rate = rdrobust(y=rdd$p_covid , x=rdd$age, c=960,
                     fuzzy = rdd$p_treat, p=1, h=60, kernel="uniform", weights = rdd$num_people)
  ## non covid
  # log_odds (for VE)
  rdd_other=rdrobust(y=rdd$log_odds_other , x=rdd$age, c=960,
                     fuzzy = rdd$p_treat, p=1, h=60, kernel="uniform", weights = rdd$num_people)
  # Rate
  rdd_other_rate=rdrobust(y=rdd$p_other , x=rdd$age, c=960,
                     fuzzy = rdd$p_treat, p=1, h=60, kernel="uniform", weights = rdd$num_people)
  
  # discontinuity in treatment 
  rdd_treat=rdrobust(y=rdd$p_treat , x=rdd$age, c=960,
                     p=1, h=60, kernel="uniform", weights = rdd$num_people)
  summary(rdd_covid)
  
  rdd_covid
  
  # to get the discontinuity in mortality rates
  rdd_covid_ITT =rdrobust::rdrobust(y=rdd$p_covid , x=rdd$age, c=960,
                                    p=1, h=60, kernel="uniform", weights = rdd$num_people)
  
  rdd_non_covid_ITT =rdrobust::rdrobust(y=rdd$p_other , x=rdd$age, c=960,
                                    p=1, h=60, kernel="uniform", weights = rdd$num_people)
  
  
  
  #  ##--------------------------
  #  ## LATE
  #  ##--------------------------
  #
  LATE_other = rdd_other_rate$coef[1]
  LATE_other_lcl = rdd_other_rate$ci[1,1]
  LATE_other_ucl = rdd_other_rate$ci[1,2]
  LATE_other_p = rdd_other_rate$pv[1]
 # other_discontinuity <- rdd_other$tau_cl[2]-rdd_other$tau_cl[1]
  #
  LATE_covid = rdd_covid_rate$coef[1]
  LATE_covid_lcl = rdd_covid_rate$ci[1,1]
  LATE_covid_ucl = rdd_covid_rate$ci[1,2]
  LATE_covid_p = rdd_covid_rate$pv[1]
  LATE_covid_se = rdd_covid_rate$se[1]
  #
  covid_discontinuity <- rdd_covid_ITT$coef[1]
  covid_LCL <-rdd_covid_ITT$ci[1,1]
  covid_UCL <-  rdd_covid_ITT$ci[1,2]
  
  other_discontinuity	 <- rdd_non_covid_ITT$coef[1]
  other_LCL <- rdd_non_covid_ITT$ci[1,1]
  other_UCL <-  rdd_non_covid_ITT$ci[1,2]
  
  treat_discontinuity <- rdd_covid$beta_p_r[1,2]-rdd_covid$beta_p_l[1,2]
  treat_upper <- rdd_covid$beta_p_l[1,2]
  treat_lower <- rdd_covid$beta_p_r[1,2]
  treat_from_treat_rdd <- rdd_treat$coef[1]
  treat_LCL <- rdd_treat$ci[1,1]
  treat_UCL <- rdd_treat$ci[1,2]
  treat_se <- rdd_treat$se[1]
  
  
  
  ##--------------------------
  ## VE
  ##--------------------------
  
  ve = 1-exp(rdd_covid$coef[1])
  ve_ci = 1-exp(rdd_covid$ci[1,])
  
  paste0(round(ve,2) , " (", round(ve_ci[[1]], 2), " - ", round(ve_ci[[2]], 2), ")")
  
  
  #
  ##--------------------------
  ## Save out parameters and results
  ##--------------------------
  
  results <- c("day number"=t_start, 
               "follow up time" = follow_up, 
               "treatment_discontinuity"=treat_discontinuity,
               "treat_upper"=treat_upper,
               "treat_lower"=treat_lower,
               "covid_discontinuity"=covid_discontinuity, 
               "covid_LCL" = covid_LCL,
               "covid_UCL" = covid_UCL,
               
               "other_discontinuity"=other_discontinuity, 
               "other_LCL" = other_LCL,
               "other_UCL" = other_UCL,
               "LATE_covid"=LATE_covid, 
               "LATE_covid_LCL"=LATE_covid_lcl, 
               "LATE_covid_UCL"=LATE_covid_ucl, 
               "LATE_covid_p_value"=LATE_covid_p,
               "LATE_other"=LATE_other, 
               "LATE_other_LCL"=LATE_other_lcl, 
               "LATE_other_UCL"=LATE_other_ucl, 
               "LATE_other_p_value"=LATE_other_p,
               "treat_from_treat_rdd" = treat_from_treat_rdd,
               "treat_LCL" = treat_LCL,
               "treat_UCL" = treat_UCL,
               "treat_se"=treat_se,
               "LATE_covid_se"=LATE_covid_se,
               "VE" = ve,
               "VE_LCL" = ve_ci[[1]],
               "VE_UCL" = ve_ci[[2]])
  
  end_time_sys <- Sys.time()
  time_to_run <- end_time_sys - start_time_sys
  print(time_to_run)
  
  return(results)
}
