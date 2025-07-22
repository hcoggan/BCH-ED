
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(testit)
library(httpgd)
library(broom)
library(stringr)
library(patchwork)
library(forcats)
library(pROC)
library(ggsignif)
hgd()

library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)
library(janitor)
library(MatchIt)
library(data.table)
library(exact2x2)
library(cobalt)
library(zipcodeR)
library(xgboost)
library(glmnet)
library(cowplot)


baseline_factors <- read.csv("preprocessed_demographics_for_epi.csv")
#this gives us race, sex, insurance, whether or not patient is trans, arrival mode, triage start and complete time
#outcome, complaint, state, SDI score, miles travelled, weight (as absolute deviation from the age-relevant mean),
#visit history, and number of other patients in the hospital at the time

#it does not give us temporal variables
#also link triage vitals- these become part of BASELINE FACTORS-- assume 'baseline' is at triage

linked_stays_for_epi <- read.csv("linked_stays_for_epi.csv")

triage_vital_names <- c("hypotensive_at_triage", "fever_at_triage", "low_temp_at_triage", "cts_mean_ox_at_triage", "cts_mean_heart_rate_before_triage", "cts_mean_respiratory_rate_before_triage", "ord_max_pain_at_triage")  


linked_stays_for_epi <- linked_stays_for_epi %>% select(all_of(c("csn", "hour_of_arrival", "day_of_arrival", "month_of_arrival", "year_of_arrival", triage_vital_names))) 

baseline_factors <- inner_join(baseline_factors, linked_stays_for_epi, by="csn")


baseline_factors <- baseline_factors %>% select(-c(admit_source, admitted_service, mrn, ed_arrival_time, ed_checkout_time, length_of_stay_in_minutes, triage_start_dt_tm, triage_complete_dt_tm,
        age_in_years, sdi_quantile, miles_travelled, weight, time_to_triage,  departure_timestamp)) #keep both age group and continuous age, also arrival timestamp; discard all other 'categorised' versions of continuous variables


write.csv(baseline_factors, "baseline-factors-for-patient-journey.csv")


baseline_factors <- read.csv("baseline-factors-for-patient-journey.csv")
baseline_factors <- baseline_factors %>% select(-c(is_discharged, X.1, X)) %>% rename(cts_age_in_days=age_in_days) #so that age_in_days is correctly labelled

#age_group will be removed later
#is_admitted is necessary for labelling



#DISCONNECT triage vitals
triage_vital_names <- colnames(baseline_factors)
triage_vital_names <- triage_vital_names[grepl("at_triage", triage_vital_names) | grepl("before_triage", triage_vital_names)]
baseline_factors <- baseline_factors %>% select(-c(all_of(triage_vital_names)))

#CONNECT original complaints, disconnect complaint columns
demographics <- read.csv("demographics.csv")
complaint_cols <- colnames(baseline_factors)[startsWith(colnames(baseline_factors), "complaint_")]
baseline_factors <- baseline_factors %>% select(-c(all_of(complaint_cols)))
chief_complaints <- demographics %>% select(CSN, ED_COMPLAINT) %>% clean_names()
baseline_factors <- inner_join(baseline_factors, chief_complaints, by="csn")
baseline_factors$ed_complaint <- tolower(baseline_factors$ed_complaint)

write.csv(baseline_factors, "baseline_factors_with_og_complaint.csv")

df <- read.csv("baseline_factors_with_og_complaint.csv") %>% select(csn, ed_complaint)


write.csv(df, "csn-to-complaint.csv")

weight <- read.csv("demographics.csv") %>% select(CSN, WEIGHT_KG) %>% rename(csn=CSN, weight=WEIGHT_KG)
write.csv(weight, "weight-to-csn.csv")
