
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

setwd("your/working/directory")
 save_filepath <- "your/save/filepath"

#Prepare data for use in XGB training model
#Attach top 100 labs (not eventually used) and medication-route combinations.
load_data <- function() {
    baseline_factors <- read.csv("preprocessed_demographics_for_epi.csv")

    linked_stays_for_epi <- read.csv("linked_stays_for_epi.csv")
    linked_stays_for_epi <- linked_stays_for_epi %>% select(all_of(c("csn", "hour_of_arrival", "day_of_arrival", "month_of_arrival", "year_of_arrival"))) 
    baseline_factors <- inner_join(baseline_factors, linked_stays_for_epi, by="csn") %>% filter(is_admitted==1 | is_discharged==1)


    labs <- read.csv("/Volumes/chip-lacava/Groups/BCH-ED/HC-or-work-PSB-run/labs.csv")
    medications <- read.csv("/Volumes/chip-lacava/Groups/BCH-ED/HC-or-work-PSB-run/medications.csv")

    scores <- read.csv("/Volumes/chip-lacava/Groups/BCH-ED/HC-or-work-PSB-run/scores.csv")
    scores <- scores %>% filter(SCORE_VARIABLE=="NRS Generalized Pain Score")
    vitals <- read.csv("/Volumes/chip-lacava/Groups/BCH-ED/HC-or-work-PSB-run/vitals.csv")
    radiology <- read.csv("/Volumes/chip-lacava/Groups/BCH-ED/HC-or-work-PSB-run/radiology_15Apr.csv")
    demographics <- read.csv("/Volumes/chip-lacava/Groups/BCH-ED/HC-or-work-PSB-run/demographics_8May.csv")

    #Link CSNs for radiology
    info_to_csn_df <- demographics %>% select(NAME, DOB, ED_ARRIVAL_TIME, CSN) %>% distinct(NAME, DOB, ED_ARRIVAL_TIME, CSN)
    radiology <- radiology %>% inner_join(info_to_csn_df, by=c("NAME", "DOB", "ED_ARRIVAL_TIME")) #now radiology has CSNs


    decisions <- data.frame(csn=character(0), event_type=character(0), event_name=character(0), event_time=character(0), lab_id=numeric(0))
    results <- data.frame(csn=character(0), event_type=character(0), event_name=character(0), event_time=character(0), event_value=numeric(0), lab_id=numeric(0))

    #First, add in arrival times.
    arrival_times <- baseline_factors %>% select(csn, ed_arrival_time) %>% rename(event_time=ed_arrival_time) %>% mutate(event_type="arrival", event_name="arrival", lab_id=NA)
    decisions <- rbind(decisions, arrival_times)


    #Get departure times.
    admission_times<- baseline_factors %>% filter(is_admitted==1) %>% select(csn, ed_checkout_time) %>% rename(event_time=ed_checkout_time) %>% mutate(event_type="endpoint", event_name="admission", lab_id=NA)
    decisions <- rbind(decisions, admission_times)

    discharge_times<- baseline_factors %>% filter(is_discharged==1) %>% select(csn, ed_checkout_time) %>% rename(event_time=ed_checkout_time) %>% mutate(event_type="endpoint", event_name="discharge", lab_id=NA)
    decisions <- rbind(decisions, discharge_times)


    #Get 100 most common labs and their results.
    freq_labs <- labs %>% group_by(LAB) %>% summarise(count=n()) %>% arrange(desc(count)) 
    labs <- labs %>% filter(LAB %in% freq_labs$LAB[1:100])

    lab_ordered_times <- labs %>% select(CSN, LAB_START_TIME, LAB) %>% rename(csn=CSN, event_time=LAB_START_TIME, event_name=LAB) %>% mutate(event_type="lab_ordered", lab_id=row_number()) #lab ID allows linking of labs and results
    decisions <- rbind(decisions, lab_ordered_times)

    lab_result_times <- labs %>% select(CSN, LAB_DATE_TIME, LAB, VALUE) %>% rename(csn=CSN, event_time=LAB_DATE_TIME, event_name=LAB, event_value=VALUE) %>% mutate(event_type="lab_result", lab_id=row_number()) #will be the same as the row number in previous step
    results <- rbind(results, lab_result_times)


    #Attach 100 most common combinations of medications and route.
    med_times <- medications %>% select(CSN, MED_DATE_TIME, MEDICATION, ROUTE) %>% mutate(med_route=paste0(MEDICATION, "_", ROUTE)) %>% rename(csn=CSN, event_time=MED_DATE_TIME, event_name=med_route) %>% select(-c(MEDICATION, ROUTE)) %>% mutate(event_type="medication_given", lab_id=NA)
    freq_meds <- med_times %>% group_by(event_name) %>% summarise(count=n()) %>% arrange(desc(count))
    med_times <- med_times %>% filter(event_name %in% freq_meds$event_name[1:100])
    decisions <- rbind(decisions, med_times)


    #Attach radiological tests.
    rad_times <- radiology %>% select(CSN, EVENT_END_DT_TM, RADIOLOGY_TEST) %>% rename(csn=CSN, event_time=EVENT_END_DT_TM, event_name=RADIOLOGY_TEST) %>% mutate(event_type="test_ordered", lab_id=NA)
    decisions <- rbind(decisions, rad_times)

    #Attach pain scores.
    pain_score_times <- scores %>% select(CSN, SCORE_DT_TIME, SCORE) %>% rename(csn=CSN, event_time=SCORE_DT_TIME, event_value=SCORE) %>% mutate(event_type="vitals_taken", event_name="PAIN", lab_id=NA)
    results <- rbind(results, pain_score_times)

    #Attach vitals.
    vital_times <- vitals %>% select(CSN, VITAL_DATE_TIME, VITAL, VITAL_VALUE) %>% rename(csn=CSN, event_time=VITAL_DATE_TIME, event_name=VITAL, event_value=VITAL_VALUE) %>% mutate(event_type="vitals_taken", lab_id=NA)
    vital_times$event_name <- toupper(vital_times$event_name)
    results <- rbind(results, vital_times)

    #Make sure we only have decisions/results taken during the visit.
    arrival_timestamps <- arrival_times %>% select(csn, event_time) %>% rename(arrival_time=event_time)

    decisions <- inner_join(decisions, arrival_timestamps, by="csn")
    decisions$event_time <- difftime(ymd_hms(decisions$event_time, truncated = 3), ymd_hms(decisions$arrival_time, truncated = 3), units="mins") #allowing for missing seconds
    endpoints <- decisions %>% filter(event_type == "endpoint") %>% select(csn, event_time) %>% rename(endpoint=event_time)
    decisions <- inner_join(decisions, endpoints, by="csn") %>% filter(event_time <= endpoint) %>% select(-endpoint)

    decisions <- decisions %>% filter(!(event_type=="arrival")) %>% group_by(csn) %>% arrange(event_time, .by_group = TRUE)

    write.csv(decisions, paste0(save_filepath, "patient-journey-decisions-raw-names.csv"))

    results <- inner_join(results, arrival_timestamps, by="csn")
    results$event_time <- difftime(ymd_hms(results$event_time, truncated = 3), ymd_hms(results$arrival_time, truncated = 3), units="mins") 
    results <- inner_join(results, endpoints, by="csn") %>% filter(event_time <= endpoint, event_time >= 0) %>% select(-endpoint) #discard vitals taken after checkout time and before timeframe of visit
    results <- results %>% select(-arrival_time) %>% group_by(csn)  %>% arrange(event_time, .by_group = TRUE)

    write.csv(results, paste0(save_filepath, "patient-journey-results-raw-names.csv"))
}





#Isolate the training data in 'baseline factors' and convert it to binary columns/continuous data
#this is done PER VISIT so there's no chance of contamination of visits
#this gets the first FOUR YEARS of data
get_train_and_test_data_by_date <- function(training_window=0) {

    converted_baseline_factors <- converted_baseline_factors %>% inner_join(arrival_timestamps, by="csn") 
    min_arrival_timestamp <- min(converted_baseline_factors$arrival_timestamp)
    training_data <- converted_baseline_factors %>% filter(arrival_timestamp < min_arrival_timestamp+training_window) %>% select(-arrival_timestamp)  #get data with the right chief complaint     
    test_data <- converted_baseline_factors %>% filter(arrival_timestamp >= min_arrival_timestamp+training_window)  %>% select(-arrival_timestamp)  #get data with the right chief complaint


    return(list(train=training_data, test=test_data))
}


#Generate timestamps every 30 mins during a visit.
generate_timestamps_at_intervals <- function(csn, max_time, interval=30) {
    max_timestamp <- floor(max_time/interval)
    timestamps <- (1:max_timestamp)*interval
    return(data.frame(csn=rep(csn, length(timestamps)), timestamp=timestamps))

}

#Generate timestamps every THIRTY MINUTES until the patient gets admitted
#Require patients to be in the ED for at least 30 mins
extract_timestamps_per_visit <- function(dataset, interval=30) {
    
    relevant_decisions <- decisions %>% filter(csn %in% dataset$csn, !(event_type=="endpoint")) 
    relevant_results <- results %>% filter(csn %in% dataset$csn)
    max_times <- decisions %>%  filter(event_type=="endpoint", csn %in% dataset$csn)
    visit_timestamps <- purrr::map2(max_times$csn, max_times$event_time, function(x, y) generate_timestamps_at_intervals(x, y, interval=interval)) %>% bind_rows()


    #Return these timestamps with baseline factors attached
    expanded_baseline_factors <- left_join(visit_timestamps, dataset, by="csn") #attach the 'point of triage completion' column along with this

    return(list(expanded_baseline_factors=expanded_baseline_factors, relevant_decisions=relevant_decisions, relevant_results=relevant_results))
}

#Attach vitals which have been taken at this point during a visit.
get_vitals_for_all_patients_and_timestamps <- function(csns_and_timestamps, relevant_results) {
    #get vitals
    relevant_vitals <- relevant_results %>% filter(event_type=="vitals_taken")
    up_to_date_vitals <- left_join(csns_and_timestamps, relevant_vitals, by="csn", relationship="many-to-many") %>% 
            filter(event_time <= timestamp)
    up_to_date_vitals$event_value <- as.numeric(up_to_date_vitals$event_value)

    up_to_date_vitals <- up_to_date_vitals %>% group_by(csn, timestamp, event_name) %>%  
            summarise(min = if (all(is.na(event_value))) NA else min(event_value, na.rm=TRUE), mean = if (all(is.na(event_value))) NA else mean(event_value, na.rm=TRUE), max = if (all(is.na(event_value))) NA else max(event_value, na.rm=TRUE)) %>% 
            tidyr::pivot_wider(id_cols=c(csn, timestamp), names_from=event_name, values_from=c(min, mean, max)) %>% clean_names()

    return(up_to_date_vitals)
}




#Attach medications dispensed at each timepoint- records a 1 if a medication was ever given, 0 otherwise
get_meds_for_all_patients_and_timestamps <- function(csns_and_timestamps, relevant_decisions, meds_to_use=c()) {

    #Identify all meds
    relevant_meds <- relevant_decisions %>% filter(event_type=="medication_given") 

    if (length(meds_to_use)>0) { #if we have been given specific meds to look for- not used.
        relevant_meds <- relevant_meds %>% filter(event_name %in% meds_to_use)
    }

    up_to_date_meds <- left_join(csns_and_timestamps, relevant_meds, by="csn", relationship="many-to-many") %>% 
            filter(event_time <= timestamp) %>% group_by(csn, timestamp, event_name) %>% summarise(given=1) %>%
            tidyr::pivot_wider(id_cols=c(csn, timestamp), names_prefix="med_", names_from=event_name, values_from=given, values_fill=0) 

    return(up_to_date_meds)
}

#Now get tests for all patients and timestamps-- records a 1 if a test was ever ordered (up to this point), 0 otherwise
get_tests_for_all_patients_and_timestamps <- function(csns_and_timestamps, relevant_decisions, tests_to_use=c()) {

    relevant_tests <- relevant_decisions %>% filter(event_type=="test_ordered")

    if (length(tests_to_use)>0) { #if we have been given specific tests to check--not used
        relevant_tests <- relevant_tests %>% filter(event_name %in% tests_to_use)
    }

    #Now get tests for each patient and timestamp
    up_to_date_tests <- left_join(csns_and_timestamps, relevant_tests, by="csn", relationship="many-to-many") %>% 
            filter(event_time <= timestamp) %>% group_by(csn, timestamp, event_name) %>% summarise(ordered=1) %>%
            tidyr::pivot_wider(id_cols=c(csn, timestamp), names_prefix="test_ever_ordered_", names_from=event_name, values_from=ordered, values_fill=0) 

    return(up_to_date_tests)
}




#Normalise HR and RR by age group
normalise_vitals <- function(linked_data, vitals_to_normalise=c("respiratory_rate", "heart_rate")) {
    all_colnames <- colnames(linked_data)
    vital_cols <- all_colnames[grepl(paste(vitals_to_normalise, collapse="|"), all_colnames)] #check if any of these vitals are in a col
    age_groups <- unique(linked_data$age_group)
    #print(age_groups)
    for (group in age_groups) {
        relevant_indices <- which(linked_data$age_group==group)
        for (col in vital_cols) {
            relevant_vitals <- as.numeric(unlist(linked_data[relevant_indices, col]))
            mu <- mean(relevant_vitals, na.rm=TRUE)
            sigma <- sd(relevant_vitals, na.rm=TRUE) #get mean and sd of each vital
            linked_data[relevant_indices, col] <- (relevant_vitals - mu)/sigma #do NOT take absolute value-- we can learn from positive and negative values here
        }
    }
    linked_data <- linked_data %>% select(-age_group)
    return(linked_data)
}


#Now put all of this together-- link meds and vitals to a training or testing dataset
#when converting the testing data, we'll only use labs, meds and tests found in both datasets
#Labs are only available in the last two years of dataset, so if we were pulling them they'd be eliminated by that step.
link_visit_data <- function(dataset, interval=30, vitals_to_normalise=c("respiratory_rate", "heart_rate"), meds_to_use=c(), tests_to_use=c()) {

    #Extract relevant data
    res <- extract_timestamps_per_visit(dataset, interval=interval)
    expanded_baseline_factors <- res$expanded_baseline_factors
    relevant_decisions <- res$relevant_decisions
    relevant_results <- res$relevant_results

    #Reduce expanded baseline factors to the needed arguments
    csns_and_timestamps <- expanded_baseline_factors %>% select(csn, timestamp)
    
    up_to_date_vitals <- get_vitals_for_all_patients_and_timestamps(csns_and_timestamps, relevant_results)
    up_to_date_meds <- get_meds_for_all_patients_and_timestamps(csns_and_timestamps, relevant_decisions, meds_to_use=meds_to_use)
    up_to_date_tests <- get_tests_for_all_patients_and_timestamps(csns_and_timestamps, relevant_decisions, tests_to_use=tests_to_use)

    linked_data <- left_join(expanded_baseline_factors, up_to_date_vitals, by=c("csn", "timestamp"))
    linked_data <- left_join(linked_data, up_to_date_meds, by=c("csn", "timestamp"))
    linked_data <- left_join(linked_data, up_to_date_tests, by=c("csn", "timestamp"))


    #Normalise vitals
    linked_data <- normalise_vitals(linked_data, vitals_to_normalise=vitals_to_normalise)

    #correct NAs to 0s 
    med_cols <- setdiff(colnames(up_to_date_meds), colnames(expanded_baseline_factors))
    ever_ordered_cols <- colnames(linked_data)[grepl("ever_ordered", colnames(linked_data))] #this will also get NAs in tests
    for (col in c(med_cols, ever_ordered_cols)) {
        linked_data[[col]] <- ifelse(is.na(linked_data[[col]]), 0, linked_data[[col]]) 
    }


    return(linked_data)

}

#Get linked datasets with top 100 meds and tests
compile_linked_datasets_by_date <- function(interval=30, training_window=0, training_window_name="", vitals_to_normalise=c("respiratory_rate", "heart_rate")) {
    
    data <- get_train_and_test_data_by_date(training_window=training_window)

    #Load both sets.
    linked_training_data <- link_visit_data(data$train, interval=interval, vitals_to_normalise=vitals_to_normalise) %>% arrange(csn, timestamp)
    linked_test_data <- link_visit_data(data$test, interval=interval, vitals_to_normalise=vitals_to_normalise) %>% arrange(csn, timestamp)

    #Keep variables which appear in both.
    cols_to_keep_in_order <- intersect(colnames(linked_training_data), colnames(linked_test_data))
    linked_test_data <- linked_test_data[cols_to_keep_in_order] #Make sure these are in the same order
    linked_training_data <- linked_training_data[cols_to_keep_in_order]


    write.csv(linked_training_data, paste0(save_filepath, "training_window_", training_window_name, "_training_data.csv"))
    write.csv(linked_test_data, paste0(save_filepath, "training_window_", training_window_name, "_test_data.csv"))

    
}

#Identify labs and med/route combinations, and attach them to patient journeys.
load_data()

#Load timestamps, and filter out visits for which no decision was made.
arrival_timestamps <- read.csv("baseline-factors-with-top-200-complaint-stems-tagged.csv") %>% filter(is_admitted==1 | is_discharged==1) %>% select(csn, arrival_timestamp)
converted_baseline_factors <- read.csv("cvtd-factors-with-tagged-complaint.csv")
decisions <- read.csv("patient-journey-decisions-raw-names.csv")
results <- read.csv("patient-journey-results-raw-names.csv")


#Render lowercase and replace so all event_names can be easily used as column names
decisions$event_name <- tolower(decisions$event_name)
decisions$event_name <- gsub("[^a-z0-9_]", "_", decisions$event_name) # replace non-alphanumeric with underscore
decisions$event_name <- gsub("_{2,}", "_", decisions$event_name) # replace multiple underscores with single
decisions$event_name <- gsub("^_|_$", "", decisions$event_name) # remove leading/trailing underscores



results$event_name <- tolower(results$event_name)
results$event_name <- gsub("[^a-z0-9_]", "_", results$event_name) # replace non-alphanumeric with underscore
results$event_name <- gsub("_{2,}", "_", results$event_name) # replace multiple underscores with single
results$event_name <- gsub("^_|_$", "", results$event_name) # remove leading/trailing underscores




#Get first 4 years of data
compile_linked_datasets_by_date(training_window = 60*24*365.25*4, training_window_name = "four_years")
