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
library(geosphere)
library(xgboost)
library(glmnet)
library(cowplot)
library(yardstick)
library(ggridges)
library(comorbidity)


# setwd("/Volumes/chip-lacava/Groups/BCH-ED/")
# loadpath <- "raw-data/"
# savepath <- "reprocessing/"


setwd("/Users/helenacoggan/Documents/OSXLAP12460/BCH-data-2025/temp-dict")
loadpath <- "raw-data/"
savepath <- "reprocessing/"



# Attach records of medications, labs, vitals, tests to BCH visits. Needs to be 
# pre-EPIC; get top 50 med/route combinations, top 50 labs (age-normalised results), 
# all tests, age-normalised vitals.

#First load visits.
visits <- as.data.frame(fread(paste0(savepath, "preprocessed-visits.csv")))

#Load original demographics files, for linkage purposes 
linkage_table <- as.data.frame(fread(paste0(loadpath, "LaCava_Demopgraphics_May1_V1.csv"))) %>%
    select(CSN, MRN, NAME, DOB, ED_ARRIVAL_TIME) %>% rename(csn=CSN, mrn=MRN, name=NAME, dob=DOB,  
    arrival_time=ED_ARRIVAL_TIME) %>% group_by(csn) %>% mutate(duplicated=n()>1) %>% 
    filter(!duplicated) %>% select(-duplicated) %>% ungroup()
    
print("Loaded first linkage table")

#Link departure times.
departure_times <- as.data.frame(fread((paste0(loadpath, "LaCava_Dispositions_Aug7.csv")))) %>%
    select(MRN, NAME, DOB, ED_CHECKIN_DT_TM, ED_CHECKOUT_DT_TM) %>%
    rename(mrn=MRN, name=NAME, dob=DOB,  
    arrival_time=ED_CHECKIN_DT_TM, departure_time=ED_CHECKOUT_DT_TM) %>%
    distinct(mrn, name, dob, arrival_time, .keep_all = TRUE) 

print("Loaded second linkage table")

linkage_table <- linkage_table %>% 
    inner_join(departure_times, by=c("mrn", "name", "dob", "arrival_time")) %>%
    filter(csn %in% visits$csn) #Keep only CSNs in final dataset.

print("Constructed linkage table")

#Add a column to the visits table indicating whether the visit takes place before the
#EPIC transition.
epic <- linkage_table %>% mutate(pre_epic=ifelse(
    year(ymd_hm(arrival_time)) < 2024 | (year(ymd_hm(arrival_time)) == 2024 & month(ymd_hm(arrival_time)) <= 5),
    1, 0
)) %>% select(csn, pre_epic)

print("Identified pre-EPIC visits")

#Now filter out all post-EPIC visits:
visits <- visits %>% inner_join(epic, by="csn") %>%
    filter(pre_epic==1) %>% select(-c(pre_epic, V1))


write.csv(visits, paste0(savepath, "prediction/preprocessed-visits-pre-epic.csv"))


#Now load in medications-- they have CSNs, med_date_time, medication and route.
#Medications can be trusted throughout 2019-2024.
medications <- as.data.frame(fread(paste0(loadpath, "LaCava_Meds_May1_V1.csv"))) %>%
    select(CSN, MED_DATE_TIME, MEDICATION, ROUTE) %>% rename(csn=CSN, event_time=MED_DATE_TIME) %>%
    inner_join(linkage_table, by="csn") %>%
    mutate(timestamp = as.numeric(difftime(ymd_hm(event_time), ymd_hm(arrival_time), units="mins"))) %>%
    filter(timestamp >= 0) %>%
    filter(as.numeric(difftime(ymd_hm(departure_time), ymd_hm(event_time), units="mins")) >= 0) #Keep only meds dispensed during the ED visit.

#Clean up medication routes: eye, ear, oral, tube, intramuscular, subcutaneous, inhalation, nasal, topical, intravenous, rectal
medications$ROUTE <- case_when(
    toupper(medications$ROUTE) %in% c("BOTH EYES", "EYE BOTH", "EYE LEFT",
         "EYE RIGHT", "INTRAOCULAR", "LEFT EYE", "OPHTHALMIC (EYE)", "OPTH", "RIGHT EYE") ~ "by_eye",
    toupper(medications$ROUTE) %in% c("EACH EAR", "EAR BOTH", "EAR LEFT",
        "EAR RIGHT", "OTIC", "OTIC (EAR)", "RIGHT EAR") ~ "by_ear",
    toupper(medications$ROUTE) %in% c("BUCCAL", "CHEWED", "ORAL", "PO",
        "SUBLINGUAL", "SL", "SWISH & SPIT", "SWISH & SWALLOW", "SWISH AND SPIT",
        "SWISH AND SWALLOW", "SWISH OR SWAB", "TO ORAL MUCOSA", "MOUTH/THROAT") ~ "oral",
    toupper(medications$ROUTE) %in% c("G-TUBE", "GTUBE", "NASOGASTRIC TUBE",
        "NASOJEJUNAL TUBE", "J-TUBE", "JTUBE", "NJ", "NG") ~ "enteral",
    toupper(medications$ROUTE) %in% c("IM", "INTRAMUSCULAR") ~ "intramuscular",
    toupper(medications$ROUTE) %in% c("ID", "INTRADERMAL") ~ "intradermal",
    toupper(medications$ROUTE) %in% c("SUBCUTANEOUS") ~ "subcutaneous",
    toupper(medications$ROUTE) %in% c("INH", "INHALATION", "NEB", "NEBULIZATION", "MDI") ~ "inhalation",
    toupper(medications$ROUTE) %in% c("EACH NOSTRIL", "LEFT NOSTRIL", "NASAL") ~ "nasal",
    toupper(medications$ROUTE) %in% c("ICU-IV", "IV", "INTRAVENOUS", "INTRAVENOUS PUSH", 
        "IV LOADING DOSE", "IV LOCK", "IV PUSH", "IV EXTEND", "INTRA-CATHETER") ~ "IV",
    toupper(medications$ROUTE) %in% c("TRANSDERMAL", "TOPICAL", "TD") ~ "topical",
    toupper(medications$ROUTE) %in% c("PR", "RECTAL") ~ "rectal",
    .default =  "other"
    
)

#Now get the top 50 med-route combinations:
medications <- medications %>% mutate(medroute=paste0(toupper(MEDICATION), "_", toupper(ROUTE)))
common_medroutes <- medications %>% group_by(medroute) %>%
    summarise(count=n()) %>% arrange(desc(count))

print(paste("Med-route combinations:", nrow(common_medroutes)))
print(paste("Number of orders of 50th most common medroute", common_medroutes$count[50]))
print(paste("Number of orders of 100th most common medroute", common_medroutes$count[100]))
print(paste("Number of orders of 200th most common medroute", common_medroutes$count[200]))
print(paste("Number of orders of 500th most common medroute", common_medroutes$count[500]))
print(paste("Number of orders of 1000th most common medroute", common_medroutes$count[1000]))

medications <- medications %>% filter(medroute %in% common_medroutes$medroute[1:50]) %>%
    select(csn, timestamp, medroute) 

write.csv(medications, paste0(savepath, "prediction/medications.csv"))

print("Linked medications")



# # Now load in labs: in the 'new file' we have CSN, LAB_DATE_TIME, LAB, and LAB VALUE.
# # Cross-comparison with the old dataframe suggests that LAB_DATE_TIME is the time at which the lab RESULT
# # is recorded, not the time at which the lab is ordered.
# # Labs can also be trusted during the whole timeframe of the study.

labs <- as.data.frame(fread(paste0(loadpath, "LaCava_Complete_Lab_Records.csv"))) %>%
    select(CSN, LAB_DATE_TIME, LAB, LAB_VALUE) %>% rename(csn=CSN, event_time=LAB_DATE_TIME) %>%
    inner_join(linkage_table, by="csn") %>%
    mutate(timestamp = as.numeric(difftime(ymd_hm(event_time), ymd_hm(arrival_time), units="mins"))) %>%
    filter(timestamp >= 0) %>%
    filter(as.numeric(difftime(ymd_hm(departure_time), ymd_hm(event_time), units="mins")) >= 0)  %>%
    mutate(numeric_conversion = as.numeric(LAB_VALUE)) %>%
    mutate(lab_value = case_when(
        is.na(numeric_conversion) & (toupper(LAB_VALUE) %in% c("NEGATIVE", "NOT DETECTED", "NORMAL", "0.2 (NORMAL)", "1 (NORMAL)", "NONE", "NO", "NON-REACTIVE", "NONREACTIVE") | 
            startsWith(toupper(LAB_VALUE), "<")) ~ "negative",
        is.na(numeric_conversion) & (toupper(LAB_VALUE) %in% c("TRACE", "YES", "HIGH", "REACTIVE", "PRESENT", "POSITIVE (ANY DEGREE OF UNIFORM PINK COLOR)", "POSITIVE", "DETECTED", "ABNORMAL") | 
            startsWith(toupper(LAB_VALUE), ">") | grepl("+", toupper(LAB_VALUE))) ~ "positive",
        .default = LAB_VALUE
    )) %>% filter(!is.na(numeric_conversion) | lab_value %in% c("positive", "negative")) #Handle common non-numeric values




# #Handle non-numeric lab values:
# non_numeric <- labs %>% mutate(non_numeric=as.numeric(LAB_VALUE)) %>%
#     filter(is.na(non_numeric)) %>% mutate(lab_value=toupper(LAB_VALUE)) %>% 
#     group_by(lab_value) %>%
#     summarise(count=n()) %>% arrange(desc(count))

# print(non_numeric)
# write.csv(non_numeric, paste0(savepath, "prediction/non-numeric-lab-values.csv"))



#Now, again, keep only the top 50 lab results.
common_labs <- labs %>% mutate(lab=toupper(LAB)) %>% group_by(lab) %>% 
    summarise(count=n()) %>% arrange(desc(count))

print(paste("Labs:", nrow(common_labs)))
print(paste("Number of orders of 50th most common lab", common_labs$count[50]))
print(paste("Number of orders of 100th most common lab", common_labs$count[100]))
print(paste("Number of orders of 200th most common lab", common_labs$count[200]))
print(paste("Number of orders of 500th most common lab", common_labs$count[500]))
print(paste("Number of orders of 1000th most common lab", common_labs$count[1000]))

#Define age groups for linkage
age_group_linkage <- visits %>% select(csn, age_group)

#Save records of 50 most common labs
#Now standardise lab values within each age group according to the max and minimum

labs <- labs %>% mutate(lab=toupper(LAB)) %>% select(-LAB) %>% filter(lab %in% common_labs$lab[1:50]) %>% 
    select(csn, timestamp, lab, lab_value) %>%
    inner_join(age_group_linkage, by="csn") %>%
    group_by(lab) %>% mutate(classification=ifelse(all(lab_value %in% 
        c("negative", "positive")), "binary", "cts")) %>% ungroup()

#Handle binary and continuous labs separately.
binary_labs <- labs %>% filter(classification=="binary") %>%
    mutate(lab_result=ifelse(lab_value=="positive", 1, 0)) %>% #Convert binary lab results to 1s and 0s
    select(csn, timestamp, lab, lab_result)

#Scale cts labs to be between 0 and 1.
cts_labs <- labs %>% filter(classification=="cts") %>%
    group_by(lab, age_group) %>% mutate(
        lower_bound = quantile(as.numeric(lab_value), na.rm=TRUE, probs=c(0.01)), #Assign 1st and 99th percentiles as lower and upper bounds.
        upper_bound = quantile(as.numeric(lab_value), na.rm=TRUE, probs=c(0.99))
    ) %>% ungroup() %>%
    mutate(rounded_value=case_when(
        as.numeric(lab_value) < lower_bound | lab_value == "negative" ~ lower_bound,
        as.numeric(lab_value) > upper_bound | lab_value == "positive" ~ upper_bound,  #Any verbal 'negative' or 'positives' are set to lower and upper bounds.      
        .default = as.numeric(lab_value)
    ),
        lab_result = (rounded_value-lower_bound)/(upper_bound-lower_bound)) %>% #Scale everything to be between 0 and 1.
    select(csn, timestamp, lab, lab_result)

labs <- rbind(binary_labs, cts_labs)

write.csv(labs, paste0(savepath, "prediction/labs.csv"))

print("Converted labs")

#Now load radiological tests: we want test, MRN,  ED_ARRIVAL_TIME (to link CSN from the dispositions file),
#RADIOLOGY_TEST for time of tests. NOTE: NO TESTS ARE RECORDED AFTER MAY 2024; we handle this using the pre-epic flag.
#Type issues mean we struggle to link on name and DOB, but MRN and arrival time should be sufficient.
rads <- as.data.frame(fread(paste0(loadpath, "LaCava_Radiology_Apr14.csv"))) %>%
    select(MRN,ED_ARRIVAL_TIME, EVENT_END_DT_TM, RADIOLOGY_TEST) %>%
    rename(mrn=MRN,  arrival_time=ED_ARRIVAL_TIME,
        event_time=EVENT_END_DT_TM, test=RADIOLOGY_TEST) 

#Convert to character form
rads$mrn <- as.character(rads$mrn)

rads <- rads %>%
    inner_join(linkage_table, by=c("mrn", "arrival_time")) %>%
    mutate(timestamp = as.numeric(difftime(ymd_hm(event_time), ymd_hm(arrival_time), units="mins"))) %>%
    filter(timestamp >= 0) %>%
    filter(as.numeric(difftime(ymd_hm(departure_time), ymd_hm(event_time), units="mins")) >= 0) %>% #Keep tests that happen during the visit.
    select(csn, timestamp, test)


write.csv(rads, paste0(savepath, "prediction/tests.csv"))

print("Linked tests")

#Link vitals type by type:
vitals <- as.data.frame(fread(paste0(savepath, "intermediate-files/vitals-with-age-and-pain.csv")))

vitals_to_take <- c("hr", "rr", "sbp", "pain", "sp_o2")
vitals_to_normalise <- c("hr", "rr") 

#Filter the values we want and convert them to numeric form.
vitals_across_stay <- vitals %>% 
    filter(measure %in% vitals_to_take) 
vitals_across_stay$value <- as.numeric(vitals_across_stay$value)



#Transform BP reading to 'distance from PALS criteria'.
sbp_readings <- which(vitals_across_stay$measure=="sbp")
sbp <- vitals_across_stay$value[sbp_readings]
days <- vitals_across_stay$age_in_days[sbp_readings]

vitals_across_stay$value[sbp_readings] <- ifelse(days <= 28 & sbp < 60, pmax(60-sbp, 0) , ifelse(
    days > 28 & 365 & sbp < 70, pmax(70-sbp, 0), ifelse(
    days > 365 & days < 365.25*10 & sbp < (70 + 2*days/(365.25)), pmax((70 + 2*days/(365.25))-sbp, 0), ifelse( 
    days >= 365.25*10 & sbp < 90, pmax(90-sbp, 0), 0))))

# #Transform temperature reading to 'positive distance from 38C'.
# temperature_readings <- which(vitals_across_stay$measure=="temp")
# vitals_across_stay$value[temperature_readings] <- ifelse(vitals_across_stay$value[temperature_readings] > 38, vitals_across_stay$value[temperature_readings]-38, 0)

#Filter out non-numeric values.
vitals_across_stay <- vitals_across_stay %>% filter(!is.na(value))

write.csv(vitals_across_stay, paste0(savepath, "prediction/processed-vitals.csv"))
    
print("Processed vitals")


#Load visits and their duration.
visits <- as.data.frame(fread(paste0(savepath, "prediction/preprocessed-visits-pre-epic.csv")))
# print(head(visits))

#Define a function to sample at the start and end of visits, and every hour during.
get_snapshots <- function(i, visits) {
    num_snapshots <- floor(visits$ed_los[i]/60)
    if(num_snapshots==0) {
        snapshots <- c(0, visits$ed_los[i])
    } else {
        snapshots <- c(0, (1:num_snapshots)*60, visits$ed_los[i])
    }
    return(list(csn=rep(visits$csn[i], length(snapshots)), snapshot=snapshots))
}

#Get snapshot times for all visits.
samples <-  as.data.frame(data.table::rbindlist(purrr::map(1:nrow(visits), ~get_snapshots(.x, visits), .progress=TRUE))) 

print(nrow(samples))
samples <- samples %>%
        distinct(csn, snapshot)

print(nrow(samples))
write.csv(samples, paste0(savepath, "prediction/sample-times.csv"))

#Load sample times for each visit.
samples <- as.data.frame(fread(paste0(savepath, "prediction/sample-times.csv")))
#samples$csn <- as.double(samples$csn)


#Link medications.
medications <- as.data.frame(fread(paste0(savepath, "prediction/medications.csv"))) %>% 
    inner_join(samples, by="csn", relationship="many-to-many") %>%
    filter(timestamp <= snapshot) %>% #Medication dispensed before the snapshot
    group_by(csn, medroute, snapshot) %>% summarise(dispensed=1) %>% ungroup() %>%
    tidyr::pivot_wider(id_cols=c(csn, snapshot), names_from=medroute, names_prefix="medroute_", 
        values_from=dispensed, values_fill = 0) %>% clean_names() #Pivot to a dataframe with 1 if a medroute has been dispensed before this point and 0 otherwise

#Link tests.
tests <- as.data.frame(fread(paste0(savepath, "prediction/tests.csv"))) %>% 
    inner_join(samples, by="csn", relationship="many-to-many") %>%
    filter(timestamp <= snapshot) %>% #Test ordered before the snapshot
    group_by(csn, test, snapshot) %>% summarise(ordered=1) %>% ungroup() %>%
    tidyr::pivot_wider(id_cols=c(csn, snapshot), names_from=test, names_prefix="test_", 
        values_from=ordered, values_fill = 0) %>% clean_names() #Pivot to a dataframe with 1 if a test has been ordered before this point and 0 otherwise

#Link labs. Keep the most recent lab result available of any type; set values to NA if the test was never ordered.
labs <-  as.data.frame(fread(paste0(savepath, "prediction/labs.csv"))) %>%
    inner_join(samples, by="csn", relationship="many-to-many") %>%
    filter(timestamp <= snapshot) %>% #Test ordered before the snapshot
    group_by(csn, lab, snapshot) %>% arrange(timestamp) %>% 
    summarise(most_recent_result=last(lab_result)) %>% ungroup() %>%
    tidyr::pivot_wider(id_cols=c(csn, snapshot), names_from=lab, names_prefix="lab_", 
        values_from=most_recent_result, values_fill = NA) %>% clean_names() #Results are 0-1 if a result has been retrieved, NA otherwise.

#Link vitals (get the mean, min, max of each visit up to each snapshot, and then normalise
#taking the LAST min/mean/max (i.e. the overall) as the focal measurements for normalisation.

vitals <- as.data.frame(fread(paste0(savepath, "prediction/processed-vitals.csv"))) %>% 
    select(-V1) %>%
    rename(event_time=time_since_arrival) %>% inner_join(samples, by="csn", relationship="many-to-many") %>%
    filter(event_time <= snapshot) %>% select(csn, snapshot, measure, value, event_time)

print(head(vitals))
vitals_dt <- as.data.table(vitals)

#Get mean, min, max at each point in the visit.
vitals_dt <- vitals_dt[, .(
    mean = mean(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE)
), by = .(csn, snapshot, measure)]

#Pivot:
vitals <- as.data.frame(dcast(vitals_dt, csn + snapshot ~ measure, 
                value.var = c("mean", "max", "min"), 
                sep = "_")) %>% inner_join(select(visits, c(csn, age_group)), by="csn")
    

#Normalise these measurements by age group.
#First, mark the final snapshot of each visit:
final_snapshot_vitals <- vitals %>% group_by(csn) %>% 
    mutate(final=(snapshot==max(snapshot))) %>% ungroup() %>%
    filter(final) %>% select(-final)

#Perform normalisation.
for (group in unique(vitals$age_group)) {
    #Find indices of visits in both dataframes corresponding to the age group
    indices_in_vitals <- which(vitals$age_group==group)
    indices_in_final_snapshot_vitals <- which(final_snapshot_vitals$age_group==group)

    for (col in colnames(vitals)[
            endsWith(colnames(vitals), "hr") |
            endsWith(colnames(vitals), "rr")]) {
        
  
        mu <- mean(final_snapshot_vitals[indices_in_final_snapshot_vitals,][[col]], na.rm=TRUE)
        sigma <- sd(final_snapshot_vitals[indices_in_final_snapshot_vitals,][[col]], na.rm=TRUE)
        
        vitals[indices_in_vitals, col] <- (vitals[indices_in_vitals, col]-mu)/sigma #Since XGBoost can handle missing values, we don't take the absolute value after normalising as we don't need to 'reserve' -1.
    }
}

#Drop unnecessary age group column. NOTE: at the end of this, pain is continuous.
vitals <- vitals %>% select(-age_group)

print("check 1")

#Now attach other information!-- Baseline factors, assumed known during or shortly after triage.

#Convert categorical variables.
factors_to_convert <- c("sex", "race", "ed_arrival_mode", "language", "triage_acuity", 
    "age_group", "insurance", "state_of_origin",
    "year_of_arrival", "season", "time_of_day") #pain-unknown is already marked, so default here is none
refs <- c("M", "Non-Hispanic White", "Walk in", "English", "unknown",
            "fifteen_and_older", "Private",
        "in-state", 2019, "winter", "morning")

for (i in 1:length(factors_to_convert)) {
    var <- factors_to_convert[i]
    ref <- refs[i]

    #convert NA to unknown in the specific case of triage acuity
    visits[[var]] <- ifelse(is.na(visits[[var]]), "unknown", visits[[var]])


    assert(ref %in% unique(visits[[var]]))
    for (val in unique(visits[[var]])) {
        if (!is.na(val) & !(val==ref)) {
            name <- paste0(var, "_", val)
            visits[[name]] <- ifelse(visits[[var]]==val, 1, 0)
            assert(!any(is.na(visits[[name]])))
        }
    }
}
print('check 2')

#Divide into test and training cohorts
month_and_year <- linkage_table %>% 
    select(csn, arrival_time) %>%
    mutate(month_of_arrival=month(arrival_time), year_of_arrival=year(arrival_time)) %>%
    select(-arrival_time)

#Delete converted factors
visits <- visits %>% select(-all_of(factors_to_convert)) %>%
    inner_join(month_and_year, by="csn") %>%
    mutate(train_or_test=ifelse((year_of_arrival < 2023) | (year_of_arrival == 2023 & month_of_arrival < 6), "train", "test")) %>%  #Take everything up to May 2023 as training data and after that as testing data.
    select(-c(month_of_arrival, year_of_arrival)) %>% clean_names()

#Attach weight: for these purposes we want the number of standard deviations from the mean, but want to allow this,
#like HR and RR, to be positive or negative, since XGBoost can handle nonlinearity and NAs
visits <- visits %>% select(-all_of(starts_with("weight")))
weight <- as.data.frame(fread(paste0(loadpath, "LaCava_Demopgraphics_May1_V1.csv"))) %>%
    select(CSN, WEIGHT_KG) %>% rename(csn=CSN, weight=WEIGHT_KG)
visits <- visits %>% left_join(weight, by="csn")

#Now normalise weight:
for (group in unique(visits$age_group)) {
    #Find indices of visits corresponding to the age group
    idx <- which(visits$age_group==group)
    mu <- mean(visits$weight[idx], na.rm=TRUE)
    sigma <- sd(visits$weight[idx], na.rm=TRUE)
    visits$weight[idx] <- (visits$weight[idx]-mu)/sigma
}

print("check 3")

#Now select baseline factors (vars-known-at-triage in bch-produce-odds-ratios-- modified so that triage acuity is included but triage vitals are NOT)
all_vars <- colnames(visits)

vars_known_at_triage <- c("sex_f", all_vars[startsWith(all_vars, "complaint_contains_")],
     all_vars[startsWith(all_vars, "triage_acuity_")], 
     "num_previous_admissions", "num_previous_visits_without_admission", 
     all_vars[startsWith(all_vars, "pre_diagnosis")], all_vars[startsWith(all_vars, "language")], 
     all_vars[startsWith(all_vars, "race_")], all_vars[startsWith(all_vars, "ed_arrival_mode_")],
     all_vars[startsWith(all_vars, "age_group_")], all_vars[startsWith(all_vars, "insurance_")],
     all_vars[startsWith(all_vars, "year_of_arrival")], all_vars[startsWith(all_vars, "season_")],
     all_vars[startsWith(all_vars, "time_of_day_")], "is_weekend", "miles_travelled", "sdi_score", "weight",
     all_vars[startsWith(all_vars, "state_of_origin")],
     "is_trans_or_nb", "pseudo_nedocs")

#Now take these factors and join them up with snapshots:
baseline_factors <- visits %>% select(all_of(c("csn", "train_or_test", "is_admitted", vars_known_at_triage))) 
visit_snapshots <- samples %>% 
    inner_join(baseline_factors, by="csn") %>%
    left_join(medications, by=c("csn", "snapshot")) %>%
    left_join(labs, by=c("csn", "snapshot")) %>%
    left_join(tests, by=c("csn", "snapshot")) %>%
    left_join(vitals, by=c("csn", "snapshot")) 

print("check 4")

#Separate test and training data
training_data <- visit_snapshots %>% filter(train_or_test=="train") %>% select(-train_or_test)
test_data <- visit_snapshots %>% filter(train_or_test=="test") %>% select(-train_or_test)

#Save both test and training data.
write.csv(training_data, paste0(savepath, "prediction/training-data.csv"))
write.csv(test_data, paste0(savepath, "prediction/test-data.csv"))

training_data <- as.data.frame(fread( paste0(savepath, "prediction/training-data.csv")))
write.csv(training_data, paste0(savepath, "prediction/to-transfer/training-data.csv"))


test_data <- as.data.frame(fread( paste0(savepath, "prediction/test-data.csv")))
write.csv(test_data, paste0(savepath, "prediction/to-transfer/test-data.csv"))

#Select 1 sample from each visit for hyperparameter tuning
training_data <- as.data.frame(fread( paste0(savepath, "prediction/training-data.csv"))) %>%
    select(-V1) %>%
    group_by(csn) %>% mutate(chosen_snapshot=sample(snapshot, 1)) %>% ungroup() %>%
    filter(snapshot==chosen_snapshot)

write.csv(training_data, paste0(savepath, "prediction/to-transfer/training-data-for-cv.csv"))