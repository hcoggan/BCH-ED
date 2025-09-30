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


#Attach records of ECGs, medications, labs and tests to BCH visits (right now: did they receive any labs, any tests, any meds/IV meds?)

#First load visits.
visits <- as.data.frame(fread(paste0(savepath, "preprocessed-visits.csv")))

print(paste("Before linking events, we have", nrow(visits), "visits."))

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

#Now load in medications-- they have CSNs, med_date_time, medication and route.
#Medications can be trusted throughout 2019-2024.
medications <- as.data.frame(fread(paste0(loadpath, "LaCava_Meds_May1_V1.csv"))) %>%
    select(CSN, MED_DATE_TIME, MEDICATION, ROUTE) %>% rename(csn=CSN, event_time=MED_DATE_TIME) %>%
    inner_join(linkage_table, by="csn") %>%
    filter(as.numeric(difftime(ymd_hm(event_time), ymd_hm(arrival_time), units="mins")) >= 0) %>%
    filter(as.numeric(difftime(ymd_hm(departure_time), ymd_hm(event_time), units="mins")) >= 0) %>% #Keep only meds dispensed during the ED visit.
    group_by(csn) %>% summarise(received_any_medications=1, received_any_iv=any(grepl("IV", toupper(ROUTE)) | grepl("INTRAVENOUS", toupper(ROUTE))))

print("Linked medications")

#Now load in labs: in the 'new file' we have CSN, LAB_DATE_TIME, LAB, and LAB VALUE.
#Cross-comparison with the old dataframe suggests that LAB_DATE_TIME is the time at which the lab RESULT
#is recorded, not the time at which the lab is ordered.
#Labs can also be trusted during the whole timeframe of the study.

labs <- as.data.frame(fread(paste0(loadpath, "LaCava_Complete_Lab_Records.csv"))) %>%
    select(CSN, LAB_DATE_TIME) %>% rename(csn=CSN, event_time=LAB_DATE_TIME) %>%
    inner_join(linkage_table, by="csn") %>%
    filter(as.numeric(difftime(ymd_hm(event_time), ymd_hm(arrival_time), units="mins")) >= 0) %>%
    filter(as.numeric(difftime(ymd_hm(departure_time), ymd_hm(event_time), units="mins")) >= 0) %>% #Keep only labs whose results were received during the ED visit.
    group_by(csn) %>% summarise(received_any_labs=1)



print("Linked labs")

#Now load radiological tests: we want test, MRN,  ED_ARRIVAL_TIME (to link CSN from the dispositions file),
#RADIOLOGY_TEST for time of tests. NOTE: NO TESTS ARE RECORDED AFTER MAY 2024; we handle this using the pre-epic flag.
#Type issues mean we struggle to link on name and DOB, but MRN and arrival time should be sufficient.
rads <- as.data.frame(fread(paste0(loadpath, "LaCava_Radiology_Apr14.csv"))) %>%
    select(MRN,ED_ARRIVAL_TIME, EVENT_END_DT_TM) %>%
    rename(mrn=MRN,  arrival_time=ED_ARRIVAL_TIME,
        event_time=EVENT_END_DT_TM) 

#Convert to character form
rads$mrn <- as.character(rads$mrn)

rads <- rads %>%
    inner_join(linkage_table, by=c("mrn", "arrival_time")) %>%
    filter(as.numeric(difftime(ymd_hm(event_time), ymd_hm(arrival_time), units="mins")) >= 0) %>%
    filter(as.numeric(difftime(ymd_hm(departure_time), ymd_hm(event_time), units="mins")) >= 0) %>% #Keep only labs whose results were received during the ED visit.
    group_by(csn) %>% summarise(received_any_tests=1)


print("Linked tests")

#Link ECGs (not reliable after EPIC transition; files trusted up to and including May 2024.)
#Again, link on MRN and arrival time.


# ecgs <- as.data.frame(fread(paste0(loadpath, "LaCava_ECG_Sep4.csv"))) %>%
#     select(MRN, ED_ARRIVAL_TIME, EVENT_END_DT_TM) %>% 
#     rename(mrn=MRN, arrival_time=ED_ARRIVAL_TIME, event_time=EVENT_END_DT_TM) 

# #Convert to amenable types.
# ecgs$mrn <- as.character(ecgs$mrn)


# ecgs <- ecgs %>%
#     inner_join(linkage_table, by=c("mrn", "arrival_time")) %>%
#     filter(as.numeric(difftime(ymd_hm(event_time), ymd_hm(arrival_time), units="mins")) >= 0) %>%
#     filter(as.numeric(difftime(ymd_hm(departure_time), ymd_hm(event_time), units="mins")) >= 0) %>% #Keep only labs whose results were received during the ED visit.
#     group_by(csn) %>% summarise(received_any_ecgs=1)

# print("Linked ECGs")

#Now link all these files to the visit table.
visits <- visits %>%
    inner_join(epic, by="csn") %>% #Variable describing whether a visit happened before the EPIC transition.
    left_join(medications, by="csn") %>%
    left_join(labs, by="csn") %>%
    left_join(rads, by="csn") 


print("Link visits to all dataframes")

#Change any NAs to 0s (as it means no tests, etc were recorded).
for (col in colnames(visits)[startsWith(colnames(visits), "received_any")]) {
    visits[[col]] <- ifelse(is.na(visits[[col]]), 0, visits[[col]])
}


print(paste("After linking events, we have", nrow(visits), "visits."))

#Save this file.
write.csv(visits, paste0(savepath, "preprocessed-visits-with-linked-events.csv"))

