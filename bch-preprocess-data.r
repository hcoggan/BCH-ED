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




# #First, load dispositions; this becomes the main dataframe.
# visits <- as.data.frame(fread((paste0(loadpath, "LaCava_Dispositions_Aug7.csv")))) %>%
#     select(MRN, NAME, DOB, ZIP_CODE, ETHNICITY,
#     ADMIT_REQ_DT_TM, ROOMING_TIME, ED_CHECKIN_DT_TM, ED_DISPOSITION, 
#     ED_CHECKOUT_DT_TM) %>% rename(ethnicity=ETHNICITY, zipcode=ZIP_CODE, mrn=MRN, name=NAME, dob=DOB,  
#     arrival_time=ED_CHECKIN_DT_TM, departure_time=ED_CHECKOUT_DT_TM, 
#     rooming_time=ROOMING_TIME, admission_request_time=ADMIT_REQ_DT_TM, disposition=ED_DISPOSITION) %>%
#     distinct(mrn, name, dob, arrival_time, .keep_all = TRUE) #Rows are duplicated whenever a patient moves between department; we drop all relevant columns and keep only
#     #the rows which uniquely identify a VISIT (combination of MRN, name, DOB, and arrival time.)

# print(paste("Initially, we have", nrow(visits), "visits from", length(unique(visits$mrn)), "unique patients."))

# #To link other variables, like CSN, race, sex, insurance, weight and primary language, we need to pull from an older dataframe.
# #Here we assume that ED_ARRIVAL_TIME and ED_CHECKIN_TIME are the same.
# #We don't lose many visits doing this, so it should be fine.
# demographics <- as.data.frame(fread(paste0(loadpath, "LaCava_Demopgraphics_May1_V1.csv"))) %>%
#     select(CSN, MRN, NAME, DOB, RACE, SEX, GENDER, WEIGHT_KG, ED_ARRIVAL_MODE, ED_ARRIVAL_TIME, PREFERRED_LANGAUAGE,
#     PRIMARY_INSURANCE_PAYORXX, ED_COMPLAINT) %>% rename(csn=CSN, mrn=MRN, name=NAME, dob=DOB, race=RACE, sex=SEX,
#     weight=WEIGHT_KG, ed_arrival_mode=ED_ARRIVAL_MODE, arrival_time=ED_ARRIVAL_TIME, language=PREFERRED_LANGAUAGE,
#      insurance=PRIMARY_INSURANCE_PAYORXX, gender=GENDER, complaint=ED_COMPLAINT) %>%
#     group_by(csn) %>% mutate(duplicated=n()>1) %>% filter(!duplicated) %>% select(-duplicated) #CSNs may be duplicated across rows; drop all duplicated CSNs, 
#     #as they are non-unique identifiers and may correspond to several visits, which then can't be uniquely identified.



# #Link by name, DOB, and arrival time (assuming that arrival_time and checkin_time are the same)
# visits <- visits %>% inner_join(demographics, by=c("mrn", "name", "dob", "arrival_time")) %>% 
#     mutate(age_in_days=as.numeric(difftime(ymd_hm(arrival_time), ymd(dob), units="days"))) %>%
#     select(-c(name, dob))


# print(head(visits))

# print(paste("After linking demographics, we have", nrow(visits), "visits from", length(unique(visits$mrn)), "unique patients."))


# #Save this linked file.
# write.csv(visits, paste0(savepath, "intermediate-files/visits-with-linked-demographics.csv"))

# #Load this file back.
# visits <- as.data.frame(fread((paste0(savepath, "intermediate-files/visits-with-linked-demographics.csv"))))


# #First, handle dispositions. Check that everyone who isn't admitted has an unavailable admission request time:
# admit_disposition_times <- visits %>% 
#     group_by(disposition) %>% summarise(
#         total=n(), 
#         blank_req_time=sum(admission_request_time==""),
#         frac_without_request_time=blank_req_time/total) %>%
#     select(disposition, total, frac_without_request_time) %>% arrange(desc(total))

# write.csv(admit_disposition_times, paste0(savepath, "intermediate-files/fraction-of-visits-with-linked-admission-request-time.csv"))

# #98% of visits which have a 'home' disposition have no admission request time; for 'discharge', it's 98.8%; for LWBS it's also close to 100%.
# #0.6% of visits which have an 'ED Patient Admitted' disposition have a request time; for Admit, it's 0.1%.

# #29,721 visits have an Unknown disposition, and of these, 51% have an admission request time.
# #It's fairly safe to assume that those without an admission request time were not admitted, and those with one were admitted.

# #Load vitals (except temperature.)
# vitals <- as.data.frame(fread(paste0(loadpath, "LaCava_Vitals_Aug7.csv"))) %>% 
#     select(CSN, VITAL_DATE_TIME, MEASUREMENT_DESC, VITAL) %>%
#     rename(csn=CSN, vital_time=VITAL_DATE_TIME,
#         measure=MEASUREMENT_DESC, value=VITAL) %>% filter(!(measure=="TEMPERATURE"))

# # #Load temperatures.
# # temps <- as.data.frame(fread(paste0(loadpath, "LaCava_Temperature.csv"))) %>% 
# #     select(CSN, VITAL_DATE_TIME, MEASUREMENT_DESC, VITAL) %>%
# #     rename(csn=CSN, vital_time=VITAL_DATE_TIME,
# #         measure=MEASUREMENT_DESC, value=VITAL)

# # print(table(month(temps$vital_time), year(temps$vital_time)))

# # #Combine vitals.
# # vitals <- rbind(vitals, temps) 

# #Save and load combined vitals.
# write.csv(vitals, paste0(savepath, "intermediate-files/combined-vitals.csv"))
# vitals <- as.data.frame(fread(paste0(savepath, "intermediate-files/combined-vitals.csv")))

# #Mark visits with no vitals (including vitals we don't use, like MAP and DBP).
# visits <- visits %>% mutate(any_recorded_vitals=(csn %in% vitals$csn))

# #Categorise dispositions.
# visits$disposition <- case_when(
#     (visits$disposition %in% c("Home", "Discharge")) | 
#     (visits$disposition == "Unknown" & visits$admission_request_time=="" & visits$any_recorded_vitals)  ~ "Discharge", # Assume a visit with recorded vitals and no admission request time is a discharge.
#     (visits$disposition %in% c("ED Patient Admitted", "Admit", "Send to OR")) | 
#     (visits$disposition == "Unknown" & !(visits$admission_request_time=="")) ~ "Admit", #We assume a visit with a recorded admit request time is an admission.
#     .default = "Other"
# )

# #Filter out those not admitted or discharged.
# visits <- visits %>% filter(disposition %in% c("Admit", "Discharge")) 

# print(paste("After filtering those not admitted or discharged, we have", nrow(visits), "visits from", length(unique(visits$mrn)), "unique patients."))
# write.csv(visits, paste0(savepath, "intermediate-files/visits-with-processed-dispositions.csv"))

# #Load visits with preprocessed dispositions.
# visits <- as.data.frame(fread((paste0(savepath, "intermediate-files/visits-with-processed-dispositions.csv")))) %>%
#     select(-c(V1, any_recorded_vitals)) %>% mutate(is_admitted=ifelse(disposition=="Admit", 1, 0))
# print(paste("After filtering those not admitted or discharged, we have", nrow(visits), "visits from", length(unique(visits$mrn)), "unique patients."))

# #A person's sex is unknowable if it is recorded as X or U, and should be corrected to 'male'
# #if recorded as female but gender is 'transgender female' (and vice versa)
# #on the basis that 'transgender' is likelier to be correct.
# visits$sex <- case_when(
#     visits$sex %in% c("X", "U") ~ "Unknown",
#     visits$sex == "M" & visits$gender == "Transgender Male" ~ "F",
#     visits$sex == "F" & visits$gender == "Transgender Female" ~ "M",
#     .default = visits$sex
# )

# #Define whether a patient is trans or NB.
# #Only count 'NB' or 'other', not 'blank' (most visits), 'choose not to answer', 'don't know' or 'unable to collect'.
# visits$is_trans_or_nb <- ifelse(
#     (visits$sex == "F" & visits$gender %in% c("Male", "Transgender Male", "Nonbinary (e.g. genderqueer/gender nonconforming)", "Other")) |
#     (visits$sex == "M" & visits$gender %in% c("Female", "Transgender Female", "Nonbinary (e.g. genderqueer/gender nonconforming)", "Other")), 1, 0 
# )

# #Filter out patients with unknowable sex.
# visits <- visits %>% filter(sex %in% c("M", "F"))

# race_by_ethn <- visits %>% group_by(ethnicity, race) %>% summarise(num_visits=n())
# write.csv(race_by_ethn, paste0(savepath, "race-ethnicity.csv"))
# assert(1==0)

# print(paste("After filtering those without a legible sex, we have", nrow(visits), "visits from", length(unique(visits$mrn)), "unique patients."))

# #Now categorise patients by ethnicity, using ethnicity to adjust 'unknown' or 'other' race/ethnicity if the ethnicity and racial groups overlap
# #(i.e. no editorialising about what counts as 'Asian'.
# #Anyone who puts their race as Hispanic/ethnicity as Hispanic is counted as Hispanic. Anyone who puts their race as Non-Hispanic Black OR their race as other/unknown and their ethnicity as Black
# #is counted as Non-Hispanic Black.
# visits$race <- case_when(
#     visits$race == "Asian, non-Hispanic" ~ "Asian",
#     visits$race == "Black, non-Hispanic" | 
#         ((visits$race=="Unknown" | visits$race == "Another Race, non-Hispanic") 
#          & visits$ethnicity=="Black") ~ "Non-Hispanic Black", #technically 'unknown' is not 'non-Hispanic' so those 107 Unknown/Black visits could all be Hispanic, but this is unlikely; we assume that everyone who has not
#                     #declared themselves to be Hispanic is not Hispanic.
#     visits$race == "White, non-Hispanic" &  visits$ethnicity=="Hispanic or Latino" ~ "Hispanic White",           
#     visits$race == "Hispanic" |  visits$ethnicity=="Hispanic or Latino" ~ "Hispanic", 
#     visits$race == "Another Race, non-Hispanic" | visits$race == "Multiracial, non-Hispanic" ~ "Other",
#     visits$race == "White, non-Hispanic" ~ "Non-Hispanic White",
#     .default = visits$race
# )

# #Handle PRIMARY LANGUAGE, in more detail this time
# visits$language <- case_when(
#     visits$language %in% c("Arabic", "Cape Verdean", "Chinese Mandarin", "English", "Haitian Creole", "Portuguese", "Spanish") ~ visits$language, #all languages with over 1000 visits
#     .default = "Other"
# )

# #Handle INSURANCE (everything's uppercase in this field)
# visits$insurance <- case_when(
#     grepl("ACO|MEDICAID|COMMUNITY|MASSHEALTH", visits$insurance) ~ "Public",
#     .default = "Private"
# )

# #Link miles travelled, state of origin and SDI index.
# #Because the libraries are only as good as the data and don't know about newer zip codes, there are 4900 visits
# #with no miles travelled, either because their zip codes are invalid or because they're new.
# #I have tried fixing this by loading latitude and longitude points for all zip codes, from
# #the 2020 census tabulation (https://www.census.gov/geographies/reference-files/time-series/geo/gazetteer-files.html)
# #and it doesn't work (only fixes about 20 visits.)


# #Identify distance travelled and state of origin.
# #First look at state directly:

# states <- as.data.frame(reverse_zipcode(visits$zipcode)) %>%
#     distinct(zipcode, state)
# visits <- visits %>% left_join(states, by="zipcode") %>%
#     rename(home_state=state) %>%
#     mutate(
#         miles_travelled = zip_distance(zipcode, "02115")$distance,
#         state_of_origin = case_when(
#             home_state == "MA" ~ "in-state",
#             !is.na(home_state) ~ "out-of-state",
#             .default = NA
#         )
#     ) %>% select(-home_state)



# #Link SDI scores to zip codes (from the Rober Graham Center).
# sdi_scores <- read.csv(paste0(savepath, "intermediate-files/rgcsdi-2015-2019-zcta.csv")) %>% 
#     rename(sdi_score=SDI_score) %>% 
#     mutate(zipcode=ifelse(nchar(as.character(ZCTA5_FIPS))==4, paste0("0", as.character(ZCTA5_FIPS)), as.character(ZCTA5_FIPS))) %>% #Attach a 0 to the front of all 4-character zip codes.
#     select(zipcode, sdi_score)
# visits <- visits %>% left_join(sdi_scores, by="zipcode")

# print(sum(is.na(visits$sdi_score)))

# #Save this checkpoint.
# write.csv(visits, paste0(savepath, "intermediate-files/visits-with-sdi-scores.csv"))
# visits <- as.data.frame(fread(paste0(savepath, "intermediate-files/visits-with-sdi-scores.csv")))

# print(head(visits))

# #To do: ED arrival mode, visit history, ED LOS, rooming time, time to admit request, time to departure (for patients),
# #weight, complaint, triage vitals, temporal variables, diagnoses.

# #Categorise ED arrival modes.
# visits$ed_arrival_mode <- case_when(
#     visits$ed_arrival_mode %in% c("Air transport",
#         "Ambulance: Other EMS", "Ambulance", 
#         "EMS", "Police") ~ "EMS",
#     visits$ed_arrival_mode %in% c("Critical Care", 
#         "Hospital Transport", "Transfer") ~ "Transfer",
#     visits$ed_arrival_mode %in% c("Other", 
#         "Unknown") ~ "Other/Unknown",   
#     .default = "Walk in"
# )

# #Find histories of prior visits.
# #Assign each visit a timestamp: the difference in minutes between its arrival time and midnight on 1 Jan 2019; this will allow us to order visits relative to each other.
# visits <- visits %>% 
#     mutate(arbitrary_timestamp=as.numeric(difftime(lubridate::ymd_hm(arrival_time), lubridate::ymd_hm("2019-01-01 00:00"), units="mins")))


# #To calculate prior visits we have to be able to uniquely identify patients, so discard all visits with no MRN.
# visits <- visits %>% filter(!(mrn==""))


# print(paste("After filtering those without an MRN, we have", nrow(visits), "visits from", length(unique(visits$mrn)), "unique patients."))


# #Calculate prior visits within the last 30 days.
# visits <- visits %>%
#     arrange(mrn, arbitrary_timestamp) %>% group_by(mrn) %>%
#     mutate(
#         num_previous_admissions = purrr::map_dbl(row_number(), function(i) {
#             current_time <- arbitrary_timestamp[i]
#             thirty_days_ago <- current_time - 24*60*30

#             sum(arbitrary_timestamp < current_time &
#                 arbitrary_timestamp > thirty_days_ago &
#                 is_admitted==1)
#         }
#             ),
#         num_previous_visits_without_admission = purrr::map_dbl(row_number(), function(i) {
#             current_time <- arbitrary_timestamp[i]
#             thirty_days_ago <- current_time - 24*60*30

#             sum(arbitrary_timestamp < current_time &
#                 arbitrary_timestamp > thirty_days_ago &
#                 is_admitted==0
#             )}),
    
#     ) %>% ungroup() 

# #Save this checkpoint dataframe.
# write.csv(visits, paste0(savepath, "intermediate-files/visits-with-prior-visits.csv"))

# #Load this checkpoint dataframe.
# visits <- as.data.frame(fread(paste0(savepath, "intermediate-files/visits-with-prior-visits.csv"))) %>% 
#     select(-c(V1))

# #To do: 
# #complaint, diagnoses, crowdedness


# #Define ED LOS, time to admit request, time from admit decision to departure.
# visits <- visits %>% mutate(
#     ed_los = as.numeric(difftime(ymd_hm(departure_time), ymd_hm(arrival_time), units="mins")), #Total length of a patient's visit.
#     time_to_room = as.numeric(difftime(ymd_hm(rooming_time), ymd_hm(arrival_time), units="mins")), #Time until a patient is assigned a room.
#     time_to_admit_request = as.numeric(difftime(ymd_hm(admission_request_time), ymd_hm(arrival_time), units="mins")), #FOR ADMITTED PATIENTS ONLY: time until a request to admit is made.
#     time_from_request_to_admission = as.numeric(difftime(ymd_hm(departure_time), ymd_hm(admission_request_time), units="mins")), #FOR ADMITTED PATIENTS ONLY: time between the request to admit and the admission of a patient.
    
#     #Now decipher the time a patient arrived at the ED.

#     year_of_arrival = year(arrival_time),
#     month_of_arrival = month(arrival_time),
#     day_of_week_of_arrival = wday(arrival_time),
#     hour_of_arrival = hour(arrival_time),

#     #Month, weekday and hour are all coarsened into easier variables
#     season = case_when(
#         month_of_arrival %in% c(12, 1, 2) ~ "winter",
#         month_of_arrival %in% c(3, 4, 5) ~ "spring",
#         month_of_arrival %in% c(6, 7, 8) ~ "summer",
#         month_of_arrival %in% c(9, 10, 11) ~ "autumn"
#     ),
#     is_weekend = ifelse(day_of_week_of_arrival==6 | day_of_week_of_arrival==7, 1, 0),
#     time_of_day = case_when(
#         hour_of_arrival < 6 ~ "small hours", #00:00 to 05:59
#         hour_of_arrival < 12 ~ "morning", #06:00 to 11:59
#         hour_of_arrival < 18 ~ "afternoon", #12:00 to 17:59
#         hour_of_arrival < 24 ~ "evening" #18:00 to 23:59
#     )
# ) %>% select(-c(month_of_arrival, day_of_week_of_arrival, hour_of_arrival))

# #Assign age groups
# visits <- visits %>% mutate(
#     age_group = case_when(
#         age_in_days < 365.25/4 ~ "under_3_months",
#         age_in_days < 365.25/2 ~ "three_to_6_months",
#         age_in_days < 365.25 ~ "six_to_12_months",
#         age_in_days < 3*365.25/2 ~ "twelve_to_18_months",
#         age_in_days < 365.25*3 ~ "eighteen_months_to_3_years",
#         age_in_days < 365.25*5 ~ "three_to_5_years",
#         age_in_days < 365.25*10 ~ "five_to_10_years",
#         age_in_days < 365.25*15 ~ "ten_to_15_years",
#         age_in_days >= 365.25*15 ~ "fifteen_and_older",
#         .default = NA
#     )
# )

# #Save this checkpoint dataframe.
# write.csv(visits, paste0(savepath, "intermediate-files/visits-with-age-groups.csv"))
# #Load this out of memory.
# visits <- as.data.frame(fread(paste0(savepath, "intermediate-files/visits-with-age-groups.csv")))

# print(paste("After filtering out those with no MRN, we have", nrow(visits), "visits from", length(unique(visits$mrn)), "unique patients."))

# #Select arrival times and age group to identify and normalise triage vitals.
# arrival_times_for_vitals <- visits %>% select(csn, arrival_time, age_group, age_in_days, ed_los)

# #Link vitals and weight.

# #Load in pain scores.
# scores <- as.data.frame(fread(paste0(loadpath, "LaCava_Scores_May7.csv"))) %>% 
#     select(CSN, SCORE_DT_TIME, SCORE_VARIABLE, SCORE) %>% 
#     rename(csn=CSN, vital_time=SCORE_DT_TIME, measure=SCORE_VARIABLE, value=SCORE) %>%
#     filter(measure=="NRS Generalized Pain Score") %>% mutate(converted_measure="pain") %>%
#     select(-measure)


# #Then correct names of each vitals and record their timestamps within a visit.
# vitals <- as.data.frame(fread(paste0(savepath, "intermediate-files/combined-vitals.csv"))) %>%
#     filter(!(measure %in% c("MEAN ARTERIAL PRESSURE (DEVICE)", "DIASTOLIC BLOOD PRESSURE"))) %>% #Filter out MAP and DBP, which we don't need.
#     mutate(upper_measure = toupper(measure), #Rename all vitals.
#         converted_measure = case_when(
#             upper_measure == "RESPIRATORY RATE" ~ "rr",
#             upper_measure == "HEART RATE" ~ "hr",
#             upper_measure == "SYSTOLIC BLOOD PRESSURE" ~ "sbp",
#             upper_measure == "DIASTOLIC BLOOD PRESSURE" ~ "dbp",
#             upper_measure == "OXYGEN SATURATION (SPO2)" ~ "sp_o2",
#             #upper_measure %in% c("TEMP", "TEMPERATURE") ~ "temp",
#         )) %>% select(csn, vital_time, converted_measure, value)

# #Add pain scores to other vitals.
# vitals <- rbind(vitals, scores)

# print(head(vitals))

# vitals <- vitals %>%
#     rename(measure=converted_measure) %>%
#     inner_join(arrival_times_for_vitals, by="csn") %>% #Link arrival times and ages.
#     mutate(time_since_arrival=as.numeric(difftime(ymd_hm(vital_time), ymd_hm(arrival_time), units="mins"))) %>%
#     filter(time_since_arrival >= 0, time_since_arrival <= ed_los) #Keep only vitals taken during the stay.

# # Save these vitals.
# write.csv(vitals, paste0(savepath, "intermediate-files/vitals-with-age-and-pain.csv"))
    

#Load this out of memory.
visits <- as.data.frame(fread(paste0(savepath, "intermediate-files/visits-with-age-groups.csv")))
print(paste("Before linking weight and vitals, we have", nrow(visits), "visits from", length(unique(visits$mrn)), "unique patients."))



#Normalise triage and overall vitals SEPARATELY (in line with CHLA).
vitals <- as.data.frame(fread(paste0(savepath, "intermediate-files/vitals-with-age-and-pain.csv")))




vitals_to_take <- c("hr", "rr", "sbp", "pain",  "sp_o2") #EXCLUDE TEMP
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


#Isolate triage vitals.
triage_vitals <- vitals_across_stay %>% mutate(event_time=as.numeric(difftime(ymd_hm(vital_time), ymd_hm(arrival_time), units="mins"))) %>% 
    group_by(csn) %>% mutate(time_of_first_vital=min(event_time)) %>% 
    ungroup() %>% filter(event_time==time_of_first_vital)

#Convert to dataframes and summarise.
triage_vitals <- setDT(triage_vitals)
vitals_across_stay <- setDT(vitals_across_stay)


#Link both of these to visits.
triage_vitals <- triage_vitals[, 
    .(mean = mean(value, na.rm=TRUE)), 
    by = .(csn, measure)
    ][, 
    dcast(.SD, csn ~ measure, 
            value.var = c("mean"),
            sep = "_")
    ] 

triage_vitals <- data.frame(triage_vitals) %>%
        rename_with(~ paste0("triage_", .), -csn)

#Now convert overall vitals.
vitals_across_stay <- vitals_across_stay[, 
    .(max = max(value, na.rm=TRUE), 
        min = min(value, na.rm=TRUE), 
        mean = mean(value, na.rm=TRUE)),
    by = .(csn, measure)
    ][, 
    dcast(.SD, csn ~ measure, 
            value.var = c("max", "min", "mean"),
            sep = "_")
    ] 

vitals_across_stay <- data.frame(vitals_across_stay) %>% rename_with(~ paste0("overall_", .), -csn) 

#Attach these to the visit table.
visits <- visits %>% 
    left_join(triage_vitals, by="csn") %>% 
    left_join(vitals_across_stay, by="csn")

print(head(visits))

#Normalise HR and RR by age group.
cols_to_normalise <- colnames(visits)[grepl("_hr|_rr", colnames(visits))]

#Copy across these 'raw'.
for (col in cols_to_normalise) {
    visits[[paste0("raw_", col)]] <- visits[[col]]
}

for (group in unique(visits$age_group)) {
    idx <- which(visits$age_group==group)
    for (col in cols_to_normalise) {
        rel_vitals <- visits[[col]][idx]
        mean <- mean(rel_vitals, na.rm=TRUE)
        sd <- sd(rel_vitals, na.rm=TRUE)
        visits[[col]][idx] <- abs(rel_vitals-mean)/sd 

    }
}

#Mark unknown values. (At this point pain is numeric.)
all_vitals <- colnames(visits)[grepl("overall|triage", colnames(visits))]

for (col in all_vitals) {
    visits[[paste0(col, "_unknown")]] <- 0
    idx_to_replace <- which(is.na(visits[[col]]))
    visits[[paste0(col, "_unknown")]][idx_to_replace] <- 1
    if (col %in% vitals_to_normalise) {
        visits[[col]][idx_to_replace] <- 0 #because these have been normalised, so the average at each age group should be 0
    } else { #Just replace with the age-specific mean for pain, BP (TEMP EXCLUDED)
        for (group in unique(visits$age_group)) {
            idx <- which(visits$age_group==group)
            mean <- mean(visits[[col]], na.rm=TRUE)
            visits[[col]][intersect(idx, idx_to_replace)] <- mean
        }
    }
}

#Finally, categorise pain levels (this is constructed of numeric values so no default is required, and indeed this runs without error.)
pain_cols <-  colnames(visits)[grepl("pain", colnames(visits)) & !grepl("unknown", colnames(visits))] 
for (col in pain_cols) {
    visits[[col]] <- case_when(
        visits[[col]]==0 ~ "none",
        visits[[col]]<4 ~ "mild",
        visits[[col]]<7 ~ "moderate",
        visits[[col]]<=10 ~ "severe",
    )
}


#Normalise weight.
visits <- visits %>% group_by(age_group) %>%
    mutate(mu=mean(weight, na.rm=TRUE), sigma=sd(weight, na.rm=TRUE)) %>%
    ungroup() %>% mutate(converted_weight=abs(weight-mu)/sigma) %>%
    select(-weight) %>% rename(weight=converted_weight)

#Save the visits.
write.csv(visits, paste0(savepath, "intermediate-files/visits-with-weight-and-vitals.csv"))

#Load visits back in.
visits <- as.data.frame(fread(paste0(savepath, "intermediate-files/visits-with-weight-and-vitals.csv")))
print(head(visits))

#Now get the 

#Generate lists of ICD prefixes within a given range.
generate_icd_prefixes <- function(starting_letter, start, end, exceptions) {
    range <- start:end
    range <- range[!(range %in% exceptions)]
    character_range <- ifelse(range>9, as.character(range), paste0("0", as.character(range)))
    codes <- paste0(starting_letter, character_range)
    return(codes)
}

#Find the condition associated with a particular ICD code.
find_condition <- function(code, paed_icd_code_dict) {
        res <- NA
        for (condition in names(paed_icd_code_dict)) {
                for (prefix in paed_icd_code_dict[[condition]]) {
                        if (startsWith(code, prefix)) {
                                res <- condition
                                break
                        }
                }
        }
        return(res)

}


#Define a dictionary linking prefixes to conditions.
paed_icd_code_dict = list(
    "alcohol_abuse" = c("F10", "Z71.41", "K70"),
    "anemia" = c("D50", "D51", "D52", "D53",
            "D55", "D66", "D57"),
    "anxiety" = c("F06.4", "F40", "F41", "F43.0", 
            "F43.22", "F43.8", "F43.9"),
    "any_malignancy" = c("C00", "C01", "C02", "C03",
            "C04", "C05", "C06", "C07", "C08", "C09",
            "C10", "C11", "C12", "C13", "C14", "C15",
            "C16", "C17", "C18", "C19", "C20", "C21",
            "C22", "C23", "C24", "C25", "C26", "C30",
            "C31", "C32", "C33", "C34", "C35", "C36",
            "C37", "C38", "C39", "C40", "C41", "C43", 
            "C45", "C46", "C47", "C48", "C49", "C50", 
            "C51", "C52", "C53", "C54", "C55", "C56",
            "C57", "C58", "C60", "C61", "C62", "C63",
            "C64", "C65", "C66", "C67", "C68", "C69",
            "C70", "C71", "C72", "C73", "C74", "C75",
            "C76", "C81", "C82", "C84", "C85", "C86",
            "C87", "C88", "C89", "C90", "C91", "C92",
            "C93", "C94.0", "C94.1", "C94.2", "C94.3", 
            "C95", "C96.0", "C96.2", "C96.4", "C96.9",
            "C96.A", "C96.Z", "D45"),
    "asthma" = c("J45"),
    "cardiovascular" = generate_icd_prefixes("I", 0, 99,
            c(3, 4, 17, 18, 19, 29, 53, 54, 55, 
            56, 57, 58, 59, 90, 91, 92, 93, 94)),
    "chromosomal_anomalies" = generate_icd_prefixes("Q", 90, 99, c(94)),
    "conduct_disorders" = c("F91"),
    "congenital_malformations" = generate_icd_prefixes("Q", 0, 89,
            c(8, 9, 19, 29, 46, 47, 48, 49, 57, 58, 59)),
    "depression" = c("F32", "F33", "F06.30", "F06.31", "F06.32",
            "F34.9", "F39", "F43.21", "F43.23"),
    "developmental_delays" = c("F81", "R48.0", "F80", "H93.25",
            "F82", "F88", "F89", "F70", "F71", "F72", "F78", "F79"),
    "diabetes_mellitus" = c("E08", "E09", "E10", "E11", "E13"),
    "drug_abuse" = c(generate_icd_prefixes("F", 11, 19, c(17)), "F55", "Z71.51"),
    "eating_disorders" = c("F50"),
    "epilepsy" = c("G40", "R56"),
    "gastrointestinal" = c(generate_icd_prefixes("K", 20, 31, c(24)), 
            "K50", "K51", "K52", "K58", "Z87.1", "K92", "K62"),
    "joint_disorders" = c("M21", "M24"),
    "menstrual_disorders" = c("N91", "N92"),
    "nausea_and_vomiting" = c("R11"),
    "pain_conditions" = c("G89", "R52", "R10", "M54", "R07",
            "M25.5", "F45.4"),
    "psychotic_disorders" = c("F30", "F31", "F06.33", "F20", "F22", 
            "F23", "F24", "F25", "F28", "F29", "R44"),
    "sleep_disorders" = c("F51", "G47.0", "G47.1", "G47.2", "G47.3", 
            "G47.4", "G47.5", "G47.6","G47.8", "G47.9"),
    "smoking" = c("F17", "T65.2", "Z87.891"),
    "weight_loss" = c(generate_icd_prefixes("E", 40, 46, c()), "E64.0", "R63.4", "R64")
)

#Load diagnoses file.
diagnoses <- as.data.frame(fread(paste0(loadpath, "LaCava_Diagnoses_May7.csv")))

#Get the conditions associated with each diagnosis code.
unique_codes <- unique(diagnoses$DIAGNOSIS_CODE)
corresponding_conditions <- purrr::map_vec(unique_codes, ~find_condition(.x, paed_icd_code_dict), .progress = TRUE)
codes_with_conditions <- which(!is.na(corresponding_conditions))
unique_codes <- unique_codes[codes_with_conditions]
corresponding_conditions <- corresponding_conditions[codes_with_conditions]
names(corresponding_conditions) <- unique_codes #Create a dictionary linking conditions to codes.

#Attach a diagnosis to each diagnosis on the diagnoses table.
diagnoses$condition <- purrr::map_vec(diagnoses$DIAGNOSIS_CODE, ~ifelse(.x %in% unique_codes, corresponding_conditions[[.x]], NA), .progress = TRUE)

#Discard all diagnoses with no PCI conditions.
diagnoses <- diagnoses %>% filter(!is.na(condition)) %>% select(CSN, MRN, DIAGNOSIS_DATE, condition) %>%
    rename(csn=CSN, mrn=MRN, diagnosis_date=DIAGNOSIS_DATE)

print(head(diagnoses))

#Save this modified diagnosis table.
write.csv(diagnoses, paste0(savepath, "intermediate-files/diagnoses-with-conditions.csv"))

# # #Load this back out of memory.
diagnoses <- as.data.frame(fread(paste0(savepath, "intermediate-files/diagnoses-with-conditions.csv")))
visits <- as.data.frame(fread(paste0(savepath, "intermediate-files/visits-with-weight-and-vitals.csv")))
print(paste("After linking weight and vitals, we have", nrow(visits), "visits from", length(unique(visits$mrn)), "unique patients."))


#Identify arrival time and link this to diagnoses.
arrival_times <- visits %>% select(csn, arrival_time, mrn)
diagnoses <- diagnoses %>% select(-mrn) %>% inner_join(arrival_times, by="csn")

#Identify diagnoses which took place on a day BEFORE the day of the visit, since some diagnoses are recorded as MIDNIGHT THE SAME DAY (i.e. a few hours before arrival)
#but don't seem pre-existing (e.g. sepsis). 

#Define weights attached to each condition.
paed_weight_dict = list(
        "any_malignancy" = 5,
        "depression" = 4,
        "diabetes_mellitus" = 4,
        "epilepsy" = 4,
        "drug_abuse" = 3, 
        "psychotic_disorders" = 3, 
        "anemia" = 2,
        "cardiovascular" = 2,
        "chromosomal_anomalies" = 2,
        "congenital_malformations" = 2,
        "menstrual_disorders" = 2,
        "smoking" = 2,
        "weight_loss" = 2, 
        "anxiety" = 1,
        "asthma" = 1,
        "conduct_disorders" = 1,
        "developmental_delays" = 1, 
        "eating_disorders" = 1,
        "gastrointestinal" = 1, 
        "joint_disorders" = 1,
        "nausea_and_vomiting" = 1,
        "pain_conditions" = 1,
        "sleep_disorders" = 1,
        "alcohol_abuse" = 1
    )

#Separate diagnoses which happen on days before the visit (at ANY POINT during the year preceding the visit) from diagnoses which happen within 1 day of the visit.
#First, get diagnoses which happen during the visit.
current_diagnoses <- diagnoses %>% filter(abs(as.numeric(difftime(ymd_hm(arrival_time), ymd_hm(diagnosis_date), units="days"))) < 1) %>% select(csn, condition) 
    

#Combine this into a pci_after_visit score-- measuring the severity of the diagnoses ASSOCIATED WITH THE VISIT 
post_visit_comorbidity_scores <- current_diagnoses %>% 
    distinct(csn, condition) %>% #Only count point score once per condition.
    mutate(condition_score= purrr::map_vec(condition, ~paed_weight_dict[[.x]])) %>% 
    group_by(csn) %>% summarise(pci_of_visit=sum(condition_score))
  
#Create a variables table describing individual conditions ASSOCIATED WITH THE VISIT
current_diagnoses <- current_diagnoses %>% group_by(csn, condition)  %>% summarise(condition_count=1) %>% 
    tidyr::pivot_wider(names_from=condition, names_prefix="current_diagnosis_", values_from=condition_count, values_fill=0)

#Attach both to the visit table.
visits <- visits %>% select(-V1) %>%
    left_join(current_diagnoses, by="csn") %>% left_join(post_visit_comorbidity_scores, by="csn")

#Now get all diagnoses within a year of that visit (per Ally.)
all_diagnoses <- diagnoses %>% select(-c(arrival_time, csn)) #drop the visit identifier
diagnoses_before_visit <- arrival_times %>% inner_join(all_diagnoses, by="mrn", relationship="many-to-many") %>% #All diagnoses corresponding to that patient.
    mutate(days_before_visit=as.numeric(difftime(ymd_hm(arrival_time), ymd_hm(diagnosis_date), units="days"))) %>%
    filter(days_before_visit >= 1, days_before_visit < 365.25) %>% select(csn, condition) #Keep diagnoses NOT made on the same day as the visit, but made within a year of the visit.
    
#Combine this into a pci_after_visit score-- measuring the severity of the diagnoses PREDATING THE VISIT 
pre_visit_comorbidity_scores <- diagnoses_before_visit %>% 
    distinct(csn, condition) %>% #Only count point score once per condition.
    mutate(condition_score= purrr::map_vec(condition, ~paed_weight_dict[[.x]])) %>% 
    group_by(csn) %>% summarise(pci_before_visit=sum(condition_score))

#Create a variables table describing individual conditions PREDATING THE VISIT
diagnoses_before_visit <- diagnoses_before_visit %>% 
    group_by(csn, condition) %>% summarise(condition_count=1) %>% #Record all conditions relvant to the visit within a year of the visit.
    tidyr::pivot_wider(names_from=condition, names_prefix="pre_diagnosis_", values_from=condition_count, values_fill=0)

visits <- visits %>%
    left_join(diagnoses_before_visit, by="csn") %>%
    left_join(pre_visit_comorbidity_scores, by="csn") 

#If a patient has no PCI score they have a PCI of 0 (because they have no observed conditions.)
visits$pci_before_visit <- ifelse(is.na(visits$pci_before_visit), 0, visits$pci_before_visit)
visits$pci_of_visit <- ifelse(is.na(visits$pci_of_visit), 0, visits$pci_of_visit)


#If a patient does not appear in these tables, they have no diagnoses.
diagnosis_cols <- colnames(visits)[grepl("diagnosis", colnames(visits))]
for (col in diagnosis_cols) {
    visits[[col]] <- ifelse(is.na(visits[[col]]), 0, visits[[col]])
}

#Save checkpoint visits
write.csv(visits, paste0(savepath, "intermediate-files/visits-with-diagnoses.csv"))


# #Link complaints and triage score from a different file;
#also record when triage was completed 
complaints <- as.data.frame(fread(paste0(loadpath, "LaCava_May8_Demographics.csv"))) %>%
    select(CSN, ED_COMPLAINT, TRIAGE_ACUITY, TRIAGE_COMPLETE_DT_TM) %>% 
    rename(csn=CSN, complaint=ED_COMPLAINT, triage_acuity=TRIAGE_ACUITY, triage_end_time=TRIAGE_COMPLETE_DT_TM) %>%
    mutate(corrected_complaint=gsub('["\r\n\032]', '', complaint)) %>% #delete some punctuation to avoid EOF errors
    select(-complaint) %>% rename(complaint=corrected_complaint)

visits <-  as.data.frame(fread(paste0(savepath, "intermediate-files/visits-with-diagnoses.csv")))
print(paste("Before discarding visits with no identifiable complaint (no CSN or multiple CSNs in the complaint table), we have", nrow(visits), "visits from", length(unique(visits$mrn)), "unique patients."))


#Discard visits with duplicate CSNs in the complaints table, and those with no recorded complaint.
complaints <- complaints %>% 
    group_by(csn) %>% mutate(count=n()) %>%
    filter(count==1) %>% select(-count)

#Join to visits.
visits <- visits %>%
    select(-complaint) %>% inner_join(complaints, by="csn")


print(paste("After discarding visits with no identifiable complaint (no CSN or multiple CSNs in the complaint table), we have", nrow(visits), "visits from", length(unique(visits$mrn)), "unique patients."))
write.csv(visits, paste0(savepath, "intermediate-files/visits-with-linked-raw-complaint.csv"))


# #TO DO: Crowdedness, complaints
#Load this out of memory:
visits <- read.csv(paste0(savepath, "intermediate-files/visits-with-linked-raw-complaint.csv")) %>%
    select(-c(V1, V1, mu, sigma))

#Define crowdedness: number of patients in the ED when you arrived
#AND a pseudo-NEDOCS: the number of patients multiplied by 6-triage acuity

print(head(visits))

#Convert triage scores to numeric form. (This will create 281 NAs.)
visits$triage_acuity <- as.numeric(visits$triage_acuity)
visits$time_to_triage <- as.numeric(difftime(ymd_hm(visits$triage_end_time), ymd_hm(visits$arrival_time), units="mins"))

#Get crowdedness scores and pseudo-NEDOCS scores: first, create specialised dataframes.
crowdedness_df <- visits %>% 
    select(csn, arbitrary_timestamp, ed_los, triage_acuity) %>%
    rename(arrival=arbitrary_timestamp) %>% mutate(departure=arrival+ed_los) %>%
    select(-ed_los) %>% arrange(arrival)

assert(nrow(crowdedness_df)==nrow(visits))
assert(length(unique(visits$csn))==nrow(visits))
assert(length(unique(crowdedness_df$csn))==nrow(visits))

#Now the dataframe above is arranged in ASCENDING ORDER OF ARRIVAL TIME,
#so visits earlier in the dataframe started BEFORE visits later in the dataframe

#We can leverage this to save ourselves 
score_results <- purrr::map(
    1:nrow(visits),
    function (i) {
        #Count visits which started earlier than the arrival time, and ended after;
        #we only need to check the departure times of visits which occur before the current one 
        #in the data frame, because of the ordering
        if (i==1) {
            return(list(0, 0))
        } else {
            present_at_arrival <- which(crowdedness_df$departure[1:(i-1)] > crowdedness_df$arrival[i])
        
            crowdedness_score <- length(present_at_arrival) #no offset, so present-at-arrival also gives indices in the original dataframe
            pseudo_nedocs <-  sum(6-crowdedness_df$triage_acuity[present_at_arrival], na.rm=TRUE)

            return(list(crowdedness_score, pseudo_nedocs)) #Crowdedness is number of patients there when you got to ED
                #pseudo-NEDOCS weights them by severity, such that a level-1 patient contributes 5 points
        }
        
    }, .progress=TRUE
)

#Set crowdedness score- remembering these are not in the saem order as the visits dataframe
crowdedness_df$crowdedness <- purrr::map_dbl(score_results, 1)
crowdedness_df$pseudo_nedocs <- purrr::map_dbl(score_results, 2)
crowdedness_df <- crowdedness_df %>% select(csn, crowdedness, pseudo_nedocs)
assert(nrow(crowdedness_df)==nrow(visits))
visits <- visits %>% left_join(crowdedness_df, by="csn")

#Save this in memory.
write.csv(visits, paste0(savepath, "intermediate-files/visits-with-linked-crowdedness.csv"))

#Read this back out of memory.
visits <- read.csv(paste0(savepath, "intermediate-files/visits-with-linked-crowdedness.csv"))

#Now deal with complaints:
print(head(visits))

#First, convert them all to lowercase

visits$complaint <- tolower(visits$complaint)

#Expand common abbreviations.
common_abbrev <- list(
    "uri" = c("upper respiratory infection"),
    "ent" = c("ear, nose or throat"),
    "si" = c("suicidal ideation"),
    "unk" = c("unknown"),
    "lwbs" = c("left without being seen"),
    "rsv" = c("respiratory syncytial virus"),
    "uti" = c("urinary tract infection"),
    "injryhead" = c("head injury"),
    "mva" = c("motor vehicle accident"),
    "sorethroat" = c("sore throat"),
    "psych" = c("psychiatric"),
    "eval" = c("evaluation"),
    "dka" = c("diabetic ketoacidosis"),
    "allerreact" = c("allergic reaction"),
    "urti" = c("upper respiratory tract infection"),
    "neuro" = c("neurological"),
    "nodm" = c("new onset diabetes mellitus"),
    "injryarm" = c("arm injury"),
    "inj" = c("injury"),
    "lac" = c("laceration"),
    "ftt" = c("failure to thrive"),
    "mvc" = c("motor vehicle collision"),
    "gi" = c("gastrointestinal"),
    "injryankle" = c("ankle injury"),
    "scrotpain" = c("scrotal pain"),
    "ed" = c("eating disorder"),
    "gtubeprob" = c("g-tube problem"),
    "prob" = c("problem"),
    "loc" = c("loss of consciousness"),
    "gu" = c("genitourinary"),
    "appy" = c("appendicitis"),
    "w/o" = c("without"),
    "rxn" = c("reaction"),
    "fb" = c("foreign body"),
    "vs" = c("versus"),
    "redeye" = c("red eye"),
    "injryabdm" = c("abdominal injury"),
    "ha" = c("headache"),
    "injryleg" = c("leg injury"),
    "injryelbow" = c("elbow injury"),
    "pna" = c("pneumonia"),
    "iv" = c("intravenous"),
    "rx" = c("reaction"),
    "obs" = c("other"), #this is a common abbrev but unclear what it means
    "swe" = c("swelling"),
    "strep" = c("streptococcal"),
    "injryshoul" = c("shoulder injury"),
    "resp" = c("respiratory"),
    "gtube" = c("g-tube"),
    "g tube" = c("g-tube"),
    "injryknee" = c("knee injury"),
    "injryeye" = c("eye injury"),
    "sbo" = c("small bowel obstruction"),
    "s/p" = c("status post"),
    "ukn" = c("unknown"),
    "fuo" = c("fever of unknown origin"),
    "injrynose" = c("nose injury"),
    "ng tube" = c("nasal gastric tube"),
    "r/o" = c("rule out"),
    "uc" = c("ulcerative colitis"),
    "vom" = c("vomiting"),
    "gen" = c("general"),
    "arfid" = c("eating disorder"), #avoidant/restrictive food intake disorder
    "gerd" = c("gastroesophageal reflux disease"),
    "ams" = c("altered mental status"),
    "aom" = c("acute otitis media"),
    "t&a" = c("tonsillectomy and adenoidectomy"),
    "t and a" = c("tonsillectomy and adenoidectomy"),
    "t/a" = c("tonsillectomy and adenoidectomy"),
    "re-eval" = c("re-evaluation"),
    "dvt" = c("deep vein thrombosis"),
    "po" = c("feeding"),
    "po's" = c("feeding"),
    "std" = c("venereal disease"),
    "svt" = c("superventricular tachycardia"),
    "abnl" = c("abnormal"),
    "abd" = c("abdominal"),
    "cf" = c("cystic fibrosis"),
    "flu" = c("influenza"),
    "diff" = c("difficulty"),
    "wob" = c("work of breathing"),
    "osa" = c("obstructive sleep apnea"),
    "asd" = c("autism spectrum disorder"),
    "nurs el" = c("nursemaid's elbow"),
    "acs" = c("acute coronary syndrome"),
    "ibd" = c("irritable bowel disorder"),
    "mis-c" = c("multisystem inflammatory syndrome in children"),
    "s/a" = c("suicide attempt"),
    "sa" = c("suicide attempt"),
    "behav" = c("behavioral"),
    "fx" = c("fracture"),
    "vp" = c("ventriculoperineal"),
    "itp" = c("immune thrombocytopenic purpura"),
    "infect" = c("infection"),
    "aki" = c("acute kidney injury"),
    "injryear" = c("ear injury"),
    "nursemaids" = c("nursemaid's"),
    "rlq" = c("right lower quadrant"),
    "all" = c("acute lymphoblastic leukemia"),
    "b-all" = c("b cell acute lymphoblastic leukemia"),
    "scfe" = c("slipped cap femoral epiphysis"),
    "scd" = c("sickle cell disease"),
    "auto" = c("motor vehicle"),
    "vs" = c("versus"),
    "v." = c("versus"),
    "ped" = c("pedestrian"),
    "gsw" = c("gunshot wound"),
    "hsv" = c("herpes simplex virus"),
    "htn" = c("hypertension"),
    "pid" = c("pelvic inflammatory disorder"),
    "fnd" = c("functional neurological disorder"),
    "s.o.b" = c("shortness of breath"),
    "sob" = c("shortness of breath"),
    "v/d" = c("venereal disease"),
    "aiha" = c("autoimmune hemolytic anemia"),
    "aml" = c("acute myeloid leukemia"),
    "pede" = c("pedestrian"),
    "atrt" = c("atypical teratoid rhabdoid tumor"),
    "egd" = c("upper endoscopy"),
    "iddm" = c("insulin-dependent diabetes mellitus"),
    "dx" = c("diagnosis"),
    "mrsa" = c("methicillin-resistant staphylococcus aureus"),
    "od" = c("overdose"),
    "w/" = c("with"),
    "r" = c("right"),
    "ro" = c("rule out"),
    "tbi" = c("traumatic brain injury"),
    "bpd" = c("bipolar disorder"),
    "cvl" = c("central venous line"),
    "shuntmalf" = c("shunt malfunction"),
    "t1dm" = c("type 1 diabetes mellitus"),
    "tib fib" = c("tibula and fibula"),
    "sx" = c("symptoms"),
    "dz" = c("disease"),
    "avm" = c("arteriovenous malformation"),
    "789.00" = c("abdominal pain"),
    "hsp" = c("henoch-schönlein purpura"),
    "wt" = c("weight")
)

#Wherever we find these abbreviations (on their own, not within a word), replace them with their expansions.
for (abbrev in names(common_abbrev)) {
    pattern <- paste0("\\b", abbrev, "\\b")
    visits$complaint <- str_replace_all(visits$complaint, regex(pattern, ignore_case = TRUE), common_abbrev[[abbrev]][1]) #   
}

#Split complaints by commas and order them in descending frequency
all_complaints <- as.data.frame(table(unlist(str_split(visits$complaint, ",")))) %>% 
    rename(ComplaintPhrase=Var1) %>% arrange(desc(Freq))
write.csv(all_complaints, paste0(savepath, "intermediate-files/complaint-phrases-by-freq.csv"))


#Now define pediatric complaints according to PERC clusters.
perc_complaints <- list(
    "abdominal_pain" = c("flank pain", "abdominal pain", "appendicitis",
        "appendix", "abdominal problem"),
    "assault" = c("assault", "abuse"),
    "allergic_reaction" = c("allergic reaction"),
    "altered_mental_status" = c("altered mental status", "confusion"), 
    "asthma_or_wheezing" = c("asthma", "wheezing", "status asthmaticus"),
    "bites_or_stings" = c("bite", "sting"),
    "burn" = c("burn"), 
    "cardiac" = c("tachycardia", "palpitations", "cardiac",
        "respiratory arrest", "brue",
        "cardiorespiratory arrest", "heart",
        "bradycardia"), #heart pain not found
    "chest_pain" = c("chest pain", "rib pain"),
    "chronic_disease" = c("sickle", "cancer", "leukemia", "diabetes",
        "diabetic", "hyperglycemia", "hypoglycemia", "bleeding disorder",
        "tumor", "neoplasm"),
    "congestion" = c("upper respiratory infection", 
        "upper respiratory tract infection", "congestion"),
    "constipation" = c("constipation", "dysuria", 
        "difficulty urinating", "urinary retention"),
    "cough" = c("cough", "hemoptysis", "coughing up blood", "bronchitis",
        "bronchiolitis"),
    "croup" = c("croup"),
    "crying_or_colic" = c("crying", "colic", "fussy",
        "fussiness", "irritable", "irritability"),
    "dental" = c("tooth", "dental", "gums", "mouth pain",
        "mouth injury"),
    "device_complication" = c("device", "equipment", "prosthetic",
        "tube", "g-tube"),
    "diarrhea" = c("diarrhea", "intestinal infection", 
        "gastroenteritis", "gastrointestinal problem"),
    "ear_complaint" = c("ear", "otitis media",
        "earache"),
    "epistaxis" = c("epistaxis", "nose", "nasal"),
    "extremity" = c("arm injury", "ankle injury",
        "leg injury", "wrist injury", "finger injury",
        "knee injury", "elbow injury", "foot injury",
        "leg pain", "hand injury", "knee pain", 
        "ankle pain", "pain in limb", "injury - arm", 
        "toe injury", "injury - hand/fingers", 
        "arm injury - major", "injury - wrist",
        "leg pain", "injury - foot/toes", "hand pain",
        "ankle sprain", "injury - ankle", "injury hand",
        "foot pain"),
    "eye_complaint" = c("conjunctivitis", "eye", "visual"),
    "syncope" = c("syncope", "fainting"),
    "foreign_body" = c("foreign body"), #combining original nasal, skin
    "fever" = c("fever", "febrile"), #also including fever in neonate here
    "follow_up" = c("follow up", "abnormal lab", "abnormal labs",
        "dressing", "cast", "check", "wound evaluation"),
    "general" = c("fatigue", "weakness", "viral illness", "influenza",
        "lethargy"), #include this high frequency pain complaint under '
    "gi_bleed" = c("gastrointestinal bleed", "hematemesis",
         "blood in vomit", "blood in stool", "bloody stool",
        "hematuria", "blood in urine"),
    "gynecologic" = c("vaginal", "amenorrhea", "gynecological"),
    "head_or_neck" = c("head injury", "injury - head",
        "injury head", "head laceration", "head trauma", "injury-head",
        "face", "face injury", "cellulitis of face", "facial laceration",
        "facial injury", "facial swelling", "neck injury"),
    "headache" = c("face pain", "facial pain", "headache", "migraine"),
    "laceration" = c("laceration", "puncture wound"),
    "lump_or_mass" = c("lump", "mass", "swelling"), #may include edema if recorded as 'swelling or edema'
    "male_genital" = c("penile", "penis", "testicular",
        "testes", "penis", "genitourinary"),
    "mvc" = c("motor vehicle", "car", "pedestrian", "versus", "collision"),
    "neck_pain" = c("neck pain"),
    #exclude other
    "neurologic" = c("walking", "weakness",
        "dizzy", "dizziness", "speech", "bell's", "neurologic problem"),
    "poisoning" = c("alcohol", "ingestion", "overdose", "substance", "poisoning"),
    "poor_feeding" = c("dehydration", "failure to thrive",
         "poor nutritional intake", "weight loss"),
    "pregnancy" = c("pregnancy", "pregnant"),
    "primary_care" = c("vaccine", "medical screening exam", "refill"),
    "psych" = c("psychiatric", "mental health", "anxiety", "eating disorder",
        "depression", "agitation", "violent", "aggressive", "aggression",
        "delusion", "hallucination", "suicidal", "suicide", "inflicted",
         "intentional", "self"),
    "rash" = c("rash", "acne", "wound infection", "cellulitis", "skin",
        "eczema", "abcess", "abscess"),
    "other_respiratory" = c("shortness of breath", "difficulty breathing", 
        "tachypnea", "stridor", "choking", "pneumonia", "respiratory"),
    "seizure" = c("seizure", "seizures", "epilepsy"),
    "sore_throat" = c("sore throat", "throat problem",
         "tonsil", "throat pain", "tonsillitis"),
    "trauma" = c("trauma", "fall", "back pain"),
    "urinary" = c("painful urination", "urinary tract infection",
        "kidney"),
    "vomiting" = c("vomiting", "vomitting")
)

#The presence of any of these tags defines a complaint;
#a visit may be described by multiple complaints
for (complaint in names(perc_complaints)) {
    complaint_name <- paste0("complaint_contains_", complaint)
    tags <- perc_complaints[[complaint]]
    visits[[complaint_name]] <- 0
    for (tag in tags) {
        visits[[complaint_name]] <- ifelse(grepl(paste0("\\b", tag, "\\b"), visits$complaint) | visits[[complaint_name]]==1, 1, 0)
    }
}

print(head(visits))

write.csv(visits, paste0(savepath, "intermediate-files/visits-with-complaints.csv"))
visits <- as.data.frame(fread(paste0(savepath, "intermediate-files/visits-with-complaints.csv")))

#Filter one visit from 2012.
visits <- visits %>% 
    filter(year_of_arrival >= 2019, year_of_arrival <= 2024) %>%
    filter(ed_los < 27*24*60) #length of stay increases more or less continuously up to 26.7 days and then jumps to 46 (and can be 10k days)


#We have 331751 visits where time-to-room is 0 (i.e. patients are roomed immediately) and 221 where time to room is negative (a maximum of 14 mins before
#the arrival time). We assume this means the arrival time was 'basically immediate' and round this up to 0.



print(paste("After filtering those in the wrong date range and those with an outlier ED LOS, we have", nrow(visits), "visits from", length(unique(visits$mrn)), "unique patients."))

visits <- visits %>%
    filter(is_admitted==0 | (is_admitted==1 & time_to_admit_request >=0)) #Discard those with nonsensical time-to-admit-request times

print(paste("After filtering those with negative admit request times (if admitted), we have", nrow(visits), "visits from", length(unique(visits$mrn)), "unique patients."))

visits <- visits %>%
    filter(is_admitted==0 | (is_admitted==1 & time_from_request_to_admission >=0)) #Discard those with nonsensical time-from-decision-to-departure times

print(paste("After filtering those with negative request-to-departure-times (if admitted), we have", nrow(visits), "visits from", length(unique(visits$mrn)), "unique patients."))

#At least one visit is registered under 2 different CSNs with the same MRN, arrival and departure time. So drop all visits registered more than once (same name and DOB, as well).
#That visit happened on 30 May 24, the day before the EPIC transfer.

visits <- visits %>% group_by(mrn, arrival_time) %>%
    mutate(count=n()) %>% ungroup() %>% filter(count==1) %>%
    select(-count)

print(paste("After filtering visits duplicated under multiple CSNs, we have", nrow(visits), "visits from", length(unique(visits$mrn)), "unique patients.")) #We don't lose anything here!

#Drop unnecessary columns.
visits <- visits %>% select(-c(X.1, mrn, zipcode, ethnicity, admission_request_time,
    rooming_time, arrival_time, disposition, departure_time, gender,
    triage_end_time, V1))


write.csv(visits, paste0(savepath, "preprocessed-visits.csv"))

# #Check how many visits are described by at least one complaint.
# complaint_cols <- colnames(visits)[startsWith(colnames(visits), "complaint_contains_")]
# described_by_any_complaint <- ifelse(rowSums(visits[complaint_cols])>0, 1, 0)
# print(paste("Fraction described by any complaint", sum(described_by_any_complaint)/length(described_by_any_complaint)))

# print(head(visits))
