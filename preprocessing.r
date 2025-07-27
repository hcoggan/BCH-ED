
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
library(vroom)
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
library(data.table)
library(readr)
library(SnowballC)


setwd("your/working/directory")
save_filepath <- "your/save/path"


#SUPPLEMENTARY FUNCTIONS:

#Convert age in days into age-groups.
get_age_group <- function(age_in_years, age_in_days) {
    if (age_in_years >= 15) {
        age_group <- "older_than_15_years"
    }
    else {
        if (age_in_years >= 10) {
            age_group <- "ten_to_15_years"
        } 
        else {
            if (age_in_years >= 5) {
            age_group <- "five_to_10_years"
        } else {
            if (age_in_years > 3) {
                age_group <- "three_to_5_years"
            }
            else {
               if (age_in_days > 365.25*3/2) {
                    age_group <- "eighteen_months_to_3_years"
               } else {
                    if (age_in_days > 365.25) {
                        age_group <- "twelve_to_18_months"
                    } else {
                        if (age_in_days > 365.25/2) {
                            age_group <- "six_to_12_months"
                        } else {
                            if (age_in_days > 365.25/4) {
                                age_group <- "three_to_6_months"
                            } else {
                                age_group <- "under_3_months"
                            }
                            }
                        }
                    }
               }
            }
        }
    
        
    }
    return(age_group)
}

#Identify state from a patient's ZIP code.
get_state <- function(zip_code) {
    ifelse(nchar(zip_code)==5, reverse_zipcode(zip_code)$state, NA)
}

#Retrieve the state from a patient's zip code (given a list of unique zip codes, and a dictionary mapping those to states.)
retrieve_state <- function(zip_code, unique_zip_codes, states) {
    ifelse(zip_code %in% unique_zip_codes, states[[zip_code]], "unknown")
}




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

#Get the condition associated with a code.
get_condition <- function(code, unique_codes, corresponding_conditions) {
        ifelse(code %in% unique_codes, corresponding_conditions[[code]], NA)
}
#Get the PCI score associated with a condition
get_score <- function(condition, paed_weight_dict) {
        paed_weight_dict[[condition]]
}

#Create a flag for PALS BP criteria-- is a patient ever hypotensive during their stay?
hypotensive <- function(sbp, age_in_months) {
        (age_in_months <= 1 & sbp < 60) | (age_in_months <= 12 & age_in_months > 1 & sbp < 70) | (age_in_months <= 120 & age_in_months > 12 & sbp < (70 + (age_in_months/6))) | (age_in_months > 120 & sbp < 90)  
    }


#Extract words from a chief complaint.
extract_words <- function(complaint, all_words) {
    list_of_words <- stringr::word(str_replace_all(complaint, "[[:punct:]]", " "))
    all_words <- c(all_words, list_of_words)

}

#Tag each complaint by the appearance of certain stems.
convert_complaints <- function(i, baseline_factors, stems_to_use) {
    stems_in_complaint <- wordStem(stringr::word(str_replace_all(baseline_factors$ed_complaint[i], "[[:punct:]]", " ")))
    appears_in <- ifelse(stems_to_use %in% stems_in_complaint, 1, 0)
    res <- as.data.frame(setNames(as.list(appears_in), paste0("complaint_contains_", stems_to_use)))
    res$csn <- baseline_factors$csn[i]
    return(res)
}


#Clean names of medications
clean_value_names <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9_]", "_") %>%
    str_replace_all("_{2,}", "_") %>%
    str_remove("^_|_$")
}

#Count the number of patients there at the start of a visit.
get_patients_present_at_arrival <- function(t, arrival_times, departure_times) {
    return(sum(arrival_times < t & departure_times > t))
}


#MAIN FUNCTION: Link demographic data from across files.
link_demographics <- function() {

	demographics <- read.csv("demographics_8May.csv")

    print(paste("OG visits", nrow(demographics)))

	#link race and keep only unique visits
	demographics_with_race <- read.csv("demographics_1May.csv")
	demographics_with_race <- demographics_with_race %>% select(CSN, RACE, AGE) %>% group_by(CSN) %>% mutate(count=n()) %>% filter(count==1) %>% select(-count) %>% ungroup()
	demographics <- demographics %>% select(-c(GENDER.1, ZIP_CODE.1)) %>% group_by(CSN) %>% mutate(count=n()) %>% filter(count==1) %>% select(-count) %>% ungroup()
	demographics <- inner_join(demographics, demographics_with_race, by="CSN")

    print(paste("Visits with unique ID", nrow(demographics)))


    #link age in days
    demographics$age_in_days <- floor(difftime(as.POSIXct(demographics$ED_ARRIVAL_TIME, format="%Y-%m-%d"), as.POSIXct(demographics$DOB, format="%Y-%m-%d"), units="days"))
    demographics <- demographics %>% select(-c(DOB, NAME))

    #GENDER
    #for consistency with MIMIC, use patient-reported sex
    #use a flag indicating whether or not someone is trans

    #there are about a hundred people whose sex is M and gender is "transgender M" and vice versa; we assume thismeans they are actually trans and their sex has been mis-recorded
    demographics$SEX <- ifelse(demographics$GENDER=="Transgender Female", "M", ifelse(demographics$GENDER=="Transgender Male", "F", demographics$SEX))


    #<20 people have no recorded sex (U or X); we discard them
    demographics <- demographics %>% filter(!(SEX=="U" | SEX =="X"))

    print(paste("Visits with identifiable sex", nrow(demographics)))


    #'Other', 'Trans', or "NB" will get you this marker, as will a reported gender other than your sex
    #We do not assume you are trans if you have "choose not to answer, unable to collect or don't know"
    demographics$is_trans_or_nb <- ifelse((demographics$GENDER=="Nonbinary" | (demographics$SEX == "M" & (demographics$GENDER == "Female" | demographics$GENDER == "Transgender Female")) |
            (demographics$SEX == "F" & (demographics$GENDER == "Male" | demographics$GENDER == "Transgender Male")) | demographics$GENDER == "Nonbinary (e.g. genderqueer/gender nonconforming)" |
            (demographics$GENDER == "Other")), 1, 0)

    #drop gender beyond this
    demographics <- demographics %>% select(-GENDER)

    #RACE AND ETHNICITY

    #if someone's race is "Another race, non-Hispanic" and their ethnicity contradicts this, we adjust it
    #we also use ethnicity to fill blanks for Unknown race

    #e.g. if someone says their race is Other but gives an Asian ethnicity, we adjust race to Asian

    demographics$RACE <- ifelse((demographics$RACE == "Another Race, non-Hispanic" | demographics$RACE == "Unknown" ) & (demographics$ETHNICITY=="African American" | demographics$ETHNICITY=="Black"), "Black, non-Hispanic", demographics$RACE)
    demographics$RACE <- ifelse((demographics$RACE == "Another Race, non-Hispanic" | demographics$RACE == "Unknown" ) & (demographics$ETHNICITY=="Asian Indian" | demographics$ETHNICITY=="Chinese" |
            demographics$ETHNICITY=="Hmong" | demographics$ETHNICITY=="Japanese" | demographics$ETHNICITY=="Korean" | demographics$ETHNICITY=="Laotian/Lao" | 
            demographics$ETHNICITY=="Okinawan" | demographics$ETHNICITY=="Pakistani" | demographics$ETHNICITY=="Sri Lankan" | demographics$ETHNICITY=="Taiwanese" |
            demographics$ETHNICITY=="Thai" |  demographics$ETHNICITY=="Vietnamese" |  demographics$ETHNICITY=="Bangladeshi" |  demographics$ETHNICITY=="Cambodian"), "Asian, non-Hispanic", demographics$RACE)
    demographics$RACE <- ifelse((demographics$RACE == "Another Race, non-Hispanic" | demographics$RACE == "Unknown") & (demographics$ETHNICITY=="Hispanic or Latino"), "Hispanic", demographics$RACE) 


    demographics$RACE <- case_match(demographics$RACE, 
            "Another Race, non-Hispanic" ~ "Other",
            "Asian, non-Hispanic" ~ "Asian",
            "Black, non-Hispanic" ~ "Non-Hispanic Black",
            "Hispanic" ~ "Hispanic",
            "Multiracial, non-Hispanic" ~ "Non-Hispanic Multiracial",
            "Unknown" ~ "Unknown",
            "White, non-Hispanic" ~ "Non-Hispanic White")


    #drop ethnicity after this point
    demographics <- demographics %>% select(-ETHNICITY)


    #Convert times to the correct format
    demographics$ED_ARRIVAL_TIME <- as.POSIXct(demographics$ED_ARRIVAL_TIME, format = "%Y-%m-%d %H:%M")
    demographics$ED_CHECKOUT_TIME <- as.POSIXct(demographics$ED_CHECKOUT_TIME, format = "%Y-%m-%d %H:%M")
    demographics$TRIAGE_START_DT_TM <- as.POSIXct(demographics$TRIAGE_START_DT_TM, format = "%Y-%m-%d %H:%M")
    demographics$ADMIT_DATE <- as.POSIXct(demographics$ADMIT_DATE, format = "%Y-%m-%d %H:%M")

    #Assign each visit an age group.
    demographics$age_group <- mapply(get_age_group, demographics$AGE, demographics$age_in_days)

    #The supplied field relevant to age is just age in years; we rename this accordingly.
    demographics <- demographics %>% rename(age_in_years=AGE)

    #MODE OF ARRIVAL
    demographics$ED_ARRIVAL_MODE <- case_match(demographics$ED_ARRIVAL_MODE, 
            c("", "Unknown", "Other")  ~ "Unknown",
            c("EMS", "Ambulance", "Air transport", "Ambulance: Other EMS", "Medical Flight", "Police") ~ "EMS",
            c("Critical Care", "Hospital Transport", "Transfer") ~ "Transfer",
            c("Walk", "Walk in", "Wheelchair", "Car", "Public Transportation", "Taxi") ~ "Walk in")

    #NOW DISPOSITION-- keep those whose endpoints are 'other' (leaving AMA for now), discard those whose outcomes are unknown
    demographics$is_admitted <- ifelse(demographics$DERIVED_DESPOSITON == "Admit" | demographics$DERIVED_DESPOSITON == "ED Patient Admitted", 1, 0)
    demographics$is_discharged <- ifelse(demographics$DERIVED_DESPOSITON == "Discharge" | demographics$DERIVED_DESPOSITON == "Home", 1, 0)

    demographics <- demographics %>% filter(!((DERIVED_DESPOSITON == "Unknown" | DERIVED_DESPOSITON == ""))) %>% select(-DERIVED_DESPOSITON)

    print(paste("Visits with identifiable disposition", nrow(demographics)))


    #PREFERRED LANGUAGE
    #keep English, Spanish and other
    demographics$primary_language <- ifelse(demographics$PREFERRED_LANGAUAGE=="English", "English", ifelse(demographics$PREFERRED_LANGAUAGE=="Spanish", "Spanish", "Other"))
    demographics <- demographics %>% select(-PREFERRED_LANGAUAGE)

    #INSURANCE
    #convert into Medicaid and non-Medicaid, which is basically all private-- an insurance plan is Medicaid if its name includes one of a number of key substrings
    medicaid_substrings <- c("MEDICAID", "PUBLIC", "MASSHEALTH", "ACO", "COMMUNITY")
    demographics$has_medicaid <- ifelse(grepl(paste(medicaid_substrings, collapse="|"), demographics$PRIMARY_INSURANCE_PAYOR), 1, 0)
    demographics <- demographics %>% select(-PRIMARY_INSURANCE_PAYOR)

    #Save a list of compliants with associated frequencies.
    tab <- table(toupper(demographics$ED_COMPLAINT))
    write.csv(tab[order(tab, decreasing = TRUE)], "complaints.csv")

    #Group complaints into rough categories. This is less granular than the formulation described in the paper, and is just to get
    #a rough sense of the distribution of complaints.
    complaint_dict <- list(
        "abdominal_pain" = c("ABD"), #these three letters are sufficient to identify any
        "fever" = c("FEVER"),
        "cough" = c("COUGH", "CROUP"),
        "resp_problem" = c("RESPIRATORY", "WHEEZING", "SHORTNESS OF BREATH", "ASTHMA", "CHOKING", "STRIDOR"),
        "nvd" = c("NAUSEA", "VOMITING", "DIARRHEA", "FLU-LIKE", "EMESIS"),
        "chi" = c("CLOSED HEAD INJURY"),
        "loc" = c(" LOC"),
        "headache" = c("HEADACHE"),
        "skin_problem" = c("RASH", "SKIN", "CELLULITIS"),
        "chest_and_cardiac" = c("CHEST PAIN", "YCARDIA", "CARDIO"),
        "injury" = c("LACERATION", "FALL", "INJURY", "ACCIDENT", "ASSAULT", "WOUND", "PEDESTRIAN", "VEHICLE", "ABUSE", "TRAUMA"),
        "seizure" = c("SEIZURE"),
        "ent" = c("EAR", "NOSE", "THROAT", "ENT PROBLEM", "EPISTAXIS", "ORAL LESIONS"),
        "neuro" = c("NEUROLOGIC", "FACIAL DROOP", "TWITCHING", "ALTERED GAIT", "PARESTHESIA"),
        "eye_problem" = c("EYE", "VISION"),
        "follow_up" = c("ABNORMAL LAB", "EQUIPMENT", "SCREENING", "TUBE", "GENERAL MEDICAL", "REEVALUATION", "MEDICAL PROBLEM", "MALFUNCTION", "MEDICATION REFILL", "EVALUATION"),
        "si_and_sh" = c("SUICIDAL", "INTENT", "SELF INFLICTED"),
        "genitourinary" = c("CONSTIPATION", "GENITOURINARY", "TESTICULAR", "DYSURIA", "PENILE", "PENIS", "VAGIN", "GROIN", "URINARY", "ANURIA", "INCONTINENCE"),
        "exposure" = c("ALLERG", "BITE", "BURN", "TOXIC", "BEE STING", "EXPOSURE"),
        "psych" = c("AGGRESSION", "BEHAVIOR", "PSYCH", "ANXIETY", "DEPRESSION", "ALTERED MENTAL STATUS", "ANOREXIA"),
        "pain_swelling" = c("DENTAL", "ABSCESS", "PAIN", "SWELLING", "CYST", "EDEMA"), #this will get all pain and swelling
        "weak_and_dizzy" = c("SYNCOPE", "DIZZINESS", "LETHARGY", "WEAKNESS", "FATIGUE"),
        "fussy" = c("FUSS", "IRRITAB", "CRYING"),
        "gastro_problem" = c("STOOL", "BLOOD IN VOMIT", "GI BLEEDING", "GASTROINTESTINAL"), 
        "foreign_body" = c("FOREIGN"),
        "swallowing_and_eating_problems" = c("NUTRITIONAL", "SWALLOWING"),
        "homeless" = c("HOMELESS")
    )


    #condense complaints down to this list-- give a binary indicator as to whether chief complaint contains ANY of these subtags
    for (tag in names(complaint_dict)) {
        name <- paste0('complaint_', tag)
        demographics[[name]] <- ifelse(grepl(paste(complaint_dict[[tag]], collapse="|"), toupper(demographics$ED_COMPLAINT)), 1, 0) 
    }


    #ZIP CODES
    #Identify unique zip codes, and map them to states.
    unique_zip_codes <- unique(demographics$ZIP_CODE)
    states <- purrr::map_vec(unique_zip_codes, get_state, .progress = TRUE)

    valid_state <- which(!is.na(states))
    unique_zip_codes <- unique_zip_codes[valid_state]
    states <- states[valid_state]

    names(states) <- unique_zip_codes
    demographics$state <- purrr::map_vec(demographics$ZIP_CODE, ~retrieve_state(.x, unique_zip_codes, states), .progress = TRUE)
    demographics$state <- ifelse(!(demographics$state == "MA"), ifelse(demographics$state == "unknown", "unknown", "out_of_state"), "in_state")

    #Use zip code to extract miles travelled
    miles_travelled <- as.numeric(zip_distance(demographics$ZIP_CODE, "02115")$distance)

    #Use the Robert Graham Center's Social Deprivation Index, linked by zip code-- we only have latest data which is 2019
    sdi_scores <- read.csv("rgcsdi-2015-2019-zcta.csv")

    sdi_scores <- sdi_scores %>% select(ZCTA5_FIPS, SDI_score)
    sdi_scores$ZCTA5_FIPS <- as.character(sdi_scores$ZCTA5_FIPS)
    sdi_scores$ZCTA5_FIPS <- ifelse(nchar(sdi_scores$ZCTA5_FIPS)==4, paste0("0", sdi_scores$ZCTA5_FIPS), sdi_scores$ZCTA5_FIPS)

    demographics <- left_join(demographics, sdi_scores, by=c("ZIP_CODE"="ZCTA5_FIPS"))
        
    #keep SDI score, lose zip code
    demographics <- demographics %>% select(-c(ZIP_CODE))

    #Add new value indicating an SDI score is unknown
    demographics$SDI_score <- ifelse(is.na(demographics$SDI_score), "unknown", demographics$SDI_score)


    #Now assign index quantiles- 1 is richest, 5 is poorest. This is not used in OR analysis and is just to get a sense of demographics.
    demographics$sdi_quantile <- case_when(
            demographics$SDI_score < 13 ~ "1",
            demographics$SDI_score >= 13 & demographics$SDI_score < 27 ~ "2",
            demographics$SDI_score >= 27 & demographics$SDI_score < 46 ~ "3",
            demographics$SDI_score >= 46 & demographics$SDI_score < 70 ~ "4",
            demographics$SDI_score >= 70  ~ "5", #weighting by patient instead of by zip code doesn't help; it's just that most visits come in the form of repeated visits from poorer patients 
            .default = "unknown"
    )

    demographics$cts_SDI_score <- demographics$SDI_score
    demographics <- demographics %>% select(-c(SDI_score))

    #Divide 'miles travelled' into broad categories, and save the continuous variable under a specific name.
    #Again, not used in the analysis.
    demographics$miles_travelled <- NA
    demographics$miles_travelled <- ifelse(miles_travelled < 5, "less_than_5_miles", demographics$miles_travelled)
    demographics$miles_travelled <- ifelse(miles_travelled >= 5 & miles_travelled < 10, "five_to_10_miles", demographics$miles_travelled)
    demographics$miles_travelled <- ifelse(miles_travelled >= 10 & miles_travelled < 20, "ten_to_20_miles", demographics$miles_travelled)
    demographics$miles_travelled <- ifelse(miles_travelled >= 20, "more_than_20_miles", demographics$miles_travelled)
    demographics$miles_travelled <- ifelse(is.na(demographics$miles_travelled), "unknown", demographics$miles_travelled)

    demographics$cts_miles_travelled <- ifelse(is.na(miles_travelled), 'unknown', miles_travelled)


    #get WEIGHT normalised by age
    for (group in unique(demographics$age_group)) {
            indices <- which(demographics$age_group==group)
            mu <- mean(as.numeric(demographics$WEIGHT_KG[indices]), na.rm=TRUE)
            sigma <- sd(as.numeric(demographics$WEIGHT_KG[indices]), na.rm=TRUE)
            demographics$WEIGHT_KG[indices] <- ifelse(is.na(demographics$WEIGHT_KG[indices]), "not_taken", abs((as.numeric(demographics$WEIGHT_KG[indices])-mu)/sigma))
    } 

    demographics$cts_weight <- demographics$WEIGHT_KG #record continuous weight, taken as number of sds from mean (so doesn't distinguish between the underweight and the overweight)

    #Now convert this into categories.
    demographics$weight <- case_when(
            demographics$WEIGHT_KG < 1 ~ "less_than_1_sd",
            demographics$WEIGHT_KG < 2 & demographics$WEIGHT_KG >=1 ~ "between_1_and_2_sd",
            demographics$WEIGHT_KG > 2 ~ "more_than_2_sd",
            .default = "unknown")


    #Drop original complaint text, and while we're here, drop death indicators, height and weight
    demographics <- demographics %>% select(-c(ED_COMPLAINT, DEATH_IND, WEIGHT_KG, WEIGHT_DATE_TIME, HEIGHT_CM, HEIGHT_DATE_TIME,BMI,BMI_DATE_TIME))

    #Transform triage acuity, providing a value indicating it is unknown. 
    demographics$TRIAGE_ACUITY <- ifelse(demographics$TRIAGE_ACUITY=="", "Unknown", demographics$TRIAGE_ACUITY)


    #Filter out visits 'shorter than 0 minutes'.
    demographics <- demographics %>% select(-ADMIT_DATE) %>% filter(LENGTH.OF.STAY..IN.MINUTES. > 0 & LENGTH.OF.STAY..IN.MINUTES. < 24*60)

    print(paste("Visits within 24 hours", nrow(demographics)))


    #Identify previous visits (without admission) within 30 days, and visits with admission within 30 days; 

    #first, convert all visits to a timestamp (where the arrival time of the earliest visit is 0)
    first_stamp <- min(demographics$ED_ARRIVAL_TIME)
    demographics$arrival_timestamp <- difftime(demographics$ED_ARRIVAL_TIME, first_stamp, units="mins")

    #Get a proxy for NEDOCS- how many other patients are there when a patient arrives?
    demographics$departure_timestamp <- difftime(demographics$ED_CHECKOUT_TIME, first_stamp, units="mins")
    demographics$ord_num_patients_at_arrival<- purrr::map_vec(demographics$arrival_timestamp, function(x) get_patients_present_at_arrival(x, demographics$arrival_timestamp, demographics$departure_timestamp), .progress = TRUE)


    #Now work out how many visits there are within 30 days that don't result in admission
    #We only need to do this for patients who visit more than once
    prior_visits <- demographics %>% group_by(MRN) %>% mutate(num_total_visits = n()) %>% ungroup() %>% filter(num_total_visits > 1)


    #Count prior visits with and without admission.
    prior_visits <- prior_visits %>%
        arrange(MRN, arrival_timestamp) %>%
        group_by(MRN) %>% #Only look within visits by same patient
        mutate(
            num_prior_visits_without_admission = purrr::map_dbl(row_number(), function(i) {
            current_time <- arrival_timestamp[i]
            current_csn <- CSN[i]
            thirty_days_ago <- current_time - 24*60*30
            
            sum(arrival_timestamp < current_time & #Look for other visits from the same patient within 30 days
                arrival_timestamp > thirty_days_ago & 
                CSN != current_csn & 
                is_admitted == 0, 
                na.rm = TRUE)
            }),
            num_prior_admissions = purrr::map_dbl(row_number(), function(i) {
            current_time <- arrival_timestamp[i]
            current_csn <- CSN[i]
            thirty_days_ago <- current_time - 24*60*30
            
            sum(arrival_timestamp < current_time & 
                arrival_timestamp > thirty_days_ago & 
                CSN != current_csn & 
                is_admitted == 1, 
                na.rm = TRUE)
            })
        ) %>%
        ungroup()

    #Link back to the original dataframe, and fill in 0s for those with no entry
    prior_visits <- prior_visits %>% select(CSN, num_prior_admissions, num_prior_visits_without_admission)
    write.csv(prior_visits, "prior_visits.csv")
    demographics <- left_join(demographics, prior_visits, by="CSN")

    demographics$num_prior_admissions <- ifelse(is.na(demographics$num_prior_admissions), 0, demographics$num_prior_admissions)
    demographics$num_prior_visits_without_admission <- ifelse(is.na(demographics$num_prior_visits_without_admission), 0, demographics$num_prior_visits_without_admission)


    
    #Finally, clean names and return the file.
    demographics <- demographics %>% clean_names()
    write.csv(demographics, paste0(save_filepath, "preprocessed_demographics_for_epi.csv"))

}


#Next, get PCI scores, and COUNT the number of labs + meds given to each visit. Vitals will be linked later.
attach_pci_score <- function() {


    labs <- read.csv("labs.csv") 
    #Link and attach CSNs to radiology
    radiology <- read.csv("radiology_15Apr.csv")
    demographics <- read.csv("demographics_8May.csv")
    info_to_csn_df <- demographics %>% select(NAME, DOB, ED_ARRIVAL_TIME, CSN) %>% distinct(NAME, DOB, ED_ARRIVAL_TIME, CSN)
    radiology <- radiology %>% inner_join(info_to_csn_df, by=c("NAME", "DOB", "ED_ARRIVAL_TIME")) #now radiology has CSNs


    medications <- read.csv("medications.csv")
    diagnoses <- read.csv("diagnoses.csv")
    print("data loaded before stays")


    stays <- read.csv(paste0(save_filepath, "preprocessed_demographics_for_epi.csv"))
    print("data loaded")

    #Rename the 'prior admissions'/'visits' column.
    for (var in c("num_prior_admissions", "num_prior_visits_without_admission")) {
                            stays[[paste0("ord_", var)]] <- stays[[var]]
    }

    #Make and save a table of labs.
    tab <- table(toupper(labs$LAB))
    write.csv(tab[order(tab, decreasing = TRUE)], paste0(save_filepath, "labs_by_freq.csv"))

    #Filter out labs which occurred after a patient was discharged from the ED, and count the labs associated with each visit.
    labs <- labs %>% group_by(CSN) %>% filter(ymd_hms(LAB_START_TIME, truncated=3) < ymd_hms(ED_DISCHARGE_TIME, truncated=3)) %>% summarise(num_labs_over_visit = n())
    stays <- left_join(stays, labs, by=c("csn"="CSN"))
    stays$num_labs_over_visit <- ifelse(is.na(stays$num_labs_over_visit), 0, stays$num_labs_over_visit)


    #Now count the number of radiological tests, assuming EVENT END DT TM is when test is ordered
    radiology <- radiology %>% select(CSN, RADIOLOGY_TEST, ED_ARRIVAL_TIME, EVENT_END_DT_TM, ED_CHECKOUT_TIME)
    radiology <- radiology %>% group_by(CSN)  %>% 
            filter(ymd_hms(EVENT_END_DT_TM, truncated=3) < ymd_hms(ED_CHECKOUT_TIME, truncated=3), ymd_hms(EVENT_END_DT_TM, truncated=3) >= ymd_hms(ED_ARRIVAL_TIME, truncated=3)) %>% 
            summarise(num_tests_over_visit= n())
    stays <- left_join(stays, radiology, by=c("csn"="CSN"))
    stays$num_tests_over_visit <- ifelse(is.na(stays$num_tests_over_visit), 0, stays$num_tests_over_visit)

    #Save lists of the most common medication types/routes.
    tab <- table(toupper(medications$ROUTE))
    write.csv(tab[order(tab, decreasing = TRUE)], paste0(save_filepath, "routes_by_freq.csv"))

    tab <- table(toupper(medications$MEDICATION))
    write.csv(tab[order(tab, decreasing = TRUE)], paste0(save_filepath, "meds_by_freq.csv"))

    print("start counting medications")
    #Count 'mild painkillers', IV meds, and other meds. (Again, this is not the distinction that will be used in later analyses; this is to get a sense of the field.)
    medications <- medications %>% group_by(CSN)  %>% filter(ymd_hms(MED_DATE_TIME, truncated=3) < ymd_hms(ED_CHECKOUT_TIME, truncated=3)) %>%
            summarise(num_mild_painkillers = sum(grepl("ACETAMENOPHEN", MEDICATION) | grepl("IBUPROFEN", MEDICATION)), num_other_meds=n()-num_mild_painkillers, num_iv_meds = sum(ROUTE=="IV" | ROUTE=="intravenous"))
    
    print("medications counted")
    stays <- left_join(stays, medications, by=c("csn"="CSN"))
    stays$num_iv_meds <- ifelse(is.na(stays$num_iv_meds), 0, stays$num_iv_meds)
    stays$num_mild_painkillers <- ifelse(is.na(stays$num_mild_painkillers), 0, stays$num_mild_painkillers)
    stays$num_other_meds <- ifelse(is.na(stays$num_other_meds), 0, stays$num_other_meds)

    #Link PEDIATRIC COMORBIDITY SCORE.

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

    #Now link conditions to point scores.
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
    print("codes defined")

    #Get the conditions associated with each patient
    unique_codes <- unique(diagnoses$DIAGNOSIS_CODE)
    corresponding_conditions <- purrr::map_vec(unique_codes, ~find_condition(.x, paed_icd_code_dict), .progress = TRUE)
    codes_with_conditions <- which(!is.na(corresponding_conditions))

    unique_codes <- unique_codes[codes_with_conditions]
    corresponding_conditions <- corresponding_conditions[codes_with_conditions]
    names(corresponding_conditions) <- unique_codes
    diagnoses$condition <- purrr::map_vec(diagnoses$DIAGNOSIS_CODE, ~get_condition(.x, unique_codes, corresponding_conditions), .progress = TRUE)

    #Find unique conditions associated with each patient, and pull the PCI score.
    diagnoses <- diagnoses %>% filter(!is.na(condition)) %>% group_by(CSN) %>% distinct(condition, .keep_all = TRUE)
    diagnoses$score <- purrr::map_vec(diagnoses$condition, ~get_score(.x, paed_weight_dict), .progress=TRUE)
    diagnoses <- diagnoses %>% group_by(CSN) %>% summarise(PCI=sum(score))


    #Link back to the original dataframe
    stays <- left_join(stays, diagnoses, by=c("csn"="CSN")) 
    stays$PCI <- ifelse(is.na(stays$PCI), 0, stays$PCI) #the absence of legible conditions is a 0

    #transform but keep continuous
    stays$ord_PCI <- stays$PCI
    stays$PCI <- ifelse(stays$PCI == 0, "none", ifelse(stays$PCI==1, "one", ifelse(stays$PCI>1 & stays$PCI <5, "two_to_four", "five_or_more")))


    #Flag ordinal data
    for (var in c("num_mild_painkillers", "num_other_meds",
                    "num_iv_meds", "num_tests_over_visit", "num_labs_over_visit")) {
                            stays[[paste0("ord_", var)]] <- stays[[var]]
                    }
                            


    #Finally, get temporal variables
    stays$year_of_arrival <- lubridate::year(stays$ed_arrival_time)
    stays$month_of_arrival <- lubridate::month(stays$ed_arrival_time)
    stays$hour_of_arrival <- lubridate::hour(ymd_hms(stays$ed_arrival_time))
    stays$day_of_arrival <- lubridate::wday(stays$ed_arrival_time)

    #Filter out 66 visits with no legible hour
    stays <- stays %>% filter(!is.na(hour_of_arrival))


    #Convert prior visits and admissions to 0, 1, multiple
    stays$num_prior_admissions <- ifelse(stays$num_prior_admissions==0, "none", ifelse(stays$num_prior_admissions==1, "one", "two_or_more"))
    stays$num_prior_visits_without_admission <- ifelse(stays$num_prior_visits_without_admission==0, "none", ifelse(stays$num_prior_visits_without_admission==1, "one", "two_or_more"))

    #Convert labs and tests to 0 (none) and 1 (any)
    stays$mild_painkillers <- ifelse(stays$num_mild_painkillers>0, 1, 0)
    stays$other_meds <- ifelse(stays$num_other_meds>0, 1, 0)
    stays$any_iv_meds <- ifelse(stays$num_iv_meds>0, 1, 0)
    stays$any_labs <- ifelse(stays$num_labs_over_visit>0, 1, 0)
    stays$any_tests <- ifelse(stays$num_tests_over_visit>0, 1, 0)

    first_arrival <- min(as.POSIXct(stays$ed_arrival_time))
    stays$timestamp <- difftime(as.POSIXct(stays$ed_arrival_time), first_arrival, units="mins")

    #Drop unnecessary columns
    stays <- stays %>% select(-c(ed_arrival_time, triage_start_dt_tm,
                    num_mild_painkillers, num_other_meds, num_iv_meds, num_tests_over_visit))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       #num_labs_over_visit

    #filter out one visit which appears to be from 2012?
    stays <- stays %>% filter(!(year_of_arrival==2012))

    write.csv(stays, paste0(save_filepath, "linked_stays_for_epi.csv"))
}

#Pull baseline factors from the above files and re-link the original chief complaint.
get_baseline_factors <- function() {
    baseline_factors <- read.csv(paste0(save_filepath, "preprocessed_demographics_for_epi.csv"))
    linked_stays_for_epi <- read.csv(paste0(save_filepath, "linked_stays_for_epi.csv"))
    linked_stays_for_epi <- linked_stays_for_epi %>% select(all_of(c("csn", "hour_of_arrival", "day_of_arrival", "month_of_arrival", "year_of_arrival"))) 
    
    baseline_factors <- inner_join(baseline_factors, linked_stays_for_epi, by="csn")
    baseline_factors <- baseline_factors %>% select(-c(admit_source, admitted_service, mrn, ed_arrival_time, ed_checkout_time, length_of_stay_in_minutes, triage_start_dt_tm, triage_complete_dt_tm,
            age_in_years, sdi_quantile, miles_travelled, weight, departure_timestamp)) #keep both age group and continuous age, also arrival timestamp; discard all other 'categorised' versions of continuous variables

    write.csv(baseline_factors, paste0(save_filepath, "baseline-factors-for-patient-journey.csv"))
    baseline_factors <- baseline_factors  %>% rename(cts_age_in_days=age_in_days) #so that age_in_days is correctly labelled


    #CONNECT original complaints, disconnect complaint columns
    demographics <- read.csv("demographics_8May.csv")
    complaint_cols <- colnames(baseline_factors)[startsWith(colnames(baseline_factors), "complaint_")]
    baseline_factors <- baseline_factors %>% select(-c(all_of(complaint_cols)))
    chief_complaints <- demographics %>% select(CSN, ED_COMPLAINT) %>% clean_names()
    baseline_factors <- inner_join(baseline_factors, chief_complaints, by="csn")
    baseline_factors$ed_complaint <- tolower(baseline_factors$ed_complaint)

    write.csv(baseline_factors, paste0(save_filepath, "baseline_factors_with_og_complaint.csv"))

    #Save a dataframe linking CSN to complaint
    df <- read.csv(paste0(save_filepath, "baseline_factors_with_og_complaint.csv")) %>% select(csn, ed_complaint)
    write.csv(df, paste0(save_filepath, "csn-to-complaint.csv"))
}



#Identify the 200 most common 'word stems' in chief complaints
tag_complaints <- function() {

    baseline_factors <- read.csv(paste0(save_filepath, "preprocessed_demographics_for_epi.csv"))
    linked_stays_for_epi <- read.csv(paste0(save_filepath, "linked_stays_for_epi.csv"))
    linked_stays_for_epi <- linked_stays_for_epi %>% select(all_of(c("csn", "hour_of_arrival", "day_of_arrival", "month_of_arrival", "year_of_arrival"))) 
    baseline_factors <- inner_join(baseline_factors, linked_stays_for_epi, by="csn")

    #Link OG complaint
    original_complaints <- read.csv(paste0(save_filepath, "csn-to-complaint.csv"))
    baseline_factors <- baseline_factors %>% select(-starts_with("complaint_")) %>% inner_join(original_complaints, by="csn") %>% select(-c(mrn, admit_source, X.x, X.y,
            admitted_service, admit_source, triage_complete_dt_tm, triage_start_dt_tm, age_in_years, miles_travelled, sdi_quantile, 
            weight, departure_timestamp))


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

    #expand common abbreviations
    for (abbrev in names(common_abbrev)) {
        pattern <- paste0("\\b", abbrev, "\\b")
        baseline_factors$ed_complaint <- str_replace_all(baseline_factors$ed_complaint, regex(pattern, ignore_case = TRUE), common_abbrev[[abbrev]]) #   
    }


    #Extract a list of word stems from each complaint.

    #Get each individual word
    all_words <- c()
    all_words <- baseline_factors$ed_complaint %>% purrr::map(~extract_words(.x, all_words), .progress=TRUE) %>% unlist()

    #Order by descending frequency of stems
    stems_by_frequency <- data.frame(stem=wordStem(all_words)) %>% group_by(stem) %>% summarise(count=n()) %>% arrange(desc(count))

    #Use top 200: 200th most common complaint appears in around 100 visits

    stems_to_use <- stems_by_frequency$stem[1:200]
    #Exclude nonsensical words
    stems_to_use <- stems_to_use[!(stems_to_use %in% c("", "g", "l", "0", "n", "000"))]

    #Link with CSN
    complaint_contains_matrix <- purrr::map_dfr(1:nrow(baseline_factors), ~convert_complaints(.x, baseline_factors, stems_to_use), .progress = TRUE) 
    ed_complaints_only <- baseline_factors %>% select(csn, ed_complaint)

    check_matrix <- complaint_contains_matrix %>% inner_join(ed_complaints_only, by="csn")
    
    baseline_factors <- baseline_factors %>% inner_join(complaint_contains_matrix, by="csn") %>% select(-ed_complaint)
    write.csv(baseline_factors, paste0(save_filepath, "baseline-factors-with-top-200-complaint-stems-tagged.csv"))

}

#Pull out key medications and the times they were administered.
link_key_meds <- function() {

    #Load data
    baseline_factors <- read.csv((paste0(save_filepath, "baseline-factors-with-top-200-complaint-stems-tagged.csv"))) %>% select(-length_of_stay_in_minutes)
    print("baseline factors")
    labs <- read.csv("labs.csv")
    print("labs")
    medications <- read.csv("medications.csv")
    print("meds")
    scores <- read.csv("scores.csv")
    print("scores")
    scores <- scores %>% filter(SCORE_VARIABLE=="NRS Generalized Pain Score")
    vitals <- read.csv("vitals.csv")
    print("vitals")

    #Link and attach CSNs to radiology
    radiology <- read.csv("radiology_15Apr.csv")
    print("rad loaded")
    demographics <- read.csv("demographics_8May.csv")
    info_to_csn_df <- demographics %>% select(NAME, DOB, ED_ARRIVAL_TIME, CSN) %>% distinct(NAME, DOB, ED_ARRIVAL_TIME, CSN)
    radiology <- radiology %>% inner_join(info_to_csn_df, by=c("NAME", "DOB", "ED_ARRIVAL_TIME")) #now radiology has CSNs



    print("loaded")
    #Compile a sequence of events that happens to each patient over their ED stays


    events <- data.frame(csn=character(0), event_type=character(0), event_name=character(0), event_time=character(0), lab_id=numeric(0))

    #First, add in arrival
    arrival_times <- baseline_factors %>% select(csn, ed_arrival_time) %>% rename(event_time=ed_arrival_time) %>% mutate(event_type="arrival", event_name="arrival", lab_id=NA, event_value=NA)
    events <- rbind(events, arrival_times) #this will end up being discarded in the end


    #Get admission times for patients admitted
    admission_times<- baseline_factors %>% filter(is_admitted==1) %>% select(csn, ed_checkout_time) %>% rename(event_time=ed_checkout_time) %>% mutate(event_type="endpoint", event_name="admission", lab_id=NA, event_value=NA)
    events <- rbind(events, admission_times)

    discharge_times<- baseline_factors %>% filter(is_discharged==1) %>% select(csn, ed_checkout_time) %>% rename(event_time=ed_checkout_time) %>% mutate(event_type="endpoint", event_name="discharge", lab_id=NA, event_value=NA)
    events <- rbind(events, discharge_times)

    departure_times<- baseline_factors %>% filter(is_discharged==0 & is_admitted==0) %>% select(csn, ed_checkout_time) %>% rename(event_time=ed_checkout_time) %>% mutate(event_type="endpoint", event_name="other departure", lab_id=NA, event_value=NA)
    events <- rbind(events, departure_times)


    #Get radiological tests
    rad_times <- radiology %>% select(CSN, EVENT_END_DT_TM, RADIOLOGY_TEST) %>% 
        rename(csn=CSN, event_time=EVENT_END_DT_TM, event_name=RADIOLOGY_TEST) %>% 
            mutate(event_type="test_ordered", lab_id=NA, event_value=NA)


    events <- rbind(events, rad_times)

    #Compile pain scores and vitals

    pain_score_times <- scores %>% select(CSN, SCORE_DT_TIME, SCORE) %>% rename(csn=CSN, event_time=SCORE_DT_TIME, event_value=SCORE) %>% mutate(event_type="vitals_taken", event_name="PAIN", lab_id=NA)
    events <- rbind(events, pain_score_times)

    vital_times <- vitals %>% select(CSN, VITAL_DATE_TIME, VITAL, VITAL_VALUE) %>% rename(csn=CSN, event_time=VITAL_DATE_TIME, event_name=VITAL, event_value=VITAL_VALUE) %>% mutate(event_type="vitals_taken", lab_id=NA)
    events <- rbind(events, vital_times)


    #Attach SPECIFIC medications indicated by clinical collaborators (or flagged by the model) as indicators of likely admission or discharge
    #We only count non-intravenous medications here
    corticosteroids <- c("DEXAMETHOSONE", "BUDESONIDE", "HYDROCORTISONE", "METHYLPREDNISOLONE", "PREDNISOLONE", "PREDNISONE", "TRIAMCINOLONE") 
    antihistamines <- c("CETIRIZINE")
    adrenaline <- c("EPINEPHRINE") #will include racepinephrine
    antiepileptic <- c("MIDAZOLAM", "PHENOBARBITAL", "BRIVARACETAM", "CLOBAZAM", "CLONAZEPAM", "DIAZEPAM", "GABAPENTIN", "LACOSAMIDE",
        "PHENYTOIN", "PREGABALIN", "PRIMIDONE", "RUFINAMIDE", "TOPIRAMATE", "VALPROIC ACID", "VIGABATRIN", "ZONISAMIDE") #taken from https://pmc.ncbi.nlm.nih.gov/articles/PMC6130739/
    pain_meds <- c("TRAMADOL", "OXYCODONE", "MORPHINE", "FENTANYL", "KETAMINE") #all but the last are opioids
    antiemetic <- c("ONDANSETRON")
    oral_ibuprofen <- c("ORAL IBUPROFEN")
    topical_let <- c("TOPICAL LET")

    
    #Also note whether anything is administered intravenously
    meds_to_count <- c(corticosteroids, antihistamines, adrenaline, antiepileptic, pain_meds, antiemetic, "IV", oral_ibuprofen, topical_let)

    #Also note the presence of oral ibuprofen
    medications$MEDICATION <- toupper(medications$MEDICATION)
    medications$MEDICATION <- ifelse(medications$ROUTE == "IV" | medications$ROUTE == "INTRAVENOUS", "IV", medications$MEDICATION)
    medications$MEDICATION <- ifelse(grepl("IBUPROFEN", medications$MEDICATION) & medications$ROUTE=="PO", "ORAL IBUPROFEN", medications$MEDICATION)
    medications$MEDICATION <- ifelse(grepl("LET", medications$MEDICATION) & medications$ROUTE=="TOP", "TOPICAL LET", medications$MEDICATION)


    med_times <- medications %>% select(CSN, MED_DATE_TIME, MEDICATION) %>% filter(grepl(paste(meds_to_count, collapse="|"), MEDICATION)) %>% 
            rename(csn=CSN, event_time=MED_DATE_TIME, event_name=MEDICATION) %>% mutate(event_type="medication_given", lab_id=NA, event_value=NA)
    events <- rbind(events, med_times)


    #Clean names
    events <- events %>% mutate(cleaned_event_name = clean_value_names(event_name)) %>% 
        select(-event_name) %>%
        rename(event_name=cleaned_event_name)

    #Now convert timestamps to "minutes since arrival"
    #Make sure we get only events taken before checkout time
    arrival_timestamps <- arrival_times %>% select(csn, event_time) %>% rename(arrival_time=event_time)

    events <- inner_join(events, arrival_timestamps, by="csn")
    events$event_time <- difftime(ymd_hms(events$event_time, truncated = 3), ymd_hms(events$arrival_time, truncated = 3), units="mins") #allowing for missing seconds
    endpoints <- events %>% filter(event_type == "endpoint") %>% select(csn, event_time) %>% rename(endpoint=event_time)
    events <- inner_join(events, endpoints, by="csn") %>% filter(event_time <= endpoint) %>% select(-endpoint)

    events <- events %>% filter(!(event_type=="arrival")) %>% group_by(csn) %>% arrange(event_time, .by_group = TRUE)

    write.csv(events, paste0(save_filepath, "focused-patient-journey-events.csv"))

}

#Preprocess main demographic file.
#link_demographics()
#Attach PCI score.
gc()
#attach_pci_score()
#Pull out key baseline factors and attach original ED complaint.
#get_baseline_factors()
#Pull out 200 most common stems from ED complaint.
#tag_complaints()
#Link key meds, and the time at which they were given to patients
link_key_meds()


