
library(lubridate)
library(broom)
library(stringr)
library(patchwork)
library(pROC)

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
library(yardstick)
library(comorbidity)
library(readr)
library(testit)


#Get odds ratios for triage acuity and admissions.

# setwd("/Volumes/chip-lacava/Groups/BCH-ED/reprocessing/")


setwd("/Users/helenacoggan/Documents/OSXLAP12460/BCH-data-2025/temp-dict/reprocessing/")



#Set seed.
set.seed(2115) #Zip code of Boston Children's Hospital.


#Binarise baseline factors.
convert_baseline_factors <- function(visits) {

    converted_factors <- data.frame(visits)

    #Replace all NAs with 'unknown'.
    for (var in colnames(converted_factors)) {
        na_idx <- which(is.na(converted_factors[[var]]))
        if(length(na_idx)>0) {
            converted_factors[na_idx, var] <- "unknown"
        }
    }

    #Convert categorical variables.
    factors_to_convert <- c("sex", "race", "ed_arrival_mode", "language", "triage_acuity", 
        "age_group", "insurance", "state_of_origin",
        "year_of_arrival", "season", "time_of_day", "triage_pain",
        "overall_mean_pain", "overall_max_pain", "overall_min_pain") #pain-unknown is already marked, so default here is none
    refs <- c("M", "Non-Hispanic White", "Walk in", "English", "unknown",
             "fifteen_and_older", "Private",
            "in-state", 2019, "winter", "morning", "none", "none", "none", "none")
    for (i in 1:length(factors_to_convert)) {
        var <- factors_to_convert[i]
        ref <- refs[i]

        assert(ref %in% unique(converted_factors[[var]]))
        for (val in unique(converted_factors[[var]])) {
            if (!is.na(val) & !(val==ref)) {
                name <- paste0(var, "_", val)
                converted_factors[[name]] <- ifelse(converted_factors[[var]]==val, 1, 0)
                assert(!any(is.na(converted_factors[[name]])))
            }
        }
    }
    #Add an extra step where unknown continuous factors are replaced with their cohort-wide mean.
    #Assume all patients with no recorded weight have normal weight, as weight is often not recorded.
    converted_factors$weight_unknown <- ifelse(converted_factors$weight=="unknown", 1, 0)
    converted_factors$weight <- ifelse(converted_factors$weight=="unknown", 0, as.numeric(converted_factors$weight))

    #For the other two it should be a specific cohort-wide mean.
    cts_factors_to_convert <- c("sdi_score", "miles_travelled")
    for (var in cts_factors_to_convert) {
        mean <- mean(as.numeric(converted_factors[[var]]), na.rm=TRUE)
        converted_factors[[paste0(var, "_unknown")]] <- ifelse(converted_factors[[var]]=="unknown", 1, 0)
        converted_factors[[var]] <- ifelse(converted_factors[[var]]=="unknown", mean, as.numeric(converted_factors[[var]]))

    }
    


    #Discard converted factors (except the continuous ones),
    converted_factors <- converted_factors %>% select(-c(all_of(c(factors_to_convert)))) %>% clean_names()

    return(converted_factors)
}

#Calculate propensity scores before and after triage, and for partial and full adjustment (after triage).

get_propensity_scores <- function(visits, exposures, covariates, type) {

    #Copy the dataframe and mark the reference category. (This requires us to have all exposures in a category listed in the exposures column.)
    df <- data.frame(visits)
    df$ref <- ifelse(rowSums(df[exposures])==0, 1, 0)

    #Exclude exposures from covariates.
    covariates <- covariates[!(covariates %in% exposures)]

    #Exclude sex-related covariates (CHECK WHICH THESE ARE.)
    if(length(covariates) > 0 & any(startsWith(exposures, "sex"))) {
        covariates <- covariates[!startsWith(covariates, "sex") & 
            !endsWith(covariates, "gynecologic") & 
            !endsWith(covariates, "male_genital") &
            !endsWith(covariates, "menstrual_disorders")]
    }


    #Include only the exposures, references, and the variables we're controlling for.
    df <- df %>% select(all_of(c("csn", exposures, "ref", covariates)))


    #Initialise propensity score dataframe.
    propensity_scores <- data.frame()
    for (exposure in exposures) {

        #Filter to only exposure and reference (e.g. Black and White).
        relevant_population <- df %>% filter(.data[[exposure]]==1 | ref==1)

        #Shuffle dataframe.
        shuffled_order <- sample(1:nrow(relevant_population), size=nrow(relevant_population), replace=FALSE)
        relevant_population <- relevant_population[shuffled_order,]


        #Convert factors to numerics.
        csns <- relevant_population$csn

        #Calculate propensity scores if type is not 'unadjusted'.
        if (!(type=="unadjusted")) {
            Xs <- relevant_population %>% select(all_of(covariates)) 
            # for (col in colnames(Xs)) {
            #     print(col)
            #     assert(!any(is.na(as.numeric(Xs[[col]]))))
            # }
            Xs <- as.matrix(as.data.frame(sapply(Xs, as.numeric)))
            ys <- factor(relevant_population[[exposure]], levels=c(0, 1))


            #Fit propensity score with GLMNET.
            ps_model <- cv.glmnet(Xs, ys, family = binomial(link="logit"), nfolds=5, trace.it = 1) #use 5fold CV
            glm_ps <- as.vector(predict(ps_model, type="response", newx=Xs, s=ps_model$lambda.1se)) #use the value of s that gives the lowest error

            #Add to overall dataframe.
            subscores <- data.frame(csn=csns, treatment=exposure, propensity_score=glm_ps, type=type)
        } else { #If type is unadjusted, just assign everyone a propensity score of 0.5, as we'll match them at random.
            subscores <- data.frame(csn=csns, treatment=exposure, propensity_score=0.5, type=type)
        }

        #Add to overall dataframe.
        propensity_scores <- rbind(propensity_scores, subscores)

    }
    return(propensity_scores)

}


#Assemble all propensity score categories.
assemble_propensity_scores <- function(visits, races) {

    #Define all sets of variables.
    adjustment_variables <- list(
        "unadjusted" = c(),
        "pre-triage" = vars_known_at_triage,
        "partial" = partial_adjustment_vars,
        "full" = full_adjustment_vars
    )

    all_propensity_scores <- data.frame()

    #Get propensity scores for race and sex-based categories.
    for (type in names(adjustment_variables)) {
        covariates <- adjustment_variables[[type]]



        race_scores <- get_propensity_scores(visits, races, covariates, type)
        sex_scores <- get_propensity_scores(visits, c("sex_f"), covariates, type)

        all_propensity_scores <- rbind(all_propensity_scores, race_scores, sex_scores)

    }

    return(all_propensity_scores)
}


#Create an over-arching function which takes a (filtered) dataset, an exposure, an outcome, 
#and a set of variables to adjust for (i.e. a propensity score keyword), and then calculates odds ratios.
get_odds_ratios <- function(filtered_visits, propensity_scores, 
    outcome_vars, exposure, ref, covariates, ps_type, threshold=0.1, savetag="") {

 

    #Define framework.
    odds_ratios <- data.frame(exposure=character(0), 
        outcome_var=character(0), 
        estimate=numeric(0),
        lower_conf=numeric(0),
        upper_conf=numeric(0), 
        standard_error=numeric(0),
        pval=numeric(0), 
        num_treated=numeric(0),
        num_control=numeric(0),
        all_treated=numeric(0),
        treated_outcome_frac=numeric(0),
        control_outcome_frac=numeric(0))

    index <- 1
    
    df <- data.frame(filtered_visits) 


    #Exclude exposures from covariates.
    #Also exclude all sex-based covariates (e.g. those related to pregnancy) 
    #when the exposure is sex.
    if(length(covariates) > 0  & startsWith(exposure, "race")) {
        covariates <- covariates[!startsWith(covariates, "race")]
    }
    if(length(covariates) > 0 & startsWith(exposure, "sex")) {
        covariates <- covariates[!startsWith(covariates, "sex") &
            !endsWith(covariates, "gynecologic") & 
            !endsWith(covariates, "male_genital") &
            !endsWith(covariates, "menstrual_disorders")]
    }

    #If examining race-specific effect of languages, exclude all relevant variables:
    if (startsWith(exposure, "language")) {
        covariates <- covariates[!(startsWith(covariates, "language"))]
    }

    #Choose the relevant propensity score and thus filter out all patients not in exposure and reference (as they will have NAs)
    ps <- propensity_scores %>% filter(treatment==exposure, type==ps_type)
    df <- df %>% inner_join(ps, by="csn") %>% filter(!is.na(propensity_score))

    for (col in colnames(df)) {
        if (any(is.na(is.numeric(df[[col]])) | is.null(is.numeric(df[[col]])))) {
            print(col)
        }
    }


    #Determine caliper threshold.
    if(ps_type=="unadjusted") {
        threshold_value <- NULL
        #Define formula for balancing
        ps_formula <- as.formula(paste(exposure, "~ propensity_score"))
    } else {
        threshold_value <- threshold*sd(df$propensity_score) #We're not using a caliper threshold if not adjusting.
        #Define formula for balancing
        ps_formula <- as.formula(paste(exposure, "~", paste(covariates, collapse = " + ")))
    }


    #Perform matching.
    tryCatch({

        match_object <- matchit(ps_formula, data=df, method = 'nearest', caliper=threshold_value, distance = df$propensity_score)

 

        if (!(ps_type=="unadjusted")) {
            #Plot only the top 25 variables in a Love plot.
            bal_stats <- bal.tab(match_object)
            ranked_variables <- rownames(bal_stats$Balance[order(abs(bal_stats$Balance$Diff.Adj), decreasing = TRUE), ])
            num_to_take <- min(25, length(ranked_variables))
            ordered_vars <- ranked_variables[1:num_to_take]
            ordered_vars <- ordered_vars[!(ordered_vars=="distance")]
            if ("cci_of_visit_2" %in% ordered_vars) {
                ordered_vars[which(ordered_vars=="cci_of_visit_2")] <- "cci_of_visit" #Sometimes baltab will change the name of this variable, so we change it back.
            }
            p <- love.plot(as.formula(paste(exposure, "~", paste(ordered_vars, collapse = " + "))), data = df, weights = match_object$weights, stats = c("mean.diffs"), 
                    thresholds = c(m = .1, v = 2), abs = TRUE, 
                    binary = "std",
                    continuous = "std", 
                    var.order = "adjusted",
                    stars = "std",
                    s.d.denom = "pooled",
                    ) + theme(axis.text.y = element_text(size = 8))
            ggsave(paste0("outputs/caliper_", threshold, "/", "balance_plots/", savetag, "_", exposure, "_", ps_type, "_balance_plot.pdf"), plot=p)
        }

        #Make two data tables for treated and untreated data (exposure and reference), rejoin on pairing index.
        matched <- as.data.table(match_data(match_object))
        treated_indices <- which(matched[[exposure]]==1)
        control_indices <- which(matched[[exposure]]==0)

        #Calculate the number of patients we could have had before matching.
        all_treated <- sum(df[[exposure]]==1)


        #Run McNemar's test on each outcome.
        for (outcome_var in outcome_vars) {
            assert(!any(is.na(df[[outcome_var]])))
            treated_df <- data.frame(pair_index=matched[treated_indices, subclass], treatment_outcome=matched[[outcome_var]][treated_indices])
            control_df <- data.frame(pair_index=matched[control_indices, subclass], control_outcome=matched[[outcome_var]][control_indices])
            overall_table <- inner_join(treated_df, control_df, by="pair_index") 

            a <- sum(overall_table$treatment_outcome==1 & overall_table$control_outcome==1)
            c <- sum(overall_table$treatment_outcome==1 & overall_table$control_outcome==0)
            b <- sum(overall_table$treatment_outcome==0 & overall_table$control_outcome==1)
            d <- sum(overall_table$treatment_outcome==0 & overall_table$control_outcome==0)


            #Conduct test.
            matrix <- matrix(c(a, b, c, d), 2, 2)
            test <- mcnemar.exact(matrix)


            #Calculate raw admission rates.
            treatment_outcome_rate <- sum(df[[exposure]]==1 & df[[outcome_var]]==1)/sum(df[[exposure]]==1)
            control_outcome_rate <- sum(df[[exposure]]==0 & df[[outcome_var]]==1)/sum(df[[exposure]]==0)

            lower_conf <- as.numeric(test$conf.int[1])
            upper_conf <- as.numeric(test$conf.int[2])

            standard_error <- (upper_conf - lower_conf)/(3.92)

            #Save in table.
            mcnemar_result <- c(exposure, outcome_var, as.numeric(test$estimate), lower_conf, 
                    upper_conf, standard_error, as.numeric(test$p.val),
                    length(treated_indices), length(control_indices), all_treated,
                    treatment_outcome_rate, control_outcome_rate)

            odds_ratios[index,] <- mcnemar_result
            index <- index + 1

        } }, error = function(e) { print(paste("Matching failed: ", e$message))})
        
    return(odds_ratios)

}


#Calculate overall admission (for overall, ICU and observation) separately at 3 levels of adjustment.
#MHC correct across all levels.
get_admission_odds_ratios <- function(visits, propensity_scores, races, threshold=0.1, savetag="") {

    #Define all exposures and references.
    exposures <- c("sex_f", races)
    refs <- c("sex_m", rep("race_non_hispanic_white", length(races)))

    #Define outcomes.
    outcome_vars <- c("is_admitted")

    #Define all sets of variables.
    adjustment_variables <- list(
        "unadjusted" = c(),
        "partial" = partial_adjustment_vars,
        "full" = full_adjustment_vars
    )

    all_results <- data.frame()

    for (type in names(adjustment_variables)) {
        covariates <- adjustment_variables[[type]]
        sub_df <- data.frame()
        for (idx in 1:length(exposures)) {
            exposure <- exposures[idx]
            ref <- refs[idx]
            ors <- get_odds_ratios(visits, propensity_scores, outcome_vars,
                exposure, ref, covariates, type, threshold=threshold, savetag=savetag)
                sub_df <- rbind(sub_df, ors)
        }
        sub_df$adjustment_level <- type
        all_results <- rbind(all_results, sub_df)
    }

    return(all_results)
}





#Create a database of all 'gateway outcomes' (ED LOS quartiles, triage, receiving each of top 20 medications, labs, tests, etc).
#This returns a database where each outcome has become a binary variable.
create_outcomes_database <- function(visits, quartiles, quartile_names) {
    #Triage categories are already binary variables.
    #Vitals are already categorised as binary variables.
    #Insurance status is already categorised as a binary variable.


    #Time-to-event quartiles need to be returned as binary variables, i.e., each visit needs to have 4s variable indicating whether it is in the 1st, 2nd, 3rd or 4th ED LOS quartile.
    for (name in names(quartiles)) {
        var_name <- quartile_names[[name]]
        visits[[var_name]] <- as.numeric(visits[[var_name]])
        qs <- quartiles[[name]] 
        if (name %in% c("admit_time_quartiles", "request_to_admission_time_quartiles")) { #Mark quartiles only for admitted patients.
            for (i in 1:4) {
                if (i==4) { #Use the upper limit only for the last quartile.
                    visits[[paste0(name, "_", i)]] <- ifelse(visits[[var_name]] >= qs[i] & visits[[var_name]] <= qs[i+1] & visits$is_admitted==1, 1, 0)
                } else {
                    visits[[paste0(name, "_", i)]] <- ifelse(visits[[var_name]] >= qs[i] & visits[[var_name]] < qs[i+1] & visits$is_admitted==1, 1, 0)
                }
            }
        } else { #Mark quartiles for all patients.
            for (i in 1:4) {
                if (i==4) {
                    visits[[paste0(name, "_", i)]] <- ifelse(visits[[var_name]] >= qs[i] & visits[[var_name]] <= qs[i+1], 1, 0)
                } else {
                    visits[[paste0(name, "_", i)]] <- ifelse(visits[[var_name]] >= qs[i] & visits[[var_name]] < qs[i+1], 1, 0)
                }
            } 
        }
    }

    return(visits)

}


#Calculate odds ratios upstream and downstream of these outcomes using overall function.

#Broadly speaking there are three classes of odds ratio, each of which we will calculate and MHC-correct (later, when plotting.)

#There is one set of ORs where the cohort is 'all visits' and the outcome is intervention: all the 'explicit' interventions (begin with 'received'), plus triage score, ED-LOS, and rooming time quartiles.
get_upstream_odds_ratios <- function(visits, propensity_scores, threshold=0.1) {

    #Pull all variable names from the dataframe.
    vars <- colnames(visits)

    #Set aside ECG and test outcomes, where we want to test only for outcomes pre-EPIC.
    pre_epic_outcomes <- c("received_any_tests")

    #Define outcomes where we want to adjust for all variables received across the visit.
    non_triage_outcome_vars <- vars[(startsWith(vars, "received_") |
        startsWith(vars, "ed_los_")) & !(vars %in% c(pre_epic_outcomes, "received_any_ecgs"))]

    
    #Triage acuity scores should be adjusted only for variables known at triage.
    triage_outcome_vars <- vars[startsWith(vars, "triage_acuity_")]
    
    #Define all exposures and references.
    exposures <- c("sex_f", races)
    refs <- c("sex_m", rep("race_non_hispanic_white", length(races)))


    #Calculate results for all exposures and outcomes. 
    all_results <- data.frame()
    for (idx in 1:length(exposures)) {
        exposure <- exposures[idx]
        ref <- refs[idx]

        #Now calculate ORs for ECGs and tests, where we're interested only in the pre-EPIC cohort.
        ors <- get_odds_ratios(filter(visits, pre_epic==1), propensity_scores, pre_epic_outcomes,
            exposure, ref, full_adjustment_vars, "full", threshold=threshold, savetag="upstream-full")
        all_results <- rbind(all_results, ors)

        
        #First calculate ORs for interventions (EXCEPT ECG and tests) and ED-LOS, where we adjust for all variables known across visit.
        ors <- get_odds_ratios(visits, propensity_scores, non_triage_outcome_vars,
            exposure, ref, full_adjustment_vars, "full", threshold=threshold, savetag="upstream-full")
        all_results <- rbind(all_results, ors)

        #Now calculate ORs for triage acuity, where we adjust for all variables known at triage.
        ors <- get_odds_ratios(visits, propensity_scores, triage_outcome_vars,
            exposure, ref, partial_adjustment_vars, "pre-triage", threshold=threshold, savetag="upstream-partial")
        all_results <- rbind(all_results, ors)

    }

    #Don't MHC correct for now.

    return(all_results)
}


#There is another set of ORs where the cohorts are filtered (i.e. all those who received a certain intervention/triage score/ED-LOS, all those with a certain form of insurance,
#all those with a certain complaint (of top 10), all those with certain overall vital categories (take worst category as maximum), age, cci-before-visit, num-previous visits, ambulance), and the outcome is admission.
get_downstream_odds_ratios <- function(visits, propensity_scores, threshold=0.1) {

    #Pull all variable names from the dataframe.
    vars <- colnames(visits)


    #Pull the filtration characteristics which are already coded in the database.
    filtration_characteristics <- c(vars[(startsWith(vars, "received_") & !(endsWith(vars, "ecgs"))) |
         startsWith(vars, "ed_los_") | startsWith(vars, "triage_acuity_") | startsWith(vars, "pre_diagnosis_") |
         startsWith(vars, "current_diagnosis_")])

    
    #Include all complaint categories EXCEPT those related to sex
    complaints <- vars[startsWith(vars, "complaint_") & !startsWith(vars, "sex") & 
            !endsWith(vars, "gynecologic") & 
            !endsWith(vars, "male_genital") &
            !endsWith(vars, "menstrual_disorders")]
    filtration_characteristics <- c(filtration_characteristics, complaints)



    #Now add back 'reference' categories for other characteristics.
    characteristics_to_add <- c("age_group", "ed_arrival_mode", "insurance", "language")
    refs_to_add <- c("fifteen_and_older", "walk-in", "private", "english")

    # characteristics_to_add <- c("language")
    # refs_to_add <- c("english")

    for (i in 1:length(characteristics_to_add)) {
        relevant_vars <- colnames(visits)[startsWith(colnames(visits), characteristics_to_add[i])]
        visits[[paste0(characteristics_to_add[i], "_", refs_to_add[i])]] <- ifelse(rowSums(visits[relevant_vars])==0, 1, 0)
        filtration_characteristics <- c(filtration_characteristics, paste0(characteristics_to_add[i], "_", refs_to_add[i]), relevant_vars)
    }





    #Divide miles travelled, SDI score into quartiles.
    names <- c("sdi_score", "miles_travelled", "pseudo_nedocs")
    for (name in names) {
        qs <- quantile(visits[[name]], probs=c(0, 0.25, 0.5, 0.75, 1), na.rm=TRUE)
        for (i in 1:4) {
            char_name <- paste0(name, "_", i)
            if (i==4) {
                visits[[char_name]] <- ifelse(visits[[name]] >= qs[i] & visits[[name]] <= qs[i+1], 1, 0)
            } else {
                visits[[char_name]] <- ifelse(visits[[name]] >= qs[i] & visits[[name]] < qs[i+1], 1, 0)
            }
            filtration_characteristics <- c(filtration_characteristics, char_name)
        } 
    }

    #Record previous admissions
    visits$ord_num_previous_admissions <- case_when(
        visits$num_previous_admissions == 0 ~ "0",
        visits$num_previous_admissions == 1 ~ "1",
        visits$num_previous_admissions >= 2 ~ "2_or_more"
    )

    #Record previous discharges
    visits$ord_num_previous_visits_without_admission <- case_when(
        visits$num_previous_visits_without_admission == 0 ~ "0",
        visits$num_previous_visits_without_admission == 1 ~ "1",
        visits$num_previous_visits_without_admission >= 2 ~ "2_or_more"
    )

    #Record deviance of HR, RR, weight.
    for (name in c("weight", paste0("overall_mean_", c("hr", "rr")))) {
        visits[[paste0("ord_", name)]] <- case_when(
            visits[[paste0(name, "_unknown")]] == 1 ~ "unknown",
            visits[[name]] < 1 ~ "less_than_1_sd",
            visits[[name]] < 2 ~ "1_to_2_sd",
            visits[[name]] >= 2 ~ "2_or_more_sd"
        )
    }

    #SPO2, fever (>0), hypotension (>0)
    #Categorise blood oxygen (90 is normal for kids).
    visits$ord_spo2 <- case_when(
        visits$overall_mean_sp_o2_unknown == 1 ~ "unknown",
        visits$overall_mean_sp_o2 < 90 ~ "low",
        visits$overall_mean_sp_o2 >= 90 ~ "normal"

    )

    # #Categorise temperature.
    # visits$ord_temp <- case_when(
    #     visits$overall_mean_temp_unknown == 1 ~ "unknown",
    #     visits$overall_mean_temp > 0 ~ "fever", #Preprocessing means this is a temperature above 100.4.
    #     visits$overall_mean_temp == 0 ~ "normal"

    # )

    #Categorise blood pressure ('sbp' has been transformed to 'distance from PALS criterion')
    visits$ord_bp <- case_when(
        visits$overall_mean_sbp_unknown == 1 ~ "unknown",
        visits$overall_mean_sbp > 0 ~ "hypotensive", #Preprocessing means this is below the PALS criterion.
        visits$overall_mean_sbp == 0 ~ "normal"
    )
    

    #Ignore pain.

    #Binarise these variables and add them to the list of filtration characteristics.
    characteristics_to_binarise <- c("ord_num_previous_admissions", "ord_num_previous_visits_without_admission",
        "ord_weight", "ord_overall_mean_hr", "ord_overall_mean_rr",  "ord_spo2", "ord_bp")

    for (i in 1:length(characteristics_to_binarise)) {
        var <- characteristics_to_binarise[i]
        vals <- unique(visits[[var]])
        for (val in vals) {
            visits[[paste0(var, "_", val)]] <- ifelse(visits[[var]]==val, 1, 0)
            filtration_characteristics <- c(filtration_characteristics, paste0(var, "_", val))
        }
    }



    #Define exposures and references. (Here just focus on Black and Hispanic.)
    exposures <- c("sex_f", "race_non_hispanic_black", "race_hispanic")
    refs <- c("sex_m", "race_non_hispanic_white", "race_non_hispanic_white")




    #Now calculate odds ratios of admission, fully adjusted, for each characteristic.
    all_results <- data.frame()
    for (idx in 1:length(exposures)) {
        exposure <- exposures[idx]
        ref <- refs[idx]

        for (characteristic in rev(filtration_characteristics)) {

            print(characteristic)

            
            #For specific tests, look only before the EPIC transition.
            if (characteristic %in% c("received_any_tests")) {
                df <- visits %>% filter(.data[[characteristic]]==1, pre_epic==1)
            } else {
                df <- visits %>% filter(.data[[characteristic]]==1)
            }

            
            ors <- get_odds_ratios(df, propensity_scores, c("is_admitted"),
                exposure, ref, full_adjustment_vars, "full", threshold=threshold, savetag=characteristic)


            if (nrow(ors) > 0) {
                ors$subset <- characteristic
                all_results <- rbind(all_results, ors)
            }
        }
    }

    return(all_results)

}




visits <- read.csv("preprocessed-visits-with-linked-events.csv") %>% select(-c(X.1, V1))



visits <- convert_baseline_factors(visits)
write.csv(visits, "intermediate-files/converted-visits.csv")


visits <- as.data.frame(fread("intermediate-files/converted-visits.csv"))



#Determine variables we should exclude before triage, and for partial adjustment.
all_vars <- colnames(visits)


#Variables we can use when calculating propensity scores before triage (will exclude exposures as needed):
vars_known_at_triage <- c("sex_f", all_vars[startsWith(all_vars, "complaint_contains_")],
     all_vars[startsWith(all_vars, "triage_") & !startsWith(all_vars, "triage_acuity_")], 
     "num_previous_admissions", "num_previous_visits_without_admission", 
     all_vars[startsWith(all_vars, "pre_diagnosis")], all_vars[startsWith(all_vars, "language")], 
     all_vars[startsWith(all_vars, "race_")], all_vars[startsWith(all_vars, "ed_arrival_mode_")],
     all_vars[startsWith(all_vars, "age_group_")], all_vars[startsWith(all_vars, "insurance_")],
     all_vars[startsWith(all_vars, "year_of_arrival")], all_vars[startsWith(all_vars, "season_")],
     all_vars[startsWith(all_vars, "time_of_day_")], "is_weekend", 
     all_vars[startsWith(all_vars, "miles_travelled")],
     all_vars[startsWith(all_vars, "state_of_origin")], all_vars[startsWith(all_vars, "sdi_score")],
     "is_trans_or_nb", all_vars[startsWith(all_vars, "weight")], "pseudo_nedocs")

#Variables we can use when performing 'partial adjustment' across the whole stay
#Defined as 'patient-level' factors
partial_adjustment_vars <- c("sex_f", 
     "num_previous_admissions", "num_previous_visits_without_admission", 
     all_vars[startsWith(all_vars, "pre_diagnosis")], all_vars[startsWith(all_vars, "language")], 
     all_vars[startsWith(all_vars, "race_")], all_vars[startsWith(all_vars, "weight")],
     all_vars[startsWith(all_vars, "age_group_")], all_vars[startsWith(all_vars, "insurance_")], 
     all_vars[startsWith(all_vars, "miles_travelled")],
     all_vars[startsWith(all_vars, "state_of_origin_")], all_vars[startsWith(all_vars, "sdi_score")],
     "is_trans_or_nb")


#Variables we can use when performing full adjustment across the whole stay (exclude triage vitals but include overall)
full_adjustment_vars <-  c("sex_f", all_vars[startsWith(all_vars, "complaint_contains_")],
     all_vars[startsWith(all_vars, "triage_") & !startsWith(all_vars, "triage_acuity_")], 
     "num_previous_admissions", "num_previous_visits_without_admission", all_vars[startsWith(all_vars, "triage_acuity_")],
     all_vars[startsWith(all_vars, "pre_diagnosis")], all_vars[startsWith(all_vars, "current_diagnosis")], all_vars[startsWith(all_vars, "language")], 
     all_vars[startsWith(all_vars, "race_")], all_vars[startsWith(all_vars, "ed_arrival_mode_")],
     all_vars[startsWith(all_vars, "age_group_")], all_vars[startsWith(all_vars, "insurance_")],
     all_vars[startsWith(all_vars, "year_of_arrival")], all_vars[startsWith(all_vars, "season_")],
     all_vars[startsWith(all_vars, "time_of_day_")], "is_weekend", all_vars[startsWith(all_vars, "miles_travelled")],
     all_vars[startsWith(all_vars, "state_of_origin")], all_vars[startsWith(all_vars, "sdi_score")],
     "is_trans_or_nb", all_vars[startsWith(all_vars, "weight")], "pseudo_nedocs", all_vars[startsWith(all_vars, "overall_")])


#Define races.
races <- all_vars[startsWith(all_vars, "race_")]

vars_used <- unique(c(vars_known_at_triage, partial_adjustment_vars, full_adjustment_vars))

#Convert all of these to numeric form
for (var in vars_used) {
    visits[[var]] <- as.numeric(visits[[var]])
}

propensity_scores <- assemble_propensity_scores(visits, races)
write.csv(propensity_scores, "propensity-scores.csv")
propensity_scores <- as.data.frame(fread("propensity-scores.csv"))

admission_ors <- get_admission_odds_ratios(visits, propensity_scores, races, threshold=0.1, savetag="admission")
write.csv(admission_ors, "outputs/caliper_0.1/admission-ors.csv")

admission_ors <- get_admission_odds_ratios(visits, propensity_scores, races, threshold=0.2, savetag="admission")
write.csv(admission_ors, "outputs/caliper_0.2/admission-ors.csv")


#Calculate time quartiles (for ED LOS,  time to admission for admitted patients).
#Rooming time is almost always recorded as the same as arrival time and so can't be usefully compared (the 75th percentile of durations is 0 minutes!)

ed_los_quartiles <- quantile(visits$ed_los, probs = c(0, 0.25, 0.5, 0.75, 1.0), na.rm = TRUE)
admit_time_quartiles <- quantile(as.numeric(visits$time_to_admit_request[which(visits$is_admitted==1)]), probs = c(0, 0.25, 0.5, 0.75, 1.0), na.rm = TRUE)
request_to_admission_time_quartiles <- quantile(as.numeric(visits$time_from_request_to_admission[which(visits$is_admitted==1)]), probs =  c(0, 0.25, 0.5, 0.75, 1.0), na.rm = TRUE)


quartiles <- list(
    "ed_los_quartiles" = ed_los_quartiles,
    "admit_time_quartiles" = admit_time_quartiles,
    "request_to_admission_time_quartiles" = request_to_admission_time_quartiles
)

quartile_names <- list(
    "ed_los_quartiles" = "ed_los",
    "admit_time_quartiles" = "time_to_admit_request",
    "request_to_admission_time_quartiles" = "time_from_request_to_admission"
)


visits <- create_outcomes_database(visits, quartiles, quartile_names)
write.csv(visits, "intermediate-files/visits-with-linked-non-intervention-outcomes.csv")
visits <- as.data.frame(fread("intermediate-files/visits-with-linked-non-intervention-outcomes.csv")) %>% select(-c(V1, V1, x))



upstream_ors <- get_upstream_odds_ratios(visits, propensity_scores, threshold=0.1)
write.csv(upstream_ors, "outputs/caliper_0.1/upstream-ors.csv")


downstream_ors <- get_downstream_odds_ratios(visits, propensity_scores, threshold=0.1)
write.csv(downstream_ors, "outputs/caliper_0.1/downstream-ors.csv")



#Additional investigation: 'who suffers for speaking another language?'
#Match patients on language, then stratify by race and match e.g. black Spanish speakers to black English speakers.


#Define languages
languages <- all_vars[startsWith(all_vars, "language_")]

#Calculate propensity scores (for each language to English.)
language_propensity_scores <- get_propensity_scores(visits, languages, full_adjustment_vars, "full")
write.csv(language_propensity_scores, "intermediate-files/language_propensity_scores.csv")
language_propensity_scores <- as.data.frame(fread("intermediate-files/language_propensity_scores.csv"))


#Investigate the impact of languages, stratified by race.


#There is another set of ORs where the cohorts are filtered (i.e. all those who received a certain intervention/triage score/ED-LOS, all those with a certain form of insurance,
#all those with a certain complaint (of top 10), all those with certain overall vital categories (take worst category as maximum), age, cci-before-visit, num-previous visits, ambulance), and the outcome is admission.
get_language_odds_ratios <- function(visits_to_examine, races, languages, propensity_scores, threshold=0.1) {


    #Copy dataframe
    visits <- as.data.frame(visits_to_examine)

    #Once again, define whites and males:
    visits$race_non_hispanic_white <- ifelse(rowSums(visits[races])==0, 1, 0)
    visits$sex_m <- ifelse(visits$sex_f==0, 1, 0)



    #Now calculate odds ratios of admission, fully adjusted, for each characteristic.
    all_results <- data.frame()
    for (group in c("sex_m", "sex_f", "race_non_hispanic_white", races)) {

        #Stratify by race/sex.
        df <- visits %>% filter(.data[[group]]==1)


        for (idx in 1:length(languages)) {
            exposure <- languages[idx]
            ref <- "language_english"

            print(group)
            print(exposure)


            ors <- get_odds_ratios(df, propensity_scores, c("is_admitted"),
                exposure, ref, full_adjustment_vars, "full", threshold=threshold, savetag=paste0(group, "_", exposure))


            if (nrow(ors) > 0) {
                ors$subset <- group
                all_results <- rbind(all_results, ors)
            }
        }
    }

    return(all_results)

}



language_ors <- get_language_odds_ratios(visits, races, languages, language_propensity_scores)
write.csv(language_ors, "outputs/caliper_0.1/language-ors.csv")

#Define function to get the continuous difference between two continuous variables


#Create an over-arching function which takes a (filtered) dataset, an exposure, an outcome, 
#and a set of variables to adjust for (i.e. a propensity score keyword), and then calculates odds ratios.
get_continuous_differences <- function(filtered_visits, propensity_scores, 
    outcome_vars, exposure, ref, covariates, ps_type, threshold=0.1, savetag="") {

 

    #Define framework.
    odds_ratios <- data.frame(exposure=character(0), 
        outcome_var=character(0), 
        estimate=numeric(0),
        lower_conf=numeric(0),
        upper_conf=numeric(0), 
        standard_error=numeric(0),
        pval=numeric(0), 
        num_treated=numeric(0),
        num_control=numeric(0),
        all_treated=numeric(0),
        treated_mean_outcome=numeric(0),
        control_mean_outcome=numeric(0),
        estimate_in_mins=numeric(0),
        lower_conf_in_mins=numeric(0),
        upper_conf_in_mins=numeric(0))

    index <- 1
    
    df <- data.frame(filtered_visits) 


    #Exclude exposures from covariates.
    #Also exclude all sex-based covariates (e.g. those related to pregnancy) 
    #when the exposure is sex.
    if(length(covariates) > 0  & startsWith(exposure, "race")) {
        covariates <- covariates[!startsWith(covariates, "race")]
    }
    if(length(covariates) > 0 & startsWith(exposure, "sex")) {
        covariates <- covariates[!startsWith(covariates, "sex") &
            !endsWith(covariates, "gynecologic") & 
            !endsWith(covariates, "male_genital") &
            !endsWith(covariates, "menstrual_disorders")]
    }

    #If examining race-specific effect of languages, exclude all relevant variables:
    if (startsWith(exposure, "language")) {
        covariates <- covariates[!(startsWith(covariates, "language"))]
    }

    #Choose the relevant propensity score and thus filter out all patients not in exposure and reference (as they will have NAs)
    ps <- propensity_scores %>% filter(treatment==exposure, type==ps_type)
    df <- df %>% inner_join(ps, by="csn") %>% filter(!is.na(propensity_score))

    for (col in colnames(df)) {
        if (any(is.na(is.numeric(df[[col]])) | is.null(is.numeric(df[[col]])))) {
            print(col)
        }
    }


    #Determine caliper threshold.
    if(ps_type=="unadjusted") {
        threshold_value <- NULL
        #Define formula for balancing
        ps_formula <- as.formula(paste(exposure, "~ propensity_score"))
    } else {
        threshold_value <- threshold*sd(df$propensity_score) #We're not using a caliper threshold if not adjusting.
        #Define formula for balancing
        ps_formula <- as.formula(paste(exposure, "~", paste(covariates, collapse = " + ")))
    }


    #Perform matching.
    tryCatch({

        match_object <- matchit(ps_formula, data=df, method = 'nearest', caliper=threshold_value, distance = df$propensity_score)

 

        if (!(ps_type=="unadjusted")) {
            #Plot only the top 25 variables in a Love plot.
            bal_stats <- bal.tab(match_object)
            ranked_variables <- rownames(bal_stats$Balance[order(abs(bal_stats$Balance$Diff.Adj), decreasing = TRUE), ])
            num_to_take <- min(25, length(ranked_variables))
            ordered_vars <- ranked_variables[1:num_to_take]
            ordered_vars <- ordered_vars[!(ordered_vars=="distance")]
            if ("cci_of_visit_2" %in% ordered_vars) {
                ordered_vars[which(ordered_vars=="cci_of_visit_2")] <- "cci_of_visit" #Sometimes baltab will change the name of this variable, so we change it back.
            }
            p <- love.plot(as.formula(paste(exposure, "~", paste(ordered_vars, collapse = " + "))), data = df, weights = match_object$weights, stats = c("mean.diffs"), 
                    thresholds = c(m = .1, v = 2), abs = TRUE, 
                    binary = "std",
                    continuous = "std", 
                    var.order = "adjusted",
                    stars = "std",
                    s.d.denom = "pooled",
                    ) + theme(axis.text.y = element_text(size = 8))
            ggsave(paste0("outputs/caliper_", threshold, "/", "balance_plots/", savetag, "_", exposure, "_", ps_type, "_balance_plot.pdf"), plot=p)
        }

        #Make two data tables for treated and untreated data (exposure and reference), rejoin on pairing index.
        matched <- as.data.table(match_data(match_object))
        treated_indices <- which(matched[[exposure]]==1)
        control_indices <- which(matched[[exposure]]==0)

        #Calculate the number of patients we could have had before matching.
        all_treated <- sum(df[[exposure]]==1)


        #Run a paired Wilcox test on each outcome.
        for (outcome_var in outcome_vars) {
            
            assert(!any(is.na(df[[outcome_var]])))

            treated_df <- data.frame(pair_index=matched[treated_indices, subclass], treatment_outcome=matched[[outcome_var]][treated_indices])
            control_df <- data.frame(pair_index=matched[control_indices, subclass], control_outcome=matched[[outcome_var]][control_indices])
            overall_table <- inner_join(treated_df, control_df, by="pair_index") 

            treated_logs <- log(pmax(overall_table$treatment_outcome, 1)) #Round up all waits at 0 mins (<0.1%) to 1
            control_logs <- log(pmax(overall_table$control_outcome, 1))


            #Conduct paired t-test on log transformed outcomes
            test <- t.test(treated_logs, control_logs, paired=TRUE, conf.int=TRUE)


            #This gives an estimate of the ratios.

            #Calculate mean outcomes
            treatment_outcome_rate <- mean(df[[outcome_var]][which(df[[exposure]]==1)])
            control_outcome_rate <- mean(df[[outcome_var]][which(df[[exposure]]==0)])

            estimate <- as.numeric(exp(test$estimate))
            lower_conf <- exp(as.numeric(test$conf.int[1]))
            upper_conf <- exp(as.numeric(test$conf.int[2]))

            standard_error <- (upper_conf - lower_conf)/(3.92)


            #Save in table.
            test_result <- c(exposure, outcome_var, estimate, lower_conf, 
                    upper_conf, standard_error, as.numeric(test$p.val),
                    length(treated_indices), length(control_indices), all_treated,
                    treatment_outcome_rate, control_outcome_rate, control_outcome_rate*(estimate-1), #estimate is a ratio of treated outcome to control outcome, so the difference is control*(ratio-1)
                    control_outcome_rate*(lower_conf-1), control_outcome_rate*(upper_conf-1))

            odds_ratios[index,] <- test_result
            index <- index + 1

        } }, error = function(e) { print(paste("Matching failed: ", e$message))})

    
        
    return(odds_ratios)

}


#Check whether distributions are Gaussian after log transformation (they are)

adm_visits <- visits %>% filter(is_admitted==1)
p <- ggplot(visits, aes(x=log(ed_los))) + geom_histogram() + labs(x="Log. visit length", y="Number of visits")
ggsave("outputs/caliper_0.1/ed_los_log_dist.pdf", plot=p)

p <- ggplot(adm_visits, aes(x=log(time_to_admit_request))) + geom_histogram() + labs(x="Log. admit request time", y="Number of visits")
ggsave("outputs/caliper_0.1/admit_request_time_log_dist.pdf", plot=p)

p <- ggplot(adm_visits, aes(x=log(time_from_request_to_admission))) + geom_histogram() + labs(x="Log. bed wait time", y="Number of visits")
ggsave("outputs/caliper_0.1/bed_wait_time_log_dist.pdf", plot=p)

#So, for all relevant groups, get estimate of differences in wait times:
get_differences_in_times <- function(visits, threshold=0.1) {




    #Define exposures and references. (Here just focus on Black and Hispanic.)
    exposures <- c("sex_f", "race_non_hispanic_black", "race_hispanic")
    refs <- c("sex_m", "race_non_hispanic_white", "race_non_hispanic_white")


    #Now calculate odds ratios of admission, fully adjusted, for each characteristic.
    all_results <- data.frame()
    for (idx in 1:length(exposures)) {
        exposure <- exposures[idx]
        ref <- refs[idx]

        #Get differences for admitted patients.

        df <- visits %>% filter(is_admitted==1)

        ors <- get_continuous_differences(df, propensity_scores, c("ed_los", "time_to_admit_request", "time_from_request_to_admission"),
            exposure, ref, full_adjustment_vars, "full", threshold=threshold, savetag="decision-admitted")
        
        if (nrow(ors) > 0) {
            ors$subset <- "admitted"
            all_results <- rbind(all_results, ors)
        }

        #Get ORs for discharged

        df <- visits %>% filter(is_admitted==0)

        ors <- get_continuous_differences(df, propensity_scores, c("ed_los"),
            exposure, ref, full_adjustment_vars, "full", threshold=threshold, savetag="decision-admitted")
        
        if (nrow(ors) > 0) {
            ors$subset <- "discharged"
            all_results <- rbind(all_results, ors)
        }
 
    }

    return(all_results)

}

cts_timing_diffs <- get_differences_in_times(visits)
write.csv(cts_timing_diffs, "outputs/caliper_0.1/cts_timing_diffs.csv")
