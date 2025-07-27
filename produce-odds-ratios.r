
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


setwd("your/working/directory")
save_filepath <- "your/save/filepath"

#Define racial categories.
races <- c("race_hispanic", "race_non_hispanic_black", "race_asian", "race_unknown", "race_non_hispanic_multiracial", "race_other")


#Assemble a list of medications of interest.
corticosteroids <- c("DEXAMETHOSONE", "BUDESONIDE", "HYDROCORTISONE", "METHYLPREDNISOLONE", "PREDNISOLONE", "PREDNISONE", "TRIAMCINOLONE") 
antihistamines <- c("CETIRIZINE")
adrenaline <- c("EPINEPHRINE") #will include racepinephrine
antiepileptic <- c("MIDAZOLAM", "PHENOBARBITAL", "BRIVARACETAM", "CLOBAZAM", "CLONAZEPAM", "DIAZEPAM", "GABAPENTIN", "LACOSAMIDE",
       "PHENYTOIN", "PREGABALIN", "PRIMIDONE", "RUFINAMIDE", "TOPIRAMATE", "VALPROIC ACID", "VIGABATRIN", "ZONISAMIDE") #taken from https://pmc.ncbi.nlm.nih.gov/articles/PMC6130739/
pain_meds <- c("TRAMADOL", "OXYCODONE", "MORPHINE", "FENTANYL", "KETAMINE") #all but the last are opioids
antiemetic <- c("ONDANSETRON")
oral_ibuprofen <- c("ORAL IBUPROFEN")
topical_let <- c("TOPICAL LET")
 

#Also note whether anything is administered intravenously.
medical_dict <- list(
    "corticosteroids"=corticosteroids,
    "antihistamines"=antihistamines,
    "adrenaline"=adrenaline, 
    "antiepileptic"=antiepileptic,
    "pain_meds"=pain_meds, 
    "antiemetic"=antiemetic, 
    "IV"=c("IV"),
    "oral_ibuprofen" =oral_ibuprofen,
    "topical_let"=topical_let)



#Analyse the intersection between race and non-English primary language
plot_race_to_language_intersect <- function(baseline_factors) {

    race_to_language <- as.data.frame(table(baseline_factors$race, baseline_factors$primary_language))  %>% 
        rename(Race=Var1, PrimaryLanguage=Var2) %>% group_by(Race) %>% mutate(TotalPop=sum(Freq)) %>% ungroup() %>%
        mutate(Proportion=Freq/TotalPop)
    write.csv(race_to_language, paste0(save_filepath, "race-language-intersect.csv"))


    race_to_language$PrimaryLanguage <- factor(race_to_language$PrimaryLanguage, levels=c("English", "Spanish", "Other"))

    p <- ggplot(race_to_language, aes(x=Race, y=Proportion, fill=PrimaryLanguage)) + 
        geom_col() + theme_minimal(base_size=18) + theme(axis.text.x = element_text(angle = 45, hjust=1)) 
    
    ggsave(paste0(save_filepath, "race-language-intersect.pdf"), plot=p)
}



#Convert categorical to binary variables for use when fitting propensity scores with glmnet
convert_categorical_to_binary <- function(relevant_baseline_factors) {

    converted_factors <- data.frame(relevant_baseline_factors) 

    #Check for infinite numerical values.
    for (col in colnames(converted_factors)) {
        if (any(is.infinite(as.numeric(converted_factors[[col]])))) {
            print("infinities found") #we find no infinities
        }
    }

    #Convert categorical variables.
    factors_to_convert <- c("sex", "race", "ed_arrival_mode", "state", "primary_language", "triage_acuity", "age_group",
        "hour_of_arrival", "day_of_arrival", "month_of_arrival", "year_of_arrival")
    refs <- c("M", "Non-Hispanic White", "Walk in", "in_state", "English", "Unknown", "older_than_15_years", 1, 1, 1, 2019)
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

    #Examine SDI score and miles travelled.
    #When a value is unknown, replace it with the cohort average, and mark the fact that it's unknown in a separate variables.
    mean_sdi <- mean(as.numeric(converted_factors$cts_sdi_score), na.rm=TRUE)
    idx_to_replace <- which(is.na(as.numeric(converted_factors$cts_sdi_score)))

    converted_factors$cts_sdi_score_unknown <- 0
    converted_factors$cts_sdi_score_unknown[idx_to_replace] <- 1
    converted_factors$cts_sdi_score[idx_to_replace] <- mean_sdi

    #Repeat for miles travelled
    mean_miles_travelled <- mean(as.numeric(converted_factors$cts_miles_travelled), na.rm=TRUE)
    idx_to_replace <- which(is.na(as.numeric(converted_factors$cts_miles_travelled)))

    converted_factors$cts_miles_travelled_unknown <- 0
    converted_factors$cts_miles_travelled_unknown[idx_to_replace] <- 1

    converted_factors$cts_miles_travelled[idx_to_replace] <- mean_miles_travelled

    idx_to_replace <- which(is.na(as.numeric(converted_factors$cts_weight)))

    converted_factors$cts_weight_unknown <- 0
    converted_factors$cts_weight_unknown[idx_to_replace] <- 1
    converted_factors$cts_weight[idx_to_replace] <- 0 #since it should be normalised by weight

    converted_factors <- converted_factors %>% select(-c(all_of(c(factors_to_convert[!(factors_to_convert=="age_group")], "age_in_days")))) %>% clean_names()
    return(converted_factors)
    
}



#Link vitals-- both triage and max/min/mean across stay. Exclude temperature, which is only available for 2023.
link_vitals <- function(converted_baseline_factors, events, age_to_csn,
    vitals_to_take=c("heart_rate", "respiratory_rate", "systolic_blood_pressure", "pain", "oxygen_saturation_spo2"),
    vitals_to_normalise=c("heart_rate", "respiratory_rate")) {
    
    vitals_across_stay <- events %>% 
        filter(event_type=="vitals_taken", event_name %in% vitals_to_take) %>% inner_join(age_to_csn, by="csn")
    vitals_across_stay$event_value <- as.numeric(vitals_across_stay$event_value)

    #Transform BP reading to 'distance from PALS criteria'.
    sbp_readings <- which(vitals_across_stay$event_name=="systolic_blood_pressure")
    sbp <- vitals_across_stay$event_value[sbp_readings]
    days <- vitals_across_stay$age_in_days[sbp_readings]

    vitals_across_stay$event_value[sbp_readings] <- ifelse(days <= 28 & sbp < 60, max(60-sbp, 0) , ifelse(
            days > 28 & 365 & sbp < 70, max(70-sbp, 0), ifelse(
            days > 365 & days < 365.25*10 & sbp < (70 + 2*days/(365.25)), max((70 + 2*days/(365.25))-sbp, 0), ifelse( 
            days >= 365.25*10 & sbp < 90, max(90-sbp, 0), 0))))

    #Isolate triage vitals.
    triage_vitals <- vitals_across_stay %>% group_by(csn) %>% mutate(time_of_first_vital=min(event_time)) %>% 
        ungroup() %>% filter(event_time==time_of_first_vital)

    triage_vitals <- setDT(triage_vitals)
    vitals_across_stay <- setDT(vitals_across_stay)


    #Link both of these to visits.
    triage_vitals <- triage_vitals[, 
        .(mean = mean(event_value)), 
        by = .(csn, event_name)
        ][, 
        dcast(.SD, csn ~ event_name, 
                value.var = c("mean"),
                sep = "_")
        ] 

    triage_vitals <- data.frame(triage_vitals) %>%
            rename_with(~ paste0("triage_vitals_", .), -csn)

    vitals_across_stay <- vitals_across_stay[, 
        .(max = max(event_value), 
            min = min(event_value), 
            mean = mean(event_value)),
        by = .(csn, event_name)
        ][, 
        dcast(.SD, csn ~ event_name, 
                value.var = c("max", "min", "mean"),
                sep = "_")
        ] 

    vitals_across_stay <- data.frame(vitals_across_stay) %>% rename_with(~ paste0("overall_vitals_", .), -csn) 

    converted_baseline_factors <- converted_baseline_factors %>% 
        inner_join(triage_vitals, by="csn") %>% 
        inner_join(vitals_across_stay, by="csn")

    #Normalise HR and RR by age group.
    cols_to_normalise <- colnames(converted_baseline_factors)[grepl(paste(vitals_to_normalise, collapse = "|"), colnames(converted_baseline_factors))]
    for (group in unique(converted_baseline_factors$age_group)) {
        idx <- which(converted_baseline_factors$age_group==group)
        for (col in cols_to_normalise) {
            rel_vitals <- converted_baseline_factors[[col]][idx]
            mean <- mean(rel_vitals, na.rm=TRUE)
            sd <- sd(rel_vitals, na.rm=TRUE)
            converted_baseline_factors[[col]][idx] <- abs(rel_vitals-mean)/sd #get absolute for use in glmnet

        }
    }

    #Mark unknown values
    all_vitals <- colnames(converted_baseline_factors)[grepl("_vitals_", colnames(converted_baseline_factors))]

    for (col in all_vitals) {
        converted_baseline_factors[[paste0(col, "_unknown")]] <- 0
        idx_to_replace <- which(is.na(converted_baseline_factors[[col]]))
        converted_baseline_factors[[paste0(col, "_unknown")]][idx_to_replace] <- 1
        if (col %in% vitals_to_normalise) {
            converted_baseline_factors[[col]][idx_to_replace] <- 0 #because these have been normalised, so the average at each age group should be 0
        } else {
            mean <- mean(converted_baseline_factors[[col]], na.rm=TRUE)
            converted_baseline_factors[[col]][idx_to_replace] <- mean
        }
    }

    #Discard age groups, which are no longer needed.
    converted_baseline_factors <- converted_baseline_factors %>% select(-c(age_group))
    return(converted_baseline_factors)

}

#Calculate propensity scores with and without variables known at triage
#using-all-vitals: FALSE if balancing variables known before triage, TRUE otherwise
#using-all-factors: FALSE if deliberately excluding variables ('partial adjustment'); TRUE otherwise ('full adjustment')
get_propensity_scores <- function(converted_baseline_factors, exposures, using_all_vitals=FALSE, using_all_factors=TRUE) {

    df <- data.frame(converted_baseline_factors)
    propensity_scores <- data.frame()

    #Mark the reference category
    df$ref <- ifelse(rowSums(df[exposures])==0, 1, 0)

    #Drop triage vitals, if assessing across the visit.
    if (using_all_vitals) {
        df <- df %>% select(-all_of(starts_with("triage_vitals_"))) 
    } else {df <- df %>% select(-all_of(c(starts_with("overall_vitals_"), starts_with("triage_acuity_")))) %>% select(-pediatric_comorbidity_score)} #otherwise drop specific triage vitals


    #Drop 'extra' factors, if using partial adjustment.
    if (!(using_all_factors)) {
        cols <- colnames(df)
        to_exclude <- c("pediatric_comorbidity_score", cols[startsWith(cols, "triage_acuity_")], cols[grepl("vital", cols)], "has_medicaid", cols[grepl("complaint", cols)])
        df <- df %>% select(-all_of(cols[cols %in% to_exclude]))
    }

    #Now calculate PS for e.g. being Black relative to being White.
    for (exposure in exposures) {

        #Filter to only exposure and reference (e.g. Black and White).
        relevant_population <- df %>% filter(.data[[exposure]]==1 | ref==1)

        #Shuffle dataframe.
        shuffled_order <- sample(1:nrow(relevant_population), size=nrow(relevant_population), replace=FALSE)
        relevant_population <- relevant_population[shuffled_order,]


        #Convert factors to numerics.
        csns <- relevant_population$csn
        Xs <- relevant_population %>% select(-all_of(c(colnames(relevant_population)[colnames(relevant_population) %in% c(exposures, "ref", "csn", "is_admitted", "x", "X", "X.1", "x_1", "x_x", "x_y")]))) 
        Xs <- as.matrix(as.data.frame(sapply(Xs, as.numeric)))
        print(colnames(Xs))
        ys <- factor(relevant_population[[exposure]], levels=c(0, 1))


        #Fit propensity score with GLMNET.
        ps_model <- cv.glmnet(Xs, ys, family = binomial(link="logit"), nfolds=5, trace.it = 1) #use 5fold CV
        glm_ps <- as.vector(predict(ps_model, type="response", newx=Xs, s=ps_model$lambda.1se)) #use the value of s that gives the lowest error

        #Add to overall dataframe.
        subscores <- data.frame(csn=csns, treatment=exposure, propensity_score=glm_ps, type=ifelse(using_all_vitals, "overall", "triage"),
                completeness=ifelse(using_all_factors, "complete", "partial"))
        propensity_scores <- rbind(propensity_scores, subscores)

    }

    return(propensity_scores)

}


#Assemble all necessary propensity scores.
assemble_propensity_scores <- function(cvtd_factors) {

    #get 12 sets of propensity scores-- predicting sex, race and from triage, across whole visit/just at triage, using partial and full adjustment.
    triage_race_ps_full <- get_propensity_scores(cvtd_factors, races,  using_all_vitals=FALSE)
    overall_race_ps_full <- get_propensity_scores(cvtd_factors, races,  using_all_vitals=TRUE)

    triage_sex_ps_full <- get_propensity_scores(cvtd_factors, c("sex_f"), using_all_vitals=FALSE)
    overall_sex_ps_full <- get_propensity_scores(cvtd_factors, c("sex_f"), using_all_vitals=TRUE)

    triage_race_ps_partial <- get_propensity_scores(cvtd_factors, races,  using_all_vitals=FALSE, using_all_factors = FALSE)
    overall_race_ps_partial <- get_propensity_scores(cvtd_factors, races,  using_all_vitals=TRUE,  using_all_factors = FALSE)

    triage_sex_ps_partial <- get_propensity_scores(cvtd_factors, c("sex_f"), using_all_vitals=FALSE,  using_all_factors = FALSE)
    overall_sex_ps_partial <- get_propensity_scores(cvtd_factors, c("sex_f"), using_all_vitals=TRUE,  using_all_factors = FALSE)


    all_propensity_scores <- rbind(triage_race_ps_full, triage_sex_ps_full,
            triage_race_ps_partial, triage_sex_ps_partial, overall_race_ps_full, overall_sex_ps_full,
            overall_race_ps_partial, overall_sex_ps_partial)

    return(all_propensity_scores)
}


#Calculate odds ratios (also return raw admission rates)
calculate_odds_ratios <- function(outcome_vars, exposure, ref, propensity_scores, baseline_factors, 
        threshold=0.1, ps_method="glm", using_all_vitals=FALSE, using_all_factors=TRUE, unadjusted=FALSE, intervention_df=NULL, target_csns=NULL, characteristic_name=NULL) {

            

    odds_ratios <- data.frame(exposure=character(0), 
        outcome_var=character(0), 
        estimate=numeric(0),
        lower_conf=numeric(0),
        upper_conf=numeric(0), 
        pval=numeric(0), 
        num_treated=numeric(0),
        num_control=numeric(0),
        all_treated=numeric(0),
        treated_outcome_frac=numeric(0),
        control_outcome_frac=numeric(0))

    index <- 1
    
    df <- data.frame(baseline_factors) 


    #Exclude factors not included in propensity scores from balancing.
    illegal_covariates <- c("csn", "is_admitted", "X", "x_1", "x", "X.1", "x_x", "x_y", "ed_los_quartile_1", "ed_los_quartile_2", "ed_los_quartile_3", "ed_los_quartile_4")
    
    #Exclude sex if it's the exposure, race if a racial category is the exposure.
    if (startsWith(exposure, "sex_")) {
        to_exclude <- colnames(df)[startsWith(colnames(df), "sex_")]
        illegal_covariates <- c(illegal_covariates, to_exclude)
    } else {
        if (startsWith(exposure, "race_")) {
            to_exclude <- colnames(df)[startsWith(colnames(df), "race_")]
            illegal_covariates <- c(illegal_covariates, to_exclude)
        }
    }
    #If partially adjusting, exclude unnecessary factors from balancing.
    if (!(using_all_factors)) {
        cols <- colnames(df)
        to_exclude <- c("pediatric_comorbidity_score", cols[startsWith(cols, "triage_acuity_")], cols[grepl("vital", cols)], "has_medicaid", cols[grepl("complaint", cols)])
        illegal_covariates <- c(illegal_covariates, to_exclude)
    }

   #If this is before triage, exclude non-triage vitals.
   if (using_all_vitals) {
        illegal_covariates <- c(illegal_covariates, colnames(df)[startsWith(colnames(df), "triage_vitals_")])
    } else {illegal_covariates  <- c(illegal_covariates, colnames(df)[startsWith(colnames(df), "overall_vitals_") | startsWith(colnames(df), "triage_acuity")], "pediatric_comorbidity_score")} #otherwise drop specific triage vitals


    #Identify the covariates we want to balance.
    covariates <- colnames(df)[!(colnames(df) %in% unique(illegal_covariates))]

    #If we need access to the decisions taken during the stay, link them here.
    if (!is.null(intervention_df)) {
        df <- df %>% inner_join(intervention_df, by="csn")
    }



    #if we want to filter the dataframe for specific CSNs, we can do so here.
    if (!is.null(target_csns)) {
        df <- df %>% filter(csn %in% target_csns)
    }


    #Link relevant propensity score, and don't use any visit which has an NA in that score
    #Now the dataframe is only the exposure and reference category.
    ps <- propensity_scores %>% filter(type==ifelse(using_all_vitals, "overall", "triage"), 
        completeness==ifelse(using_all_factors, "complete", "partial"), treatment==exposure) %>% select(csn, propensity_score)
    df <- df %>% inner_join(ps, by="csn") %>% filter(!is.na(propensity_score)) 



    #If unadjusted==TRUE then this step will have selected a random propensity score based on default covariates
    #and used that to choose patients who are part of either the exposure or reference group.
    #Here we override it and set propensity score to 0.5 for all patients.

    if(unadjusted) {
        df$propensity_score <- 0.5
        covariates <- c("propensity_score")
        threshold_value <- NULL
    } else {
        threshold_value <- threshold*sd(df$propensity_score) #We're not using a caliper threshold if not adjusting.
    }

    #Define formula for balancing
    ps_formula <- as.formula(paste(exposure, "~", paste(covariates, collapse = " + ")))


    #Perform matching.
    tryCatch({
        match_object <- matchit(ps_formula, data=df, method = 'nearest', caliper=threshold_value, distance = df$propensity_score)


        if (!(unadjusted)) {
            #Plot only the top 25 variables in a Love plot.
            bal_stats <- bal.tab(match_object)
            ordered_vars <- rownames(bal_stats$Balance[order(abs(bal_stats$Balance$Diff.Adj), decreasing = TRUE), ][1:25, ])
            ordered_vars <- ordered_vars[!(ordered_vars=="distance")]
            p <- love.plot(as.formula(paste(exposure, "~", paste(ordered_vars, collapse = " + "))), data = df, weights = match_object$weights, stats = c("mean.diffs"), 
                    thresholds = c(m = .1, v = 2), abs = TRUE, 
                    binary = "std",
                    continuous = "std", 
                    var.order = "adjusted",
                    stars = "std"
                    ) + theme(axis.text.y = element_text(size = 8))
            ggsave(paste0(save_filepath, "/", ps_method, "/caliper_", threshold, "/", "balance_plots/", ifelse(is.null(characteristic_name), "", characteristic_name), "_", exposure, "_" , ifelse(using_all_vitals, "overall", "triage"), "_", ifelse(using_all_factors, "complete", "partial"), "_balance_plot.pdf"), plot=p)
        }



        #Make two data tables for treated and untreated data (exposure and reference), rejoin on pairing index.
        matched <- as.data.table(match_data(match_object))
        treated_indices <- which(matched[[exposure]]==1)
        control_indices <- which(matched[[exposure]]==0)

        #Calculate the number of patients we could have had before matching.
        all_treated <- sum(df[[exposure]]==1)


        #Run McNemar's test on each outcome.
        for (outcome_var in outcome_vars) {
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

            #Save in table.
            mcnemar_result <- c(exposure, outcome_var, as.numeric(test$estimate), as.numeric(test$conf.int[1]), 
                    as.numeric(test$conf.int[2]), as.numeric(test$p.val),
                    length(treated_indices), length(control_indices), all_treated,
                    treatment_outcome_rate, control_outcome_rate)

            odds_ratios[index,] <- mcnemar_result
            index <- index + 1

        } }, error = function(e) { print(paste("Matching failed: ", e$message))})
        
    return(odds_ratios)

}

#Get odds ratios describing the likelihood that a patient will be triaged at a certain level, or admitted.
get_triage_and_admission_ors <- function(cvtd_factors, propensity_scores) {

    exposures <- c(races, "sex_f")
    refs <- c(rep("race_non_hispanic_white", length(races)), "sex_m")


    for (threshold in c(0.1, 0.2)) {
        for (method in c("glm")) {

            overall_triage_ors <- data.frame()
            overall_admission_ors <- data.frame()

            for (j in 1:length(exposures)) {
                exposure <- exposures[j]
                ref <- refs[j]

                 #calculate admission_ors with and without all factors
                adm_ors_full <- calculate_odds_ratios(c("is_admitted"), exposure, ref, propensity_scores, cvtd_factors,
                        threshold=threshold, ps_method=method, using_all_vitals=TRUE, using_all_factors = TRUE)
                adm_ors_full$completeness <- "complete"

                adm_ors_partial <- calculate_odds_ratios(c("is_admitted"), exposure, ref, propensity_scores, cvtd_factors,
                        threshold=threshold, ps_method=method, using_all_vitals=TRUE, using_all_factors = FALSE)
                adm_ors_partial$completeness <- "partial"
                
                 adm_ors_unadjusted <- calculate_odds_ratios(c("is_admitted"), exposure, ref, propensity_scores, cvtd_factors,
                        threshold=threshold, ps_method=method, unadjusted=TRUE)
                adm_ors_unadjusted$completeness <- "unadjusted"
                


                triage_ors <- calculate_odds_ratios(paste0(c("triage_acuity_"), c(1, 2, 3, 4, 5)), exposure, ref, propensity_scores, cvtd_factors,
                        threshold=threshold, ps_method=method, using_all_vitals=FALSE, using_all_factors = TRUE)

                
                overall_triage_ors <- rbind(overall_triage_ors, triage_ors)
                overall_admission_ors <- rbind(overall_admission_ors, adm_ors_full, adm_ors_partial, adm_ors_unadjusted)
                
            }

            overall_triage_ors <- overall_triage_ors %>% filter(!is.na(exposure))
            overall_admission_ors <- overall_admission_ors %>% filter(!is.na(exposure))

            write.csv(overall_triage_ors, paste0(save_filepath, "/", method, "/caliper_", threshold, "/triage-ors.csv"))
            write.csv(overall_admission_ors, paste0(save_filepath, "/", method, "/caliper_", threshold, "/adm-ors.csv"))
        }
    }
}


#Calculate the odds ratios for being 'assigned' to each quartile of LOS.
get_length_of_stay_ors <- function(cvtd_factors, events, propensity_scores) {

    exposures <- c(races, "sex_f")
    refs <- c(rep("race_non_hispanic_white", length(races)), "sex_m")
    bf <- data.frame(cvtd_factors)

    #Extract length of stay.
    length_of_stay <- events %>% filter(event_type=="endpoint") %>% select(csn, event_time) %>% rename(ed_los=event_time)
    bf <- bf %>% inner_join(length_of_stay, by="csn")

    #Divide LOS into quartiles.
    quartile_thresholds <- quantile(bf$ed_los, probs=c(0.25, 0.5, 0.75, 1))
    for (j in 1:4) {
        lower_threshold <- ifelse(j==1, 0, quartile_thresholds[j-1])
        upper_threshold <- quartile_thresholds[j]
        quartile_colname <- paste0("ed_los_quartile_", j)
        if (j==4) {
            bf[[quartile_colname]] <- ifelse(bf$ed_los >= lower_threshold, 1, 0)
        } else {
            bf[[quartile_colname]] <- ifelse(bf$ed_los >= lower_threshold
                    & bf$ed_los < upper_threshold, 1, 0)
        }
    }
    bf <- bf %>% select(-ed_los)

    #For both caliper thresholds, calculate ORs for being assigned to each quartile.
    for (threshold in c(0.1, 0.2)) {
        for (method in c("glm")) {

        overall_los_ors <- data.frame()

        for (j in 1:length(exposures)) {
            exposure <- exposures[j]
            ref <- refs[j]

            los_ors <- calculate_odds_ratios(paste0("ed_los_quartile_", c(1, 2, 3, 4)), exposure, ref, propensity_scores, bf,
                    threshold=threshold, ps_method=method, using_all_vitals=TRUE, using_all_factors = TRUE)
            overall_los_ors <- rbind(overall_los_ors, los_ors)

        }

            #Save these.
            overall_los_ors <- overall_los_ors %>% filter(!is.na(exposure))
            write.csv(overall_los_ors, paste0(save_filepath, "/", method, "/caliper_", threshold, "/ed-los-ors.csv"))

        }
    }
}




#Count the number of instances of each medication (medications are only counted specifically if non-IV;
#all IV meds are grouped together as IV.)
assemble_med_data <- function(cvtd_factors, events, medical_dict) {

    overall_med_counts <- data.frame()

    #Count instances of each med.
    for (med_type in names(medical_dict)) {
        med_counts <- events %>% 
            filter(event_name %in% make_clean_names(medical_dict[[med_type]])) %>% #get only corticosteroids etc
            group_by(csn) %>% summarise(count=n()) %>% mutate(type=med_type)
        overall_med_counts <- rbind(overall_med_counts, med_counts)
    }

    #Pivot the table so medications are columns.
    overall_med_counts <- overall_med_counts %>% tidyr::pivot_wider(names_from=type, values_from=count, values_fill = 0)


    #Create extra columns for corticosteroids, counting dosage.
    overall_med_counts$corticosteroids_given_once <- ifelse(overall_med_counts$corticosteroids==1, 1, 0)
    overall_med_counts$corticosteroids_given_two_to_three_times <- ifelse(overall_med_counts$corticosteroids==2 | overall_med_counts$corticosteroids==3, 1, 0)
    overall_med_counts$corticosteroids_given_more_than_three_times <- ifelse(overall_med_counts$corticosteroids>3, 1, 0)

    #Convert NAs to 0s.
    cols_to_convert <- colnames(overall_med_counts)[!(colnames(overall_med_counts)=="csn")]
    for (col in cols_to_convert) {
        overall_med_counts[[col]] <- ifelse(overall_med_counts[[col]]>0, 1, 0)
    }

    #Attach this to the full list of CSNs, in case visits received no meds
    received_meds <- cvtd_factors %>% select(csn) %>% left_join(overall_med_counts, by="csn")
    received_meds[is.na(received_meds)] <- 0 #if NA, that's because it has no entry in med_counts and thus received no meds

    return(received_meds)

}

#Now count radiological tests.
assemble_test_data <- function(cvtd_factors, events) {

    overall_test_counts <- events %>% filter(event_type=="test_ordered") %>%
            group_by(csn, event_name) %>% summarise(ordered=1) %>% 
            tidyr::pivot_wider(names_from=event_name, values_from=ordered, values_fill = 0)
    
    #Attach this to the full list of CSNs, in case patients received no tests
    received_tests <- cvtd_factors %>% select(csn) %>% left_join(overall_test_counts, by="csn")
    received_tests[is.na(received_tests)] <- 0 #if NA, that's because it has no entry in test_counts and thus received no test

    return(received_tests)

}

#For each 'intervention' (meds and tests), calculate the odds ratio 
#that each group will receive it (relative to the reference group).
get_intervention_ors <- function(intervention_outcome_vars, propensity_scores, cvtd_factors, intervention_df, thresholds=c(0.1)) {

    exposures <- c(races, "sex_f")
    refs <- c(rep("race_non_hispanic_white", length(races)), "sex_m")

    #Run this for two different caliper thresholds.
    for (threshold in thresholds) {
        for (method in c("glm")) {

            overall_intervention_ors <- data.frame()

            for (j in 1:length(exposures)) {
                exposure <- exposures[j]
                ref <- refs[j]

                intv_ors <- calculate_odds_ratios(intervention_outcome_vars, exposure, ref, propensity_scores, cvtd_factors,
                        threshold=threshold, ps_method=method, using_all_vitals=TRUE, intervention_df = intervention_df)
                
                overall_intervention_ors <- rbind(overall_intervention_ors, intv_ors)
                
            }

            overall_intervention_ors <- overall_intervention_ors %>% filter(!is.na(exposure))
            write.csv(overall_intervention_ors, paste0(save_filepath, "/", method, "/caliper_", threshold, "/gets-intervention-ors.csv"))
        }
    }
}

#Assemble a dataframe which links 'characteristics' (e.g. SDI scores in a certain quartile, receiving a med) to visits
assemble_filtration_characteristics <- function(baseline_factors, intervention_df, events, num_complaints=20) {


    bf <- data.frame(baseline_factors)

    #Link length of stay.
    length_of_stay <- events %>% filter(event_type=="endpoint") %>% select(csn, event_time) %>% rename(ed_los=event_time)
    filtered_df <- length_of_stay %>% select(csn)
    bf <- bf %>% inner_join(length_of_stay, by="csn")

    #Link age in days, for PALS criterion.
    csn_to_age <- read.csv("baseline-factors-with-top-200-complaint-stems-tagged.csv") %>% select(csn, age_in_days)
    bf <- bf %>% inner_join(csn_to_age, by="csn")

    #Link the 20 most frequent complaints.
    all_complaints <- bf %>% select(starts_with("complaint_contains_"))
    complaints_by_freq <- data.frame(complaint=character(0), count=numeric(0))
    complaint_idx <- 1
    for (complaint in (colnames(all_complaints)[!(colnames(all_complaints)=="X")])) {
        complaints_by_freq[complaint_idx,] <- c(complaint, sum(all_complaints[[complaint]]))
        complaint_idx <- complaint_idx + 1
    }
    complaints_by_freq <- complaints_by_freq %>% arrange(desc(count))
    complaints_to_keep <- complaints_by_freq$complaint[1:num_complaints]


    #Link triage acuity and insurance status.
    non_intv <- bf %>% select(all_of(c("csn", "has_medicaid", "age_in_days", complaints_to_keep)), starts_with("triage_acuity_"))
    non_intv$non_medicaid_insurance <- ifelse(non_intv$has_medicaid==0, 1, 0)

    #Record whether patients have certain variables (weight, min/max HR, min/max RR) a certain number of SDs beyond mean.
    for (col in c("cts_weight", "overall_vitals_max_heart_rate", "overall_vitals_min_heart_rate", "overall_vitals_max_respiratory_rate", "overall_vitals_min_respiratory_rate")) {
        non_intv[[paste0(col, "_less_than_1_sd")]] <- ifelse(bf[[col]] < 1, 1, 0)
        non_intv[[paste0(col, "_1_to_2_sd")]] <- ifelse(bf[[col]] >= 1 & bf[[col]] < 2, 1, 0)
        non_intv[[paste0(col, "_more_than_2_sd")]] <- ifelse(bf[[col]] > 2, 1, 0)

    }

    #Record whether patients are ever hypotensive, by PALS criterion.
    non_intv$hypotensive <- ifelse(non_intv$age_in_days <= 28 & bf$overall_vitals_min_systolic_blood_pressure < 60 |
            non_intv$age_in_days > 28 & 365 & bf$overall_vitals_min_systolic_blood_pressure < 70 |
            non_intv$age_in_days > 365 & non_intv$age_in_days < 365.25*10 & bf$overall_vitals_min_systolic_blood_pressure < (70 + 2*non_intv$age_in_days/(365.25)) |
            non_intv$age_in_days >= 365.25*10 & bf$overall_vitals_min_systolic_blood_pressure < 90, 1, 0)
    non_intv <- non_intv %>% select(-age_in_days)

    #Group patients by expressed pain levels
    non_intv$max_pain_none_or_mild <- ifelse(bf$overall_vitals_max_pain < 4, 1, 0)
    non_intv$max_pain_moderate <- ifelse(bf$overall_vitals_max_pain >= 4 & bf$overall_vitals_max_pain < 7, 1, 0)
    non_intv$max_pain_severe <- ifelse(bf$overall_vitals_max_pain >= 8, 1, 0)

    #Identify patients O2 saturation above and below 90
    non_intv$ox_sat_below_90 <- ifelse(bf$overall_vitals_min_oxygen_saturation_spo2 < 90, 1, 0)
    non_intv$ox_sat_above_90 <- ifelse(bf$overall_vitals_min_oxygen_saturation_spo2 >= 90, 1, 0)

    #Divide contunuous factors into quartiles, excluding unknown values.
    factors_to_quartile <- c("cts_sdi_score", "cts_miles_travelled", "ord_num_patients_at_arrival", "ed_los")

    assert(!any(is.na(bf$ed_los)))
    assert(!any(is.na(bf$ord_num_patients_at_arrival)))

    for (colname in factors_to_quartile) {
        if(paste0(colname, "_unknown") %in% colnames(bf)) {
            quartile_thresholds <- quantile(bf[[colname]][bf[[paste0(colname, "_unknown")]]==0], probs=c(0.25, 0.5, 0.75, 1))
        } else {
            quartile_thresholds <- quantile(bf[[colname]], probs=c(0.25, 0.5, 0.75, 1))
        }
        print(quartile_thresholds)
        for (j in 1:4) {
            lower_threshold <- ifelse(j==1, 0, quartile_thresholds[j-1])
            upper_threshold <- quartile_thresholds[j]
            quartile_colname <- paste0(colname, "_quartile_", j)
            print(lower_threshold)
            print(upper_threshold)
            if(paste0(colname, "_unknown") %in% colnames(bf)) { #if there is an unknown category, exclude those from quartile division
                if (j==4) {
                    non_intv[[quartile_colname]] <- ifelse(bf[[colname]] >= lower_threshold & bf[[paste0(colname, "_unknown")]]==0, 1, 0)
                } else {
                    non_intv[[quartile_colname]] <- ifelse(bf[[colname]] >= lower_threshold
                            & bf[[colname]] < upper_threshold & bf[[paste0(colname, "_unknown")]]==0, 1, 0)
                }
            } else { #e.g. number of patients at arrival- this is always known
                if (j==4) {
                    non_intv[[quartile_colname]] <- ifelse(bf[[colname]] >= lower_threshold, 1, 0)
                } else {
                    non_intv[[quartile_colname]] <- ifelse(bf[[colname]] >= lower_threshold
                            & bf[[colname]] < upper_threshold, 1, 0)
                }
            }
        }
    }



    #Join everything together.
    filtered_df <- inner_join(filtered_df, non_intv, by="csn")
    filtered_df <- inner_join(filtered_df, intervention_df, by="csn")


    return(filtered_df)
}

#Calculate the odds ratios that visits with certain characteristics will then result in an admission
#(e.g. are Black patients with Medicaid less likely to be admitted than NHW patients with Medicaid?)
get_admission_ors_from_specific_populations <- function(filtered_df, propensity_scores, cvtd_factors, 
        threshold=0.1, ps_method="glm") {

            all_ors <- data.frame()

            exposures <- c("sex_f", races)
            refs <- c("sex_m", rep("race_non_hispanic_white", length(races)))

            characteristics <- colnames(filtered_df)[!(colnames(filtered_df) == "csn") & !(startsWith(colnames(filtered_df), "x")) & !(startsWith(colnames(filtered_df), "X"))]
            
            #Examine associations for female, Hispanic and NHB patients (the first 3 exposures.)
            for (exposure_index in 1:3) {
                exposure <- exposures[exposure_index]
                ref <- refs[exposure_index]

                for (characteristic in characteristics) {
                    idx <- which(filtered_df[[characteristic]]==1)
                    csns <- filtered_df$csn[idx]
                    ors <- calculate_odds_ratios(c("is_admitted"), exposure, ref, propensity_scores, cvtd_factors, 
                        threshold=threshold, ps_method=ps_method, using_all_vitals=TRUE, intervention_df=NULL, target_csns=csns, characteristic_name = characteristic) 
                    if (nrow(ors) > 0) {
                        ors$subset <- characteristic
                        all_ors <- rbind(all_ors, ors)
                    }
                }
            }

            write.csv(all_ors, paste0(save_filepath, "/", ps_method, "/caliper_", threshold, "/filtered-subset-admission-ors.csv"))
        }  

#MAIN FUNCTION
get_odds_ratios <- function() {

    #Load data, and limit the analyses to only patients with a recorded admission or discharge.
    baseline_factors <- read.csv("baseline-factors-with-top-200-complaint-stems-tagged.csv", encoding = "UTF-8") 


    to_count <- read.csv("/load/directory/demographics_8May.csv") %>% rename(csn=CSN, mrn=MRN) %>% select(csn, mrn) %>% inner_join(baseline_factors, by="csn")
    print(paste(length(unique(to_count$csn)), "visits from", length(unique(to_count$mrn)), "patients"))

     #Load data, and limit the analyses to only patients with a recorded admission or discharge.
    baseline_factors <- baseline_factors %>% 
        filter(is_admitted==1 | is_discharged==1) %>%
        select(-c(length_of_stay_in_minutes, ed_arrival_time, ed_checkout_time, arrival_timestamp, is_discharged)) 

    to_count <- read.csv("/load/directory/demographics_8May.csv") %>% rename(csn=CSN, mrn=MRN) %>% select(csn, mrn) %>% inner_join(baseline_factors, by="csn")
    print(paste(length(unique(to_count$csn)), "visits from", length(unique(to_count$mrn)), "patients"))


    #Link pediatric comorbidity score
    pci_scores <- read.csv("linked_stays_for_epi.csv") %>% select(csn, ord_PCI) %>% rename(pediatric_comorbidity_score=ord_PCI)
    baseline_factors <- baseline_factors %>% inner_join(pci_scores, by="csn")


    # #Load all events (e.g. administration of meds)
    events <- read.csv("focused-patient-journey-events.csv")

    # #Convert categorical to binary variables.
    cvtd_factors <- convert_categorical_to_binary(baseline_factors)
    write.csv(cvtd_factors, "cvtd-factors-with-tagged-complaint.csv")


    # #Link age in days for PALS criteria.
    age_to_csn <- baseline_factors %>% select(csn, age_in_days) 
    cvtd_factors <- link_vitals(cvtd_factors, events, age_to_csn)



    write.csv(cvtd_factors, paste0(save_filepath, "converted-baseline-factors-with-linked-vitals.csv"))
    cvtd_factors <- read.csv(paste0(save_filepath, "converted-baseline-factors-with-linked-vitals.csv"))



    #Calculate propensity scores
    propensity_scores <- assemble_propensity_scores(cvtd_factors)
    write.csv(propensity_scores, paste0(save_filepath, "propensity-scores.csv"))
    propensity_scores <- read.csv(paste0(save_filepath, "propensity-scores.csv"))

    #Get odds ratios for ESI score and overall admission.
    get_triage_and_admission_ors(cvtd_factors, propensity_scores)

    #Get odds ratios for length of stay.
    get_length_of_stay_ors(cvtd_factors, events, propensity_scores)

    #Load the interventions received by each patient
    med_df <- assemble_med_data(cvtd_factors, events, medical_dict)
    test_df <- assemble_test_data(cvtd_factors, events)

    intervention_df <- inner_join(med_df, test_df, by="csn")
    intervention_outcome_vars <- colnames(intervention_df)[!(colnames(intervention_df) %in% c("X", "X.1", "X.2", "X.x", "X.y", "csn"))]


    #Get odds ratios for other interventions (e.g. meds, tests).
    get_intervention_ors(intervention_outcome_vars, propensity_scores, cvtd_factors, intervention_df)

    filtered_df <- assemble_filtration_characteristics(cvtd_factors, intervention_df, events)
    write.csv(filtered_df, paste0(save_filepath, "filtered_df.csv"))

    #Print how many people have Medicaid at each SDI quartile.
    print(paste("Medicaid uptake in SDI Q1:", sum(filtered_df$has_medicaid[which(filtered_df$cts_sdi_score_quartile_1==1)])/length(which(filtered_df$cts_sdi_score_quartile_1==1))))
    print(paste("Medicaid uptake in SDI Q2:", sum(filtered_df$has_medicaid[which(filtered_df$cts_sdi_score_quartile_2==1)])/length(which(filtered_df$cts_sdi_score_quartile_2==1))))
    print(paste("Medicaid uptake in SDI Q3:", sum(filtered_df$has_medicaid[which(filtered_df$cts_sdi_score_quartile_3==1)])/length(which(filtered_df$cts_sdi_score_quartile_3==1))))
    print(paste("Medicaid uptake in SDI Q4:", sum(filtered_df$has_medicaid[which(filtered_df$cts_sdi_score_quartile_4==1)])/length(which(filtered_df$cts_sdi_score_quartile_4==1))))



    #Get odds ratios for admission post-intervention.
    get_admission_ors_from_specific_populations(filtered_df, propensity_scores, cvtd_factors,
        threshold=0.1, ps_method="glm") 
    # get_admission_ors_from_specific_populations(filtered_df, propensity_scores, cvtd_factors,
    #     threshold=0.2, ps_method="glm") 

}



get_odds_ratios()
