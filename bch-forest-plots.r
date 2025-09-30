

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
library(metafor)

# setwd("/Volumes/chip-lacava/Groups/BCH-ED/reprocessing")
setwd("/Users/helenacoggan/Documents/OSXLAP12460/BCH-data-2025/temp-dict/reprocessing/")


#Make forest plots to describe the BCH results.

#FIRST: a forest plot to describe the odds of admission.
admission_ors <- read.csv("outputs/caliper_0.1/admission-ors.csv")


#Create forest plot using metafor
create_single_outcome_metafor_plot <- function(data) {

#Convert OR to log OR for metafor
data$log_or <- log(data$estimate)
data$log_lower <- log(data$lower_conf)
data$log_upper <- log(data$upper_conf)

#label each exposure accordingly
all_exposures <- c("sex_f", "race_asian", "race_non_hispanic_black", "race_hispanic", "race_hispanic_white", "race_other", "race_unknown")
all_names <- c("Sex: F", "Race: Asian", "Race: NH Black", "Race: Hispanic", "Race: Hispanic White", "Race: Other", "Race: Unknown")


data$exposure <- factor(data$exposure, levels=all_exposures, labels=all_names)
data$Adjustment <- factor(data$adjustment, levels=c('unadjusted',  'partial', 'full'))

data <- data %>% arrange(exposure, Adjustment) %>% mutate(
        significance = case_when(
        pval < 0.001 ~ "***",
        pval < 0.01 ~ "**", 
        pval < 0.05 ~ "*",
        .default = ""
        ))

#Allow a blank label either side of the main label
labels <- c()
pch <- c()
adjustment_level <- c()
dividers <- c()
for (i in 1:length(all_exposures)) {
        labels <- c(labels, "", all_names[i], "")
        pch <- c(pch, 15, 16, 17)
        adjustment_level = c(adjustment_level, "Unadjusted", "Partial", "Full")
        dividers <- c(dividers, 3*i+0.5)

}
dividers <- c(dividers, 0.5)

forest(x=data$log_or, 
        ci.lb=data$log_lower,
        ci.ub=data$log_upper,
        slab = labels,
        header = "Demographic",
        xlab = "Odds Ratio",
        pch = pch,
        refline = 0,
        psize = 1,
        atransf= function(x) exp(x),
        ilab = data.frame(adjustment_level=adjustment_level, num_treated=data$num_treated, pval=data$significance),
        ilab.lab = c("Adjustment Level", "N", ""),
        ilab.xpos = log(c(0.18, 0.36, 1.7)), 
        xlim = log(c(0.1, 3.5)),
        ilab.pos = c(4, 2, 2),
        annotate=TRUE,
        cex=1.2) 
abline(h=dividers, lty=2,col='gray70')
} 



#Save forest plot.
save_forest_plots <- function(data, plot_function, 
                                filename = "admission_plot", 
                                width = 12, height = 8,
                                save_filepath="outputs/caliper_0.1") {


postscript(file = paste0(save_filepath, "/forest-plot/", filename, ".eps"), 
        width = width, 
        height = height, 
        horizontal = FALSE, 
        onefile = FALSE, 
        paper = "special",
        family = "Helvetica")

plot_function(data)

dev.off()

#also save as pdf
pdf(file = paste0(save_filepath, "/forest-plot/", filename, ".pdf"), 
        width = width, 
        height = height, 
        family = "Helvetica")

plot_function(data)

dev.off()

cat("Plot saved as:", filename, "\n")
}

#INITIAL ANALYSIS


#MHC correct across adjustment levels and outcomes.
admission_ors$pval <- p.adjust(admission_ors$pval, method="bonferroni")

#Make a plot for overall admission. (We exclude unknowns as it is a much smaller category and we can't calculate it for all adjustment levels.)
adm_ors <- admission_ors %>% filter(outcome_var=="is_admitted")
save_forest_plots(adm_ors, create_single_outcome_metafor_plot)

#SENSITIVITY-ANALYSIS 2: REPEAT FOR CALIPER THRESHOLD 0.2

admission_ors <- read.csv("outputs/caliper_0.2/admission-ors.csv") #This has already been MHC-corrected.


#MHC correct across adjustment levels and outcomes.
admission_ors$pval <- p.adjust(admission_ors$pval, method="bonferroni")


#Make a plot for overall admission. (We exclude unknowns as it is a much smaller category and we can't calculate it for all adjustment levels.)
adm_ors <- admission_ors %>% filter(outcome_var=="is_admitted")
save_forest_plots(adm_ors, create_single_outcome_metafor_plot, save_filepath="outputs/caliper_0.2")



#Now plot rooming time quartiles, ED-los quartiles, and triage acuity.

#Load upstream ORs.
upstream_ors <- read.csv("outputs/caliper_0.1/upstream-ors.csv") %>% 
        filter(exposure %in% c("sex_f", "race_hispanic", "race_non_hispanic_black")) %>%
        filter(startsWith(outcome_var, "ed_los") | (startsWith(outcome_var, "received_any_") & !(endsWith(outcome_var, "ecgs")))|
        startsWith(outcome_var, "triage_acuity") | startsWith(outcome_var, "rooming_time_")) %>% rename(characteristic=outcome_var)



upstream_ors$pval <- p.adjust(upstream_ors$pval, method="bonferroni")

downstream_ors <- read.csv("outputs/caliper_0.1/downstream-ors.csv") %>% 
        filter(exposure %in% c("sex_f", "race_hispanic", "race_non_hispanic_black")) %>% rename(characteristic=subset) %>%
        filter(!(endsWith(characteristic, "gynecologic") | endsWith(characteristic, "male_genital") | endsWith(characteristic, "menstrual_disorders") |
        endsWith(characteristic, "ecgs") ))

downstream_ors$pval <- p.adjust(downstream_ors$pval, method="bonferroni")

all_characteristics <- unique(downstream_ors$characteristic)



# #Now compare pairs of odds ratios, wherever a `normal` category exists.

# #Create a dictionary of 'reference categories':
# ref_dict <- list(
#     "hr_admission" = "max_hr_normal",
#     "rr_admission" = "max_RR_normal",
#     "bp_admission" = "max_BP_normal",
#     "temp_admission" = "max_Temp_normal",
#     "o2_admission" = "mean_Spo2_normal",
#     "prev_admit_admission" = "ord_num_previous_admissions_0",
#     "prev_visits_admission" = "ord_num_previous_visits_without_admission_0",
#     "cci_before_admission" = "ord_cci_before_visit_0",
#     "age_admission"="age_group_18-29",
#     "insurance_admission" = "insurance_other",
#     "ems_admission" = "ed_arrival_mode_self"
# )

# #Get a list of references
# all_refs <- unlist(ref_dict, use.names=FALSE)
# all_vars <- unique(downstream_ors$characteristic)

# #Initialise the column that will contain the p-values.
# downstream_ors$comparison_to_ref <- 1 #by default

# #Identify variables in the same category
# for (ref in all_refs) {


#     #Remove everything after the last underscore to find the stem
#     stem <- str_split(ref, "_")[[1]]
#     stem <- stem[1:length(stem)-1] 
#     stem <- paste(stem, collapse="_")

#     #Identify all variables starting with the same stems
#     related_vars <- all_vars[startsWith(all_vars, stem) & !(all_vars==ref)]

#     #Now compare each related variable to the reference variable at each exposure
#     for (var in related_vars) {
#         for (exp in unique(downstream_ors$exposure)) {
#             #print(paste(var, ref, exp))
#             var_index <- which(downstream_ors$characteristic==var & 
#                 downstream_ors$exposure==exp)
        
#             ref_index <- which(downstream_ors$characteristic==ref & 
#                 downstream_ors$exposure==exp)

#             #Get log odds and standard errors
#             var_log_or <- log(downstream_ors$estimate[var_index])
#             ref_log_or <- log(downstream_ors$estimate[ref_index])

#             var_se <- downstream_ors$standard_error[var_index]
#             ref_se <- downstream_ors$standard_error[ref_index]

#             #Get the test statistic for the difference
#             statistic <- (var_log_or-ref_log_or)/sqrt(var_se**2 + ref_se**2)

#             #Get a p-value.
#             pval <- 2*pnorm(q=statistic, lower.tail=FALSE)

#             downstream_ors$comparison_to_ref[var_index] <- pval

#         }

#     }


# }

# #Do not bother MHC correcting the comparisons of pairs of ORs.


# Create forest plot using metafor
create_multiple_outcome_metafor_plot <- function(df, xlab, characteristic_name, characteristic_levels, characteristic_labels, xpos=NULL, xlim=NULL, cex=1.2) {

    data <- df %>% filter(characteristic %in% characteristic_levels)

    #label each exposure accordingly
    all_exposures <- c("sex_f",  "race_non_hispanic_black", "race_hispanic")
    all_names <- c("Sex: F", "Race: NH Black",  "Race: Hispanic")

    data$log_or <- log(data$estimate)
    data$log_lower <- log(data$lower_conf)
    data$log_upper <- log(data$upper_conf)



    data$exposure <- factor(data$exposure, levels=all_exposures, labels=all_names)
    data$characteristic <- factor(data$characteristic, levels=characteristic_levels, labels=characteristic_labels)

    data <- data %>% arrange(characteristic, exposure) %>% mutate(
            significance = case_when(
            pval < 0.001 ~ "***",
            pval < 0.01 ~ "**", 
            pval < 0.05 ~ "*",
            .default = ""
            ))

    # if("comparison_to_ref" %in% colnames(data)) {
    #         #Highlight an OR in red if it is significantly different for this category than for the 'normal' category.
    #         data$color <- ifelse(data$comparison_to_ref < 0.05, "#880808", "#000000")
    # } else {
    #     data$color <- "#000000"
    # }



    #Allow a blank label either side of the main label
    labels <- c()
    pch <- c()
    exposure_names <- c()
    dividers <- c()
    num_blanks <- floor(length(characteristic_levels)/2)
    for (i in 1:length(characteristic_levels)) {
            labels <- c(labels, "", characteristic_labels[i], "")
            pch <- c(pch, 15, 16, 17)
            exposure_names <- c(exposure_names, all_names)
            dividers <- c(dividers, 3*i+0.5)

    }



    forest(x=data$log_or, 
            ci.lb=data$log_lower,
            ci.ub=data$log_upper,
            slab = labels,
            header = characteristic_name,
            xlab = paste0("Odds Ratio (", xlab, ")"),
            pch = pch,
            psize = 1,
            ilab = data.frame(exposure=exposure_names, num_treated=data$num_treated, pval=data$significance),
            ilab.lab = c("Demographic", "N", ""),
            ilab.xpos= xpos,
            ilab.pos = c(4, 2, 2),
            #col=data$color,
            annotate=TRUE,
            refline = 0,
            efac=0.5,
            atransf = function(x) exp(x),
            cex=cex) 

    abline(h=dividers, lty=2,col='gray70')
    } 


    # Function to save any of the forest plot methods as EPS
    save_multiple_forest_plot <- function(data, plot_function, xlab, characteristic_name, characteristic_levels, characteristic_labels, xpos, xlim,
                                    filename = "admission_plot", 
                                    width = 12, height = 8, cex = 1.2, save_filepath="outputs/caliper_0.1") {



    # Open EPS device
    postscript(file = paste0(save_filepath, "/forest-plot/", filename, ".eps"), 
            width = width, 
            height = height, 
            horizontal = FALSE, 
            onefile = FALSE, 
            paper = "special",
            family = "Helvetica")

    # Execute the plot function

    plot_function(data, xlab, characteristic_name, characteristic_levels, characteristic_labels, xpos, xli, cex=cex)

    # Close the device
    dev.off()

    #also save as pdf
    # Open EPS device
    pdf(file = paste0(save_filepath, "/forest-plot/", filename, ".pdf"), 
            width = width, 
            height = height, 
            family = "Helvetica")
    # Execute the plot function
    plot_function(data, xlab, characteristic_name, characteristic_levels, characteristic_labels, xpos, xlim, cex=cex)

    dev.off()

    cat("Plot saved as:", filename, "\n")
}

#TRIAGE ACUITY SCORE
save_multiple_forest_plot(upstream_ors, create_multiple_outcome_metafor_plot, "ESI score", "Triage acuity", paste0("triage_acuity_", 1:4), paste("ESI", 1:4), filename="triage_assign",
        log(c(0.2, 0.6, 2.75)), cex=1.4)
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital, given ESI score", "Assigned triage acuity", paste0("triage_acuity_", 1:4), paste("ESI", 1:4), filename="triage_admission",
        log(c(0.08, 0.32, 2.4)), cex=1.4)


#ED LOS
save_multiple_forest_plot(upstream_ors, create_multiple_outcome_metafor_plot, "LOS quartile", "Length of stay", paste0("ed_los_quartiles_", 1:4), c("Shortest quartile", "Quartile 2", "Quartile 3", 'Longest quartile'), filename="ed_los_assign",
        log(c(0.6, 0.8, 1.22)), log(c(0.15, 5)), cex=1.4)
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital, given LOS quartile", "Length of stay", paste0("ed_los_quartiles_", 1:4), c("Shortest quartile", "Quartile 2", "Quartile 3", 'Longest quartile'), filename="ed_los_admission",
        log(c(0.3, 0.5, 1.1)), log(c(0.2, 2)), cex=1.4)

#INTERVENTIONS
save_multiple_forest_plot(upstream_ors, create_multiple_outcome_metafor_plot, "Interventions", "Received...", 
        c("received_any_labs", "received_any_tests", "received_any_medications", "received_any_iv"),
        c( "Lab", "Imaging", "Med.", "IV med."), filename="intervention_assign",
        log(c(0.4, 0.65, 1.3)), log(c(0.15, 5)), cex=1.4)
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital, given interventions", "Received...",
        c( "received_any_labs", "received_any_tests", "received_any_medications", "received_any_iv"),
        c( "Lab", "Imaging", "Med.", "IV med."), filename="intervention_admission",
        log(c(0.4, 0.67, 1.48)), log(c(0.2, 2)), cex=1.4)


#CHARACTERISTICS- go one by one

#Heart rate
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        c("ord_overall_mean_rr_2_or_more_sd", "ord_overall_mean_hr_1_to_2_sd", "ord_overall_mean_hr_less_than_1_sd"),
        c("HR: >2sd", "HR: 1-2sd", "HR: <1sd"), filename="hr_admission",
        log(c(0.4, 0.62, 1.2)), log(c(0.2, 2)), cex=1.4)

#Resp rate
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        c("ord_overall_mean_rr_2_or_more_sd", "ord_overall_mean_rr_1_to_2_sd", "ord_overall_mean_rr_less_than_1_sd"),
        c("RR: >2sd", "RR: 1-2sd", "RR: <1sd"), filename="rr_admission",
        log(c(0.45, 0.67, 1.2)), log(c(0.2, 2)), cex=1.4)
# #Temp
# save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
#         c( "ord_temp_fever", "ord_temp_normal"),
#         c("Fever", "No fever"), filename="temp_admission",
#         log(c(0.35, 0.55, 1.05)), log(c(0.2, 2)), cex=1.4)

#O2
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        c( "ord_spo2_low", "ord_spo2_normal"),
        c( "O2 < 90", "O2 > 90"), filename="o2_admission",
        log(c(0.02, 0.1, 1.65)), log(c(0.2, 2)), cex=1.4)

#bp
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        c( "ord_bp_hypotensive", "ord_bp_normal"),
        c("Hypotensive", "Normal"), filename="bp_admission",
        log(c(0.20, 0.4, 1.2)), log(c(0.2, 2)), cex=1.4)

#age,
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        c( "age_group_under_3_months", "age_group_three_to_6_months", "age_group_six_to_12_months", "age_group_twelve_to_18_months",
         "age_group_eighteen_months_to_3_years", "age_group_three_to_5_years", "age_group_five_to_10_years",
         "age_group_ten_to_15_years", "age_group_fifteen_and_older"),
        c("0-3 mo", "3-6 mo", "6-12 mo", "12-18 mo", "18-36 mo", "3-5 y", "5-10 y", "10-15 y", "15+ y"), filename="age_admission",
        log(c(0.25, 0.5, 1.3)), log(c(0.2, 2)), cex=1.4)

#insurance
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        c( "insurance_public", "insurance_private"),
        c("Insurance: Public", "Insurance: Private"), filename="insurance_admission",
        log(c(0.42, 0.63, 1.1)), log(c(0.2, 2)), cex=1.4)

# ed arrival mode
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        c( "ed_arrival_mode_ems", "ed_arrival_mode_transfer", "ed_arrival_mode_walk-in"),
        c("EMS", "Transfer", "Walk-in"), filename="ems_admission",
        log(c(0.35, 0.6, 1.2)), log(c(0.2, 2)), cex=1.4)


#Weight
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        c("ord_weight_2_or_more_sd", "ord_weight_1_to_2_sd", "ord_weight_less_than_1_sd"),
        c("Weight: >2sd", "Weight: 1-2sd", "Weight: <1sd"), filename="weight_admission",
        log(c(0.2, 0.4, 1.3)), log(c(0.2, 2)), cex=1.4)


#number of admissions
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        c( "ord_num_previous_admissions_0", "ord_num_previous_admissions_1", "ord_num_previous_admissions_2_or_more"),
        c("No prev. admits", "1 prev. admit", "2+ prev. admits"), filename="prev_admits_admission",
        log(c(0.02, 0.3, 7)), log(c(0.2, 2)), cex=1.4)

#number of visits without admission
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        c( "ord_num_previous_visits_without_admission_0", "ord_num_previous_visits_without_admission_1", "ord_num_previous_visits_without_admission_2_or_more"),
        c("No prev. discharges", "1 prev. discharge", "2+ prev. discharges"), filename="prev_visits_admission",
        log(c(0.2, 0.6, 2)), log(c(0.2, 2)), cex=1.4)

#language
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        c( "language_english",  "language_portuguese", "language_spanish", "language_other"),
        c("English", "Portuguese", "Spanish", "Other"),
         filename="language_admission",
        log(c(0.015, 0.1, 2.5)), log(c(0.2, 2)), cex=1.4)

#SDI score
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        c( "sdi_score_1", "sdi_score_2", "sdi_score_3", "sdi_score_4"),
        c("Least deprived quartile", "Quartile 2", "Quartile 3", "Most deprived quartile"),
         filename="sdi_admission",
        log(c(0.41, 0.6, 1.15)), log(c(0.2, 2)), cex=1.4)

#Miles travelled
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        c( "miles_travelled_1", "miles_travelled_2", "miles_travelled_3", "miles_travelled_4"),
        c("Most local quartile", "Quartile 2", "Quartile 3", "Most distant quartile"),
         filename="distance_admission",
        log(c(0.3, 0.65, 1.9)), log(c(0.2, 2)), cex=1.4)

#Crowding
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        c( "pseudo_nedocs_1", "pseudo_nedocs_2", "pseudo_nedocs_3", "pseudo_nedocs_4"),
        c("Least crowded quartile", "Quartile 2", "Quartile 3", "Most crowded quartile"),
         filename="crowding_admission",
        log(c(0.52, 0.7, 1.11)), log(c(0.2, 2)), cex=1.4)

vars <- c("insurance_public", "insurance_private", "ord_weight_2_or_more_sd", "ord_weight_1_to_2_sd", "ord_weight_less_than_1_sd", "ord_overall_mean_hr_2_or_more_sd", "ord_overall_mean_hr_1_to_2_sd", "ord_overall_mean_hr_less_than_1_sd",  "ord_overall_mean_rr_2_or_more_sd", "ord_overall_mean_rr_1_to_2_sd", "ord_overall_mean_rr_less_than_1_sd",
                   "sdi_score_1", "sdi_score_2", "sdi_score_3", "sdi_score_4", 
                         "miles_travelled_1", "miles_travelled_2", "miles_travelled_3", "miles_travelled_4",  "pseudo_nedocs_1", "pseudo_nedocs_2", "pseudo_nedocs_3", "pseudo_nedocs_4")
labels <- c("Insurance: Public", "Insurance: Private", "Weight: >2sd", "Weight: 1-2sd", "Weight: <1sd", "HR: >2sd", "HR: 1-2sd", "HR: <1sd", "RR: >2sd", "RR: 1-2sd", "RR: <1sd", "SDI: Q1", "SDI: Q2", "SDI: Q3", "SDI: Q4", 
                "Miles travelled: Q1", "Miles travelled: Q2", "Miles travelled: Q3", "Miles travelled: Q4", "Crowding: Q1", "Crowding: Q2", "Crowding: Q3", "Crowding: Q4")

#Overall
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        vars,
        labels,
         filename="characteristic_admission",
        log(c(0.165, 0.48, 1.84)), log(c(0.2, 2)), cex=1.4, width=10, height=18)





#List pre-diagnoses in descending order of visits, excluding vaginal complaints and pregnancy, and shingles/other chest/pelvic (for lack of sample size/no ORs)
visits <- read.csv("preprocessed-visits-with-linked-events.csv") 
comorbidities <- colnames(visits)[startsWith(colnames(visits), "pre_diagnosis_")]
comorbidities <-  comorbidities[!(comorbidities %in% c("pre_diagnosis_menstrual_disorders"))]

comorbidities_by_freq <- data.frame(comorb=character(0), num_visits=numeric(0))
index <- 1
for (comorb in comorbidities) {
        comorbidities_by_freq[index,] <- c(comorb, sum(visits[[comorb]]))
        index <- index + 1
}


comorbidities_by_freq <- comorbidities_by_freq %>% arrange(desc(as.numeric(num_visits)))
write.csv(comorbidities_by_freq, "comorbidities-by-freq.csv")

comorbidities_by_freq <- read.csv("comorbidities-by-freq.csv")
comorbidities <- comorbidities_by_freq$comorb

#comorbidities (top 10 in descending order of frequency)
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        c("pre_diagnosis_pain_conditions", "pre_diagnosis_nausea_and_vomiting", "pre_diagnosis_asthma", "pre_diagnosis_congenital_malformations",
        "pre_diagnosis_gastrointestinal", "pre_diagnosis_epilepsy", "pre_diagnosis_developmental_delays", "pre_diagnosis_cardiovascular", 
        "pre_diagnosis_anxiety", "pre_diagnosis_depression"),
        c("Pain", "NVD", "Asthma", "Congenital", "GI", "Epilepsy", "Developmental", "CV", "Anxiety", "Depression"),
        filename="pre_conditions_1_admission",
        log(c(0.12, 0.35, 1.8)), log(c(0.2, 2)), cex=1.4)


#comorbidities
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        c( "pre_diagnosis_any_malignancy", "pre_diagnosis_anemia",  "pre_diagnosis_conduct_disorders", "pre_diagnosis_chromosomal_anomalies",
         "pre_diagnosis_diabetes_mellitus", "pre_diagnosis_sleep_disorders" ,  "pre_diagnosis_psychotic_disorders", "pre_diagnosis_eating_disorders",  
        "pre_diagnosis_weight_loss", "pre_diagnosis_drug_abuse"  ),
        c("Cancer", "Anaemia", "Behavioural", "Chromosomal", "Diabetes mel.", "Sleep", "Psychosis", "ED", "Weight loss", "Drug abuse"),
        filename="pre_conditions_2_admission",
        log(c(0.005, 0.08, 7)), log(c(0.2, 2)), cex=1.4)




#List complaints in descending order of visits, excluding vaginal complaints and pregnancy, and shingles/other chest/pelvic (for lack of sample size/no ORs)
visits <- read.csv("preprocessed-visits.csv") 

complaints <- colnames(visits)[startsWith(colnames(visits), "complaint_")]


complaint_by_freq <- data.frame(complaint=character(0), num_visits=numeric(0))
index <- 1
for (complaint in complaints) {
        complaint_by_freq[index,] <- c(complaint, sum(visits[[complaint]]))
        index <- index + 1
}

write.csv(complaint_by_freq, "complaint-by-freq.csv")

complaint_by_freq <- read.csv("complaint-by-freq.csv")
complaint_by_freq <- complaint_by_freq %>% arrange(desc(as.numeric(num_visits)))
complaints <- complaint_by_freq$complaint


complaint_names <- c( "Fever", "Abdominal pain", "Extremity", "Other respiratory",
        "Vomiting", "Cough", "Rash", "Psychiatric", "Laceration", "Head or neck", "Trauma",
        "Ear", "Seizure", "Asthma/wheezing", "Headache", "Lump/mass", "Eye", "Sore throat", "Chest pain", "Dental")

#all complaints-- order them by frequency
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        complaints[1:10], complaint_names[1:10], filename="complaints_top_10_admission",
        log(c(0.05, 0.3, 4.4)), log(c(0.2, 2)), cex=1.4)

#next complaints-- order them by frequency
save_multiple_forest_plot(downstream_ors, create_multiple_outcome_metafor_plot, "Admission to hospital", "",
        complaints[11:20], complaint_names[11:20], filename="complaints_next_10_admission",
        log(c(0.0008, 0.03, 3)), log(c(0.2, 2)), cex=1.4)







language_ors <-  read.csv("outputs/caliper_0.1/language-ors.csv") %>% 
        filter(subset %in% c("race_non_hispanic_white", "race_non_hispanic_black", "race_hispanic")) %>% rename(characteristic=subset)
        
language_ors$pval <- p.adjust(language_ors$pval, method="bonferroni")

#Make language-specific plot.

print(head(filter(language_ors, pval < 0.05)))

# Open EPS device
pdf(file = paste0("outputs/caliper_0.1/forest-plot/language-ors.pdf"), 
         width = 12, 
        height = 8,
        family = "Helvetica")


exposure_names <- c("NHW/Arabic", "NHW/Spanish", "NHW/Portuguese", "NHW/Other",
        "NHB/Cape Verdean", "NHB/Haitian Creole", "NHB/Spanish", "NHB/Portuguese", "NHB/Other",
        "Hispanic/Cape Verdean", "Hispanic/Haitian Creole", "Hispanic/Spanish", "Hispanic/Portuguese", "Hispanic/Other")

races <- c(rep("race_non_hispanic_white", 4), rep("race_non_hispanic_black", 5), rep("race_hispanic", 5))
languages <- c("language_arabic", "language_spanish", "language_portuguese", "language_other",
        rep(c("language_cape_verdean", "language_haitian_creole", "language_spanish", "language_portuguese", "language_other"), 2))

data <- data.frame(exposure=languages, characteristic=races) %>% inner_join(language_ors, by=c("exposure", "characteristic"))


data$log_or <- log(data$estimate)
data$log_lower <- log(data$lower_conf)
data$log_upper <- log(data$upper_conf)


data <- data %>% mutate(
        significance = case_when(
        pval < 0.001 ~ "***",
        pval < 0.01 ~ "**", 
        pval < 0.05 ~ "*",
        .default = ""
        ))


forest(x=data$log_or, 
        ci.lb=data$log_lower,
        ci.ub=data$log_upper,
        header = "Race/Language",
        xlab = "Odds Ratio (admission, relative to same-race English speakers)",
        slab=exposure_names,
        psize = 1,
        ilab = data.frame(num_treated=data$num_treated, pval=data$significance),
        ilab.lab = c("N", ""),
        ilab.xpos= log(c(0.1, 7)),
        ilab.pos= c(2, 2),
        #col=data$color,
        annotate=TRUE,
        refline = 0,
        atransf = function(x) exp(x),
        cex=1.4)






#also save as pdf


dev.off()

