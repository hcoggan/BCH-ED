

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
library(ggforce)
library(forestploter)
library(meta)
library(grid)
library(gridExtra)
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
library(ggforestplot)
library(meta)


setwd("your/working/directory")
 save_filepath <- "your/save/filepath"

#Make forest plots describing the odds ratios in the paper.

method <- "glm"
caliper_threshold <- 0.1

#here are the overall odds for admission from a particular group, at multiple levels of adjustment
adm_ors <- read.csv(paste0(save_filepath, "/", method, "/caliper_", caliper_threshold, "/adm-ors.csv"))

#correct across all levels of adjustment
adm_ors$pval <- p.adjust(adm_ors$pval, method="bonferroni")

#Load this at the second caliper threshold
adm_ors_check <- read.csv(paste0(save_filepath, "/", method, "/caliper_0.2/adm-ors.csv"))

#correct across all levels of adjustment
adm_ors_check$pval <- p.adjust(adm_ors_check$pval, method="bonferroni")


#here are the odds, with corrected p-values, that patients will be 'assigned' to a particular category
subset_assignment_ors <- read.csv(paste0(save_filepath, "/", method, "/caliper_", caliper_threshold, "/gets-intervention-ors.csv")) %>% 
    filter(exposure %in% c("sex_f", "race_hispanic", "race_non_hispanic_black"))
ed_los_ors <- read.csv(paste0(save_filepath, "/", method, "/caliper_", caliper_threshold, "/ed-los-ors.csv")) %>% 
    filter(exposure %in% c("sex_f", "race_hispanic", "race_non_hispanic_black"))
triage_ors <- read.csv(paste0(save_filepath, "/", method, "/caliper_", caliper_threshold, "/triage-ors.csv")) %>%
    filter(exposure %in% c("sex_f", "race_hispanic", "race_non_hispanic_black"))

subset_assignment_ors <- rbind(subset_assignment_ors, triage_ors, ed_los_ors) %>% rename(characteristic=outcome_var)

#Correct across all categories
subset_assignment_ors$pval <- p.adjust(subset_assignment_ors$pval, method="bonferroni")

#Save the corrected results
write.csv(subset_assignment_ors, paste0(save_filepath, "/", method, "/caliper_", caliper_threshold, "/filtered-subset-assignment-corrected-pvals.csv"))

#here are the odds, with corrected p-values, that patients will be admitted *from* that category
subset_admission_ors <- read.csv(paste0(save_filepath, "/", method, "/caliper_", caliper_threshold, "/filtered-subset-admission-ors.csv")) %>%
    filter(exposure %in% c("sex_f", "race_hispanic", "race_non_hispanic_black")) %>% rename(characteristic=subset)

#Correct across all categories
subset_admission_ors$pval <- p.adjust(subset_admission_ors$pval, method="bonferroni")

#Save the corrected results
write.csv(subset_admission_ors, paste0(save_filepath, "/", method, "/caliper_", caliper_threshold, "/filtered-subset-admission-corrected-pvals.csv"))

#PLOT 1: OVERALL ADMISSION RATES AT DIFFERENT ADJUSTMENT LEVELS

adm_ors <- adm_ors %>% rename(Adjustment=completeness)
adm_ors_check <- adm_ors_check %>% rename(Adjustment=completeness)


# Create forest plot using metafor
create_single_outcome_metafor_plot <- function(data) {

    # Convert OR to log OR for metafor
    data$log_or <- log(data$estimate)
    data$log_lower <- log(data$lower_conf)
    data$log_upper <- log(data$upper_conf)

    #label each exposure accordingly
    all_exposures <- c("sex_f", "race_asian", "race_hispanic", "race_non_hispanic_black", "race_non_hispanic_multiracial", "race_other", "race_unknown")
    all_names <- c("Sex: F", "Race: Asian", "Race: Hispanic", "Race: NH Black", "Race: NH Multiracial", "Race: Other", "Race: Unknown")


    data$exposure <- factor(data$exposure, levels=all_exposures, labels=all_names)
    data$Adjustment <- factor(data$Adjustment, levels=c('unadjusted',  'partial', 'complete'))

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
        atransf=exp,
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
                                width = 12, height = 8) {

    
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

save_forest_plots(adm_ors, create_single_outcome_metafor_plot)

save_forest_plots(adm_ors_check, create_single_outcome_metafor_plot, filename="adm_plot_0.2")

#PLOT 2: ADMISSION AND ASSIGNMENT for ESI and LOS

assignment_ors_here <- subset_assignment_ors %>% filter(startsWith(characteristic, "triage") | startsWith(characteristic, "ed_los"))
admission_ors_here <- subset_admission_ors %>% filter(startsWith(characteristic, "triage") | startsWith(characteristic, "ed_los"))

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
    dividers <- c(dividers, 0.5)

    forest(x=data$log_or, 
        ci.lb=data$log_lower,
        ci.ub=data$log_upper,
        slab = labels,
        header = characteristic_name,
        xlab = paste0("Odds Ratio (", xlab, ")"),
        pch = pch,
        refline = 0,
        psize = 1,
        atrans = exp,
        ilab = data.frame(exposure=exposure_names, num_treated=data$num_treated, pval=data$significance),
        ilab.lab = c("Demographic", "N", ""),
        ilab.xpos= xpos,
        ilab.pos = c(4, 2, 2),
        annotate=TRUE,
        cex=cex) 
    abline(h=dividers, lty=2,col='gray70')
} 

# Function to save any of the forest plot methods as EPS
save_multiple_forest_plot <- function(data, plot_function, xlab, characteristic_name, characteristic_levels, characteristic_labels, xpos, xlim,
                                filename = "admission_plot", 
                                width = 12, height = 8, cex = 1.2) {

    

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

save_multiple_forest_plot(assignment_ors_here, create_multiple_outcome_metafor_plot, "ESI score", "Triage acuity", paste0("triage_acuity_", 1:4), paste("ESI", 1:4), filename="triage_assign",
        log(c(0.15, 0.4, 2.7)), log(c(0.15, 5)))
save_multiple_forest_plot(admission_ors_here, create_multiple_outcome_metafor_plot, "Admission to hospital, given ESI score", "Assigned triage acuity", paste0("triage_acuity_", 1:4), paste("ESI", 1:4), filename="triage_admission",
        log(c(0.035, 0.18, 3)), log(c(0.2, 2)))

save_multiple_forest_plot(assignment_ors_here, create_multiple_outcome_metafor_plot, "LOS quartile", "Length of stay", paste0("ed_los_quartile_", 1:4), c("Shortest quartile", "Quartile 2", "Quartile 3", 'Longest quartile'), filename="ed_los_assign",
        log(c(0.55, 0.75, 1.3)), log(c(0.15, 5)))
save_multiple_forest_plot(admission_ors_here, create_multiple_outcome_metafor_plot, "Admission to hospital, given LOS quartile", "Length of stay", paste0("ed_los_quartile_", 1:4), c("Shortest quartile", "Quartile 2", "Quartile 3", 'Longest quartile'), filename="ed_los_admission",
        log(c(0.17, 0.33, 1.2 )), log(c(0.2, 2)))
#PLOT 3: ADMISSION for MAX RESP RATE, HEART RATE, WEIGHT (INSURANCE, SDI, MILES TRAVELLED, CROWDEDNESS?)




admission_ors_here <- subset_admission_ors %>% 
    filter(exposure %in% c("sex_f", "race_non_hispanic_black",  "race_hispanic"), 
        characteristic %in% c("has_medicaid", "non_medicaid_insurance",
            "ox_sat_below_90", "ox_sat_above_90") | startsWith(characteristic, "cts_miles_travelled_") |  startsWith(characteristic, "ord_num_patients_at_arrival_") |
            startsWith(characteristic, "cts_weight_") | startsWith(characteristic, "cts_sdi_") |
            startsWith(characteristic, "overall_vitals_max_heart_") | startsWith(characteristic, "overall_vitals_max_respiratory_"))




print(unique(admission_ors_here$characteristic))

slab <- c("Insurance: Medicaid", "Insurance: Other", paste("Weight:", c("<1sd", "1-2sd", ">2sd")),
    paste("HR:", c("<1sd", "1-2sd", ">2sd")), paste("RR:", c("<1sd", "1-2sd", ">2sd")), "O2: <90", "O2: >90",
    paste0("SDI: Q", 1:4), paste0("Miles travelled: Q", 1:4), paste0("Crowding: Q", 1:4))


save_multiple_forest_plot(admission_ors_here, create_multiple_outcome_metafor_plot, "Admission to hospital", "Characteristic", 
    c("has_medicaid", "non_medicaid_insurance", paste0("cts_weight_", c("less_than_1_sd", "1_to_2_sd", "more_than_2_sd")),
    paste0("overall_vitals_max_heart_rate_", c("less_than_1_sd", "1_to_2_sd", "more_than_2_sd")),
    paste0("overall_vitals_max_respiratory_rate_", c("less_than_1_sd", "1_to_2_sd", "more_than_2_sd")),
    "ox_sat_below_90", "ox_sat_above_90", paste0("cts_sdi_score_quartile_", 1:4), paste0("cts_miles_travelled_quartile_", 1:4),
    paste0("ord_num_patients_at_arrival_quartile_", 1:4)),
    c("Insurance: Medicaid", "Insurance: Other", paste("Weight:", c("<1sd", "1-2sd", ">2sd")),
    paste("HR:", c("<1sd", "1-2sd", ">2sd")), paste("RR:", c("<1sd", "1-2sd", ">2sd")), "O2: <90", "O2: >90",
    paste0("SDI: Q", 1:4), paste0("Miles travelled: Q", 1:4), paste0("Crowding: Q", 1:4)), filename="characteristic_admission",
        log(c(0.12, 0.37, 1.5)), log(c(0.2, 2)), width=8, height=14, cex=1)
