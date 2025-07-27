
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
library(yardstick)
library(ggridges)


setwd("/Volumes/chip-lacava/Groups/BCH-ED/HC-or-work-PSB-run/")
save_filepath <- "/Users/helenacoggan/Documents/OSXLAP12460/BCH-data-2025/HC-or-work-PSB-run/"




#1. RANK SHAP VALUES

#Rename first timestamp column (which registers its actual value, not its SHAP values) to timepoint
#Timestamp column then records SHAP value of timestamp itself

shap_values <- read.csv("shapley-values/time_based_shap_values.csv") %>% rename(timepoint=timestamp, timestamp=timestamp.1)

#rank the importance of all variables, and then of meds

feature_cols <- colnames(shap_values)[!(colnames(shap_values) %in% c("X", "csn", "timepoint"))]

rank_abs_shap_values <- function(shap_values, colnames) {
    df <- shap_values %>% select(all_of(c("csn", "timepoint", colnames))) %>% 
        tidyr::pivot_longer(colnames, names_to="feature", values_to="shap_value") %>%
        group_by(feature) %>% summarise(av_abs_value=mean(abs(shap_value)), label=signif(av_abs_value, 3)) %>%
        arrange(desc(av_abs_value))
    return(df)

}

all_features <- rank_abs_shap_values(shap_values, feature_cols)
write.csv(all_features, "time-based-av-abs-shap-values-all-features.csv")

#LINK PREDICTIONS TO RACE/SEX

# #link sex and race
baseline_factors <- read.csv("baseline-factors-with-top-200-complaint-stems-tagged.csv", encoding = "UTF-8") %>% 
    filter(is_admitted==1 | is_discharged==1) %>%
    select(csn, sex, race, is_admitted)

predictions <- read.csv("shapley-values/time_based_predictions_on_test_data.csv") %>% inner_join(baseline_factors, by="csn")
write.csv(predictions, "shapley-values/time_based_predictions-on-test-data-with-linked-demographics.csv")

predictions <- read.csv("shapley-values/time_based_predictions-on-test-data-with-linked-demographics.csv")

#3. BOOTSTRAP OVERALL AUROCS-- at specific timepoints

bootstrapped_aurocs <- data.frame(timepoint=character(0), mean=numeric(0), lower_conf=numeric(0), upper_conf=numeric(0))
bootstrapped_auprcs <- data.frame(timepoint=character(0), mean=numeric(0), lower_conf=numeric(0), upper_conf=numeric(0))

index <- 1
for (timepoint in c(30, 120, "last")) {
    if (timepoint=="last") { #choose the last timepoint of each visit
        preds <- predictions %>% group_by(csn) %>% mutate(last_timepoint=max(timestamp)) %>% 
            ungroup() %>% filter(timestamp==last_timepoint)
    }
    else {
        preds <- predictions %>% filter(timestamp==timepoint)
    }
    aurocs <- c()
    auprcs <- c()
    
    for (sample in 1:1000) {
        idx <- sample(1:nrow(predictions), size=nrow(predictions), replace=TRUE)
        auprc <- yardstick::pr_auc_vec(truth= factor(preds$is_admitted[idx], levels=c(0, 1)), estimate=preds$ys_pred[idx], case_weights=preds$weight[idx], event_level="second")
        auroc <- yardstick::roc_auc_vec(truth= factor(preds$is_admitted[idx], levels=c(0, 1)), estimate=preds$ys_pred[idx], case_weights=preds$weight[idx], event_level="second")
        auprcs <- c(auprcs, auprc)
        aurocs <- c(aurocs, auroc)
    }

    aurocs <- sort(aurocs)
    bootstrapped_aurocs[index,] <- c(timepoint, mean(aurocs), aurocs[25], aurocs[975])
    auprcs <- sort(auprcs)
    bootstrapped_auprcs[index,] <- c(timepoint, mean(auprcs), auprcs[25], auprcs[975])
    index <- index + 1


}

write.csv(bootstrapped_aurocs, "shapley-values/time_based_boostrapped_aurocs.csv")
write.csv(bootstrapped_auprcs, "shapley-values/time_based_boostrapped_auprcs.csv")

# 4. LINK SHAPLEY VALUES TO RACE/SEX

shap_values <- shap_values %>% inner_join(baseline_factors, by="csn")
write.csv(shap_values, "shapley-values/shap_values_with_race_and_sex.csv")

# 5. GET AUROC/AUPRC VALUES BY TIMEPOINT AND LENGTH OF VISIT

metrics_by_timestamp <- predictions %>% group_by(csn) %>% mutate(max_timestamp = max(timestamp)) %>%
    filter(timestamp > 0, max_timestamp > 0) %>% ungroup() %>%
    group_by(timestamp, max_timestamp) %>% summarise(AUPRC=yardstick::pr_auc_vec(truth= factor(is_admitted, levels=c(0, 1)), estimate=ys_pred, case_weights=weight, event_level="second"),
         AUROC=yardstick::roc_auc_vec(truth= factor(is_admitted, levels=c(0, 1)), estimate=ys_pred, case_weights=weight, event_level="second"))

admission_rate_by_duration <- predictions %>% group_by(csn) %>% mutate(max_timestamp = max(timestamp)) %>%
    filter(timestamp > 0, max_timestamp > 0) %>% ungroup() %>% group_by(max_timestamp) %>%
    summarise(num_visits=n(), AdmissionRate=sum(is_admitted)/num_visits)

write.csv(metrics_by_timestamp, "shapley-values/metrics-by-timestamp.csv")
write.csv(admission_rate_by_duration, "shapley-values/admission-rate-by-timestamp.csv")

#Now plot these.
metrics_by_timestamp <- read.csv("shapley-values/metrics-by-timestamp.csv")
admission_rate_by_duration <- read.csv("shapley-values/admission-rate-by-timestamp.csv")

auroc_plot <- ggplot(metrics_by_timestamp, aes(x=timestamp/60, y=max_timestamp/60, fill=AUROC)) +
    geom_tile() + labs(x="Time since arrival at ED (in hours)", y = "Length of visit (in hours)") + scale_fill_binned(type="viridis") +
    theme_minimal(base_size=19) + theme(legend.position = "bottom", legend.text = element_text(size=13))

auprc_plot <- ggplot(metrics_by_timestamp, aes(x=timestamp/60, y=max_timestamp/60, fill=AUPRC)) +
    geom_tile() + labs(x="Time since arrival at ED (in hours)", y = "Length of visit (in hours)") + scale_fill_binned(type="viridis") +
    theme_minimal(base_size=19) + theme(legend.position = "bottom", legend.text = element_text(size=13))


adm_rate_plot <- ggplot(admission_rate_by_duration, aes(x= AdmissionRate, y = max_timestamp/60)) + ylim(0, max(admission_rate_by_duration$max_timestamp/60)) +
    geom_smooth(color="black") + labs(x="Admission Rate", y = "")  + theme_minimal(base_size=20) + theme(axis.text.y = element_blank())


auroc_plot <- auroc_plot + adm_rate_plot + plot_layout(ncol=2, widths = c(2, 1))
auprc_plot <- auprc_plot + adm_rate_plot + plot_layout(ncol=2, widths = c(2, 1))

ggsave("shapley-values/aurocs-by-timestamp.pdf", plot=auroc_plot)
ggsave("shapley-values/auprcs-by-timestamp.pdf", plot=auprc_plot)

# 6. GET DEMO-SPECIFIC SHAPLEY VALUES FOR NS, WEIGHT (FOR >2SD), HR (FOR >2SD), DISTANCE in UPPER QUARTILE, FEVER,  TIME SINCE ARRIVAL


# link data
test_data <- read.csv("shapley-values/training_fraction_four_years_test_data.csv") 
shap_values <- read.csv("shapley-values/shap_values_with_race_and_sex.csv") 

vars_to_select <- c("sex", "sex", "race", "race", "race")
vals_to_select <- c("M", "F", "Non-Hispanic White", "Non-Hispanic Black", "Hispanic")

#features to analyse
binary_features <- c("triage_acuity_2", "triage_acuity_4", "med_ns_bolus_iv", "complaint_contains_fever", "complaint_contains_abdomin")
cts_features <- c("cts_weight", "max_heart_rate", "cts_miles_travelled", "timestamp", "cts_sdi_score", "ord_num_patients_at_arrival")
cts_feature_thresholds <- c(2, 2, 20.1, 355, 92, 44)

shapley_values_to_plot <- data.frame()

for (i in 1:length(vars_to_select)) {
    var <- vars_to_select[i]
    val <- vals_to_select[i]
    df <- shap_values %>% filter(.data[[var]]==val) #get, e.g. sex = female
    for (bin in binary_features) {
        #get values from test data
        ids <- test_data %>% filter(.data[[bin]]==1) %>% select(csn, timestamp) %>% rename(timepoint=timestamp) #get the identifiers of visits for whom a variable is true
        shapvals <- df %>% inner_join(ids, by=c("csn", "timepoint")) 
        shapley_values_to_plot <- rbind(shapley_values_to_plot,
                data.frame(Demographic=val, 
                           Feature=bin,
                           Value=shapvals[[bin]]))
    }
    #now for continuous features
    for (j in 1:length(cts_features)) {
        #get values from test data
        cts_feat <- cts_features[j]
        cts_thresh <- cts_feature_thresholds[j]
        ids <- test_data %>% filter(.data[[cts_feat]]>=cts_thresh) %>%  select(csn, timestamp) %>% rename(timepoint=timestamp) #get the identifiers of visits for whom a variable is true
        shapvals <- df %>% inner_join(ids, by=c("csn", "timepoint")) 
        shapley_values_to_plot <- rbind(shapley_values_to_plot,
                data.frame(Demographic=val, 
                           Feature=cts_feat,
                           Value=shapvals[[cts_feat]]))
    }

}

write.csv(shapley_values_to_plot, "shapley-values/shap_values_to_plot.csv")

shapley_values_to_plot <- read.csv("shapley-values/shap_values_to_plot.csv") %>% filter(Feature %in% c(c("triage_acuity_2", "cts_miles_travelled",  "timestamp", 
                 "max_heart_rate", "med_ns_bolus_iv")), !(Demographic=="F"| Demographic=="M"))

shapley_values_to_plot$Demographic <- factor(shapley_values_to_plot$Demographic,
        levels=c("M", "F", "Non-Hispanic White", "Non-Hispanic Black", "Hispanic"),
        labels=c("Sex: M", "Sex: F", "Race: NHW", "Race: NHB", "Race: Hispanic"))
shapley_values_to_plot$Feature <- factor(shapley_values_to_plot$Feature,
        levels=c("triage_acuity_4", "triage_acuity_2",  "complaint_contains_fever", "cts_miles_travelled",  "timestamp", 
                 "max_heart_rate",   "med_ns_bolus_iv"),
        labels=c("ESI 4", "ESI 2",  "Complaint relates to fever", "Lives distant from hospital", "Long time since arrival", "Abnormal max. heart rate",  "Received normal saline"))

shapley_labels <- shapley_values_to_plot %>% filter(!(Demographic=="Sex: F"| Demographic=="Sex: M")) %>% group_by(Feature, Demographic) %>% summarise(MeanValue=signif(mean(Value), 2), MaxValue=max(Value))


p <- ggplot(shapley_values_to_plot, aes(x=Feature, y=Value, fill=Demographic)) +
        geom_text(data=shapley_labels, aes(x=Feature, y=1.3, label = MeanValue, fill=Demographic), 
                        position = position_dodge(width=0.8), 
                        hjust = 0.5, vjust = -1, size = 4, color="black") + coord_cartesian(clip="off") +
        geom_boxplot(position = position_dodge(width = 0.8)) + labs(x="Feature", y="SHAP") + theme_minimal(base_size=20)  + theme(axis.text.x = element_text(angle = 45, hjust=1)) 
ggsave("shapley-values/shap-values.pdf", plot=p)