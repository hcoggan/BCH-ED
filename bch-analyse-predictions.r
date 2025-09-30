
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
library(ggallin)


#Ask whether different groups have different assigned Shapley values for each of the top X factors predicting admission;
#whether they have different F1 scores, and whether certainty is achieved more quickly for any given group.

# setwd("/Volumes/chip-lacava/Groups/BCH-ED/reprocessing/")
setwd("/Users/helenacoggan/Documents/OSXLAP12460/BCH-data-2025/temp-dict/reprocessing")
set.seed(2115)

# #Load shapley values
# shapley_values <- read.csv("prediction/time_based_shap_values.csv") %>%
#     rename(timepoint=snapshot, snapshot=snapshot.1) #'snapshot' refers to its Shapley value

# #Load visits
# visits <- read.csv("preprocessed-visits.csv")


# #Load predictions
# test_data <- read.csv("prediction/test-data.csv")
# #rain_data <- read.csv("prediction/training-data.csv")
# predictions <- read.csv("prediction/time_based_predictions_on_test_data.csv")
# predictions$snapshot <- as.numeric(predictions$snapshot)




# # print(head(predictions))


#In the order of the beeswarm plot:
features_to_plot <- c('triage_acuity_4', 'triage_acuity_2', 'miles_travelled', 'snapshot', 'complaint_contains_fever', 'lab_potassium', 'lab_hematocrit', 'ed_arrival_mode_transfer', 'triage_acuity_3', 'pseudo_nedocs', 'complaint_contains_extremity', 'weight', 'min_sp_o2', 'complaint_contains_psych', 'complaint_contains_abdominal_pain', 'lab_influenza_a_poct', 'medroute_d5w_ns_1_000_ml_iv',
     'complaint_contains_vomiting', 'lab_c_reactive_protein', 'ed_arrival_mode_ems', 'max_hr', 'year_of_arrival_2022', 'lab_hemoglobin', 'max_rr')

# #First get the actual Shapley values:
# shap_by_race <- shapley_values %>% inner_join(select(visits, c(csn, race)), by="csn") %>%
#     filter(race %in% c("Hispanic", "Non-Hispanic White", "Non-Hispanic Black")) %>%
#     tidyr::pivot_longer(!c(csn, timepoint, race), names_to="feature", values_to="shap_value") %>%
#     filter(feature %in% features_to_plot)

# #Now BOOTSTRAP the mean absolute shapley value for each feature
# bootstrapped_shap_values <- data.frame()

# bootstrap_shap_values <- function(i) {
#     idx <- sample(1:nrow(shap_by_race), size=nrow(shap_by_race), replace=TRUE)
#     to_add <- shap_by_race[idx,]
#     to_add$estimate <- i 
#     return(to_add)

# }
# bootstrapped_shap_values <- purrr::map_df(1:1000, bootstrap_shap_values, .progress=TRUE) %>% rbind()


# bootstrapped_shap_values <- bootstrapped_shap_values %>%
#     group_by(feature, race, estimate) %>% mutate(importance=mean(abs(shap_value))) %>%
#     group_by(feature, race) %>% summarise(mean=mean(importance, na.rm=TRUE), lower_conf=quantile(importance, 0.025, na.rm=TRUE), upper_conf=quantile(importance, 0.975, na.rm=TRUE))
# write.csv(bootstrapped_shap_values, "prediction/bootstrapped-shap-value.csv")

bootstrapped_shap_values <- read.csv("prediction/bootstrapped-shap-value.csv")



#Make sure values are in the right order
bootstrapped_shap_values$feature <- factor(bootstrapped_shap_values$feature, levels=rev(features_to_plot))
bootstrapped_shap_values$race <- factor(bootstrapped_shap_values$race, levels=rev(c("Hispanic", "Non-Hispanic Black", "Non-Hispanic White")))


#Make a plot with bars horizontally aligned
p <- ggplot(bootstrapped_shap_values) +
    geom_col(aes(x=mean, y=feature, fill=race), position=position_dodge(width=0.8)) +
    geom_errorbarh(aes(y=feature, group=race, xmin=lower_conf, xmax=upper_conf), color = "black", height = 0.2, alpha = 0.4,
        position=position_dodge(width=0.8)) +
    labs(x="Average feature importance",
        y="", fill="Race") + theme_minimal(base_size=20) + theme(axis.ticks.y = element_blank(),
      axis.text.y = element_blank()) +  scale_fill_brewer(palette = "Dark2") + 
      guides(fill = guide_legend(reverse=TRUE))
ggsave("prediction/shap-by-race.pdf", plot=p, width=10, height=15)



# Get the boostrapped AUROCs and AUPRCs at each timepoint

get_auroc <- function(preds) {
    idx <- sample(1:nrow(preds), size=nrow(preds), replace=TRUE)
    auroc <- yardstick::roc_auc_vec(truth= factor(preds$ys_true[idx], levels=c(0, 1)), estimate=preds$ys_pred[idx], case_weights=preds$weight[idx], event_level="second")
    return(auroc)
}

get_auprc <- function(preds) {
    idx <- sample(1:nrow(preds), size=nrow(preds), replace=TRUE)
    auprc <- yardstick::pr_auc_vec(truth= factor(preds$ys_true[idx], levels=c(0, 1)), estimate=preds$ys_pred[idx], case_weights=preds$weight[idx], event_level="second")
    return(auprc)
}


# bootstrapped_aurocs <- data.frame(timepoint=character(0), mean=numeric(0), lower_conf=numeric(0), upper_conf=numeric(0))
# bootstrapped_auprcs <- data.frame(timepoint=character(0), mean=numeric(0), lower_conf=numeric(0), upper_conf=numeric(0))


# index <- 1
# for (timepoint in c(0, 60, "last")) {
#     if (timepoint=="last") { #choose the last timepoint of each visit
#         preds <- predictions %>% group_by(csn) %>% mutate(last_timepoint=max(snapshot)) %>% 
#             ungroup() %>% filter(snapshot==last_timepoint)
#     }
#     else {
#         assert(timepoint %in% predictions$snapshot)
#         preds <- predictions %>% filter(snapshot==timepoint)
#         print(head(preds))
#     }
#     aurocs <- purrr::map_vec(1:1000, ~get_auroc(preds), .progress=TRUE)
#     auprcs <- purrr::map_vec(1:1000, ~get_auprc(preds), .progress=TRUE)

#     aurocs <- sort(aurocs)
#     bootstrapped_aurocs[index,] <- c(timepoint, mean(aurocs), aurocs[25], aurocs[975])
#     auprcs <- sort(auprcs)
#     bootstrapped_auprcs[index,] <- c(timepoint, mean(auprcs), auprcs[25], auprcs[975])
#     index <- index + 1
# }

# write.csv(bootstrapped_aurocs, "prediction/time_based_bootstrapped_aurocs.csv")
# write.csv(bootstrapped_auprcs, "prediction/time_based_bootstrapped_auprcs.csv")

# #Get bootstrapped AUROCs and AUPRCs by race

# bootstrapped_aurocs <- data.frame(race=character(0), mean=numeric(0), lower_conf=numeric(0), upper_conf=numeric(0))
# bootstrapped_auprcs <- data.frame(race=character(0), mean=numeric(0), lower_conf=numeric(0), upper_conf=numeric(0))

# #Attach predictions to race
# predictions <- predictions %>% inner_join(select(visits, c(csn, race)), by="csn")


# index <- 1
# for (race_ in c("Hispanic", "Non-Hispanic White", "Non-Hispanic Black")) {
#     preds <- predictions %>% filter(snapshot==60, race==race_)
#     print(head(preds))

#     aurocs <- purrr::map_vec(1:1000, ~get_auroc(preds), .progress=TRUE)
#     auprcs <- purrr::map_vec(1:1000, ~get_auprc(preds), .progress=TRUE)

#     print(length(aurocs))
#     print(unique(aurocs))
#     print(unique(preds$ys_true))


#     aurocs <- sort(aurocs)
#     bootstrapped_aurocs[index,] <- c(race_, mean(aurocs), aurocs[25], aurocs[975])
#     auprcs <- sort(auprcs)
#     bootstrapped_auprcs[index,] <- c(race_, mean(auprcs), auprcs[25], auprcs[975])
#     index <- index + 1
# }

# write.csv(bootstrapped_aurocs, "prediction/time_based_bootstrapped_aurocs_by_race_at_1_hour.csv")
# write.csv(bootstrapped_auprcs, "prediction/time_based_bootstrapped_auprcs_by_race_at_1_hour.csv")

