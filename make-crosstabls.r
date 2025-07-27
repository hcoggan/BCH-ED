
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

#make a table of crosstabs to describe patient characteristics
setwd("your/working/directory")
 save_filepath <- "your/save/filepath"


baseline_factors <- read.csv("baseline-factors-with-top-200-complaint-stems-tagged.csv", encoding = "UTF-8") %>% 
    filter(is_admitted==1 | is_discharged==1) %>%
    select(-c(ed_arrival_time, ed_checkout_time)) %>%
     #include only those for whom a decision was made- admitted or discharged
    mutate(insurance=ifelse(has_medicaid==1, "Medicaid", "Non-Medicaid"), disposition=ifelse(is_admitted==1, "Admitted", "Discharged"))


###COLUMNS: RAW NUMBER OF VISITS (and percentage of total), ADMISSION RATE

print("Raw admission rate")
print(sum(baseline_factors$is_admitted)/nrow(baseline_factors))

vars_to_include <- c("sex", "race", "age_group", "ed_arrival_mode", "insurance", "primary_language", "triage_acuity")

all_tabs <- data.frame()

for (i in 1:length(vars_to_include)) {
    tab <- baseline_factors %>% group_by(.data[[vars_to_include[i]]]) %>% 
        summarise(count=n(), admission_rate=signif(sum(is_admitted)*100/count, 3)) %>% ungroup() %>% mutate(perc=signif(count*100/sum(count), 3)) %>%
        mutate(label=paste0(count, " (", perc, ")")) %>% rename(Characteristic= !!vars_to_include[i]) %>% 
        select(Characteristic, label, admission_rate)
    all_tabs <- rbind(all_tabs, tab)
    

}


write.csv(all_tabs, paste0(save_filepath, "visit_characteristics.csv"))
