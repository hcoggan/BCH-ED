
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


#Make a descriptive table of the Stanford cohort, including admission rates.

setwd("/Volumes/chip-lacava/Groups/BCH-ED/reprocessing")
visits <- read.csv("preprocessed-visits.csv")


vars_to_include <- c("sex", "race", "age_group", "ed_arrival_mode", "insurance",  "triage_acuity", "language")

all_tabs <- data.frame()

for (i in 1:length(vars_to_include)) {
    tab <- visits %>% group_by(.data[[vars_to_include[i]]]) %>% 
        summarise(count=n(), admission_rate=signif(sum(is_admitted)*100/count, 3)) %>% ungroup() %>% mutate(perc=signif(count*100/sum(count), 3)) %>%
        mutate(Number=paste0(count, " (", perc, ")")) %>% rename(Characteristic= !!vars_to_include[i]) %>% 
        select(Characteristic, Number, admission_rate)
    all_tabs <- rbind(all_tabs, tab)
    

}

print(paste("Admission rate:", sum(visits$is_admitted)/nrow(visits)))

write.csv(all_tabs, "visit_characteristics.csv")


print(paste("ED LOS", quantile(visits$ed_los, probs=c(0.25, 0.5, 0.75))))
print(paste("SDI", quantile(visits$sdi_score, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)))
print(paste("Miles travelled", quantile(visits$miles_travelled, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)))