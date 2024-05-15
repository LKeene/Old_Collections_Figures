# This script can be used to perform the multiple linear regression analysis on 
# the RNA concentration and RNA length for experimentally dried and frozen insects.

library(tidyverse)
library(emmeans)
library(performance)
library(broom)
library(car)

# Read in the data
conc <- read_csv("metadata/OverTimeRNAConcentrations.csv")
fly_length <- read_csv("metadata/fly_length_data.csv")
mos_length <- read_csv("metadata/mos_length_data.csv")

# Cleanup & convert predictors to factors
conc <- conc %>% 
  select(week, storage, organism, concentration) %>% 
  mutate(week = as.factor(week),
         storage = as.factor(storage))

# Fly concentration summary table
conc %>% 
  filter(organism == "fly") %>% 
  select(week, storage, concentration) %>% 
  group_by(week, storage) %>% 
  summarize(n = n(),
            Mean = mean(concentration),
            SD = sd(concentration),
            SE_mean = SD/sqrt(n)) 

# Mosquito concentration summary table
conc %>% 
  filter(organism == "mosquito") %>% 
  select(week, storage, concentration) %>% 
  group_by(week, storage) %>% 
  summarize(n = n(),
            Mean = mean(concentration),
            SD = sd(concentration),
            SE_mean = SD/sqrt(n))

fly_length <- fly_length %>% 
  select(weeks, group, mean_length, sample_name) %>% 
  group_by(weeks, group) %>% 
  filter(group != "old") %>% 
  mutate(weeks = as.factor(weeks),
         group = as.factor(group))

# Fly length summary table
fly_length %>% 
  select(weeks, group, mean_length) %>% 
  group_by(weeks, group) %>% 
  summarize(n = n(),
            Mean = mean(mean_length),
            SD = sd(mean_length),
            SE_mean = SD/sqrt(n))

mos_length <- mos_length %>% 
  select(weeks, group, mean_length, sample_name) %>% 
  group_by(weeks, group) %>% 
  mutate(weeks = as.factor(weeks),
         group = as.factor(group))

# Mosquito length summary table
mos_length %>% 
  select(weeks, group, mean_length) %>% 
  group_by(weeks, group) %>% 
  summarize(n = n(),
            Mean = mean(mean_length),
            SD = sd(mean_length),
            SE_mean = SD/sqrt(n)) 

# Get averages and standard deviations
conc_mean <- conc %>% 
  group_by(week, storage, organism) %>% 
  mutate(mean_conc = mean(concentration),
         sd_conc = sd(concentration))

df_mean_fly <- fly_length %>% 
  group_by(group, weeks) %>% 
  mutate(mean_length_g = mean(mean_length),
         sd_mean_length_g = sd(mean_length))

df_mean_mos <- mos_length %>% 
  group_by(group, weeks) %>% 
  mutate(mean_length_g = mean(mean_length),
         sd_mean_length_g = sd(mean_length))

## Fly Concentration
fly_conc <- conc_mean %>% 
  filter(organism == "fly")

# Fly concentration One-Way Model
conc_fly_mlr2 <- lm(concentration~week+storage, data = fly_conc)
#tidy(conc_fly_mlr2)
Anova(conc_fly_mlr2)

check_model(conc_fly_mlr2, check = c("qq", "linearity", "homogeneity"))

emm_fly_conc <- emmeans(conc_fly_mlr2, ~week)
pwpp(emm_fly_conc)

## Mosquito Concentration
mos_conc <- conc_mean %>% 
  filter(organism == "mosquito")

# Mosquito concentration One-Way Model
conc_mos_mlr2 <- lm(concentration~week+storage, data = mos_conc)
#tidy(conc_mos_mlr2)
Anova(conc_mos_mlr2)

check_model(conc_mos_mlr2, check = c("qq", "linearity", "homogeneity"))

emm_mos_conc <- emmeans(conc_mos_mlr2, ~week)
pwpp(emm_mos_conc)

df_mean_fly <- df_mean_fly %>% 
  filter(group != "Fresh")

# Fly length One-Way Model
length_fly_mlr2 <- lm(mean_length_g~weeks+group, data = df_mean_fly)
#tidy(length_fly_mlr2)
Anova(length_fly_mlr2)

check_model(length_fly_mlr2, check = c("qq", "linearity", "homogeneity"))

emm_fly_leng <- emmeans(length_fly_mlr2, ~weeks)
pwpp(emm_fly_leng)

# Mosquito length One-Way Model
df_mean_mos <- df_mean_mos %>% 
  filter(group != "Fresh")
length_mos_mlr2 <- lm(mean_length_g~weeks+group, data = df_mean_mos)
#tidy(length_mos_mlr2)
Anova(length_mos_mlr2)

check_model(length_mos_mlr2, check = c("qq", "linearity", "homogeneity"))

emm_mos_leng <- emmeans(length_mos_mlr2, ~weeks)
pwpp(emm_mos_leng)

