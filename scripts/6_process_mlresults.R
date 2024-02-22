## Process results machine learning

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(stringr)
library(forcats)

source("scripts/functions.R")

path_true <- 'rural_urban/output_XGB_class_rural_urban_2024_02_05__21-24-45'
path_permuted <- 'rural_urban/output_XGB_class_rural_urban_2024_02_05__21-37-26_PERMUTED'
data_path <- 'rural_urban/input_data'
labels <- c("Urban", "Rural")

plot_feature_importance_class(path_true, 20)
plot_feature_importance_color_microbiome(path_true, 20)
plot_features_tests_class(data_path, path_true, top_n=20, labels)
plot_features_tests_top(data_path, path_true, top_n=20, nrow=4, labels)

path_true <- 'urban_ams/output_XGB_class_urban_ams_2024_02_21__14-01-19'
path_permuted <- 'urban_ams/output_XGB_class_urban_ams_2024_02_21__14-10-44_PERMUTED'
data_path <- 'urban_ams/input_data'
labels <- c("Amsterdam", "Urban Ghana")

plot_feature_importance_class(path_true, 20)
plot_feature_importance_color_microbiome(path_true, 20)
plot_features_tests_class(data_path, path_true, top_n=20, labels)
plot_features_tests_top(data_path, path_true, top_n=20, nrow=4, labels)

path_true <- 'salt/output_XGB_reg_salt_2024_01_15__00-08-54'
# path_permuted <- 'salt/output_XGB_class_menopause_2023_12_21__20-51-36_PERMUTED'
data_path <- 'salt/input_data'

# compared_to_permuted_class(path_true, path_permuted)
plot_feature_importance_class(path_true, 20)
plot_feature_importance_color_microbiome(path_true, 20)
plot_features_tests_reg(data_path, path_true, top_n=20, 
                        outcome_name = "Sodium intake (mg)",
                        x_lab = "clr(counts)")
plot_features_top_n_reg(data_path, path_true, top_n=20, nrow=4, 
                        outcome_name = "Sodium intake (mg)",
                        x_lab = "clr(counts)")
