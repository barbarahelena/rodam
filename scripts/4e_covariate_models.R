# create XGB input files for covariates

library(dplyr)
library(phyloseq)
rm(list=ls())

make_input_folder_cont <- function(clindata, otutable, namevar, namefolder){
    write_data <- function(x, data_path){
        x <- as.matrix(x)
        if(any(is.na(x))){
            cat('There are missing values in the input data!\n')
        }
        write.table(x, file.path(data_path, 'X_data.txt'), row.names = F, col.names = F, 
                    sep = '\t', quote = F)
        write.table(colnames(x), file.path(data_path,'feat_ids.txt'), row.names = F, 
                    col.names = F, sep = '\t', quote = F)
        write.table(rownames(x), file.path(data_path,'subject_ids.txt'), row.names = F, 
                    col.names = F, sep = '\t', quote = F)
    }
    write_y <- function(x, name_y, data_path){
        if(missing(name_y)){
            cat('\n\nYou need to provide a name for the y data file!\n')
        }
        if(!name_y %in% c('y_binary.txt', 'y_reg.txt')){
            cat('\nThe file name is not compatible with XGBeast!\n' )
        }
        if(any(is.na(x))){
            cat('\nThere are missing values in the outcome data!\n')
        }
        write.table(x, file = file.path(data_path, name_y), row.names = F, col.names = F, sep = '\t', quote = F)
    }
    df$var <- df[[namevar]]
    df <- df %>% dplyr::select(ID, var) %>% filter(!is.na(var))
    mb <- prune_samples(sample_names(mb) %in% df$ID, mb)
    otu <- t(as(mb@otu_table, "matrix"))
    tk <- apply(otu, 2, function(x) sum(x > 10) > (0.3*length(x)))
    mbdf <- otu[,tk]
    mbdf <- as.data.frame(mbdf)
    df <- df[match(rownames(mbdf), df$ID), ]
    print("all IDs in same sequence for predictors and outcome:")
    print(all(df$ID == rownames(mbdf))) # TRUE
    print("make folder for input data ML model")
    path <- namefolder
    dir.create(path)
    dir.create(file.path(path, "input_data"))
    write_data(mbdf, file.path(path, 'input_data'))
    y <- as.data.frame(df$var)
    y
    write_y(y, name_y = 'y_reg.txt', file.path(path, 'input_data'))
}

make_input_folder_bin <- function(clindata, otutable, namevar, namefolder){
    write_data <- function(x, data_path){
        x <- as.matrix(x)
        if(any(is.na(x))){
            cat('There are missing values in the input data!\n')
        }
        write.table(x, file.path(data_path, 'X_data.txt'), row.names = F, col.names = F, 
                    sep = '\t', quote = F)
        write.table(colnames(x), file.path(data_path,'feat_ids.txt'), row.names = F, 
                    col.names = F, sep = '\t', quote = F)
        write.table(rownames(x), file.path(data_path,'subject_ids.txt'), row.names = F, 
                    col.names = F, sep = '\t', quote = F)
    }
    write_y <- function(x, name_y, data_path){
        if(missing(name_y)){
            cat('\n\nYou need to provide a name for the y data file!\n')
        }
        if(!name_y %in% c('y_binary.txt', 'y_reg.txt')){
            cat('\nThe file name is not compatible with XGBeast!\n' )
        }
        if(any(is.na(x))){
            cat('\nThere are missing values in the outcome data!\n')
        }
        write.table(x, file = file.path(data_path, name_y), row.names = F, col.names = F, sep = '\t', quote = F)
    }
    df$var <- df[[namevar]]
    df <- df %>% dplyr::select(ID, var) %>% filter(!is.na(var)) %>% 
                mutate(var = case_when(
                    var == "Yes" ~ 1,
                    var == "No" ~ 0
                    )
                ) %>% filter(!is.na(var))
    print(df$var)
    mb <- prune_samples(sample_names(mb) %in% df$ID, mb)
    otu <- t(as(mb@otu_table, "matrix"))
    tk <- apply(otu, 2, function(x) sum(x > 10) > (0.3*length(x)))
    mbdf <- otu[,tk]
    mbdf <- as.data.frame(mbdf)
    df <- df[match(rownames(mbdf), df$ID), ]
    print("all IDs in same sequence for predictors and outcome:")
    print(all(df$ID == rownames(mbdf))) # TRUE
    print("make folder for input data ML model")
    path <- namefolder
    dir.create(path)
    dir.create(file.path(path, "input_data"))
    write_data(mbdf, file.path(path, 'input_data'))
    y <- as.data.frame(df$var)
    y
    write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))
}

## Open dataframe
df <- readRDS('data/clinicaldata.RDS') %>% 
    filter(Site == "Amsterdam") %>% 
    mutate(
        Women = case_when(
            Sex == "Female" ~ "Yes",
            Sex == "Male" ~ "No"
        ),
        Active_binary = case_when(
            Active == "Low level" ~ "No",
            Active == "Moderate level" | Active == "High level" ~ "Yes"
        ),
        EarlyLifeUrban_binary = case_when(
            LocEarlyLife == "I lived in a village (rural)" ~ "No",
            LocEarlyLife == "I lived in a city (urban)" ~ "Yes",
            LocEarlyLife == "I lived both in a city and a village" ~ "Yes"
        )
    )
mb <- readRDS('data/phyloseq_sampledata.RDS')

# Continuous variables (regression models)
make_input_folder_cont(df, mb, "Age", "covariatemodels_amsterdam/age")
make_input_folder_cont(df, mb, "BMI", "covariatemodels_amsterdam/bmi")
make_input_folder_cont(df, mb, "TotalCalories", "covariatemodels_amsterdam/kcal")
make_input_folder_cont(df, mb, "Carbohydrates", "covariatemodels_amsterdam/carbohydrates")
make_input_folder_cont(df, mb, "Proteins", "covariatemodels_amsterdam/proteins")
make_input_folder_cont(df, mb, "Fibre", "covariatemodels_amsterdam/fibre")
make_input_folder_cont(df, mb, "Fat", "covariatemodels_amsterdam/fat")
make_input_folder_cont(df, mb, "SodiumInt", "covariatemodels_amsterdam/sodium")
make_input_folder_cont(df, mb, "CRP", "covariatemodels_amsterdam/crp")
make_input_folder_cont(df, mb, "LDL", "covariatemodels_amsterdam/ldl")
make_input_folder_cont(df, mb, "GFR", "covariatemodels_amsterdam/gfr")
make_input_folder_cont(df, mb, "PhysAct", "covariatemodels_amsterdam/physicalactivity")
make_input_folder_cont(df, mb, "FramRisk", "covariatemodels_amsterdam/framingham")

# Binary variables (classification models)
make_input_folder_bin(df, mb, "Women", "covariatemodels_amsterdam/women")
make_input_folder_bin(df, mb, "HT", "covariatemodels_amsterdam/hypertension")
make_input_folder_bin(df, mb, "DM", "covariatemodels_amsterdam/diabetes")
make_input_folder_bin(df, mb, "Hyperchol", "covariatemodels_amsterdam/hyperchol")
make_input_folder_bin(df, mb, "AntiHT", "covariatemodels_amsterdam/bpmed")
make_input_folder_bin(df, mb, "Active_binary", "covariatemodels_amsterdam/active")
make_input_folder_bin(df, mb, "EarlyLifeUrban_binary", "covariatemodels_amsterdam/earlylifeurban")

