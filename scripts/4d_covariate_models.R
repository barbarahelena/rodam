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
                )
    # print(df$var)
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
    mutate(
        Women = case_when(
            Sex == "Female" ~ "Yes",
            Sex == "Male" ~ "No"
        ),
        Active_binary = case_when(
            Active == "Low level" ~ "No",
            Active == "Moderate level" | Active == "High level" ~ "Yes"
        ),
        Occupation_manual = case_when(
            Occupation_binary == "manual" ~ "Yes",
            Occupation_binary == "non manual" ~ "No"
        ),
        EarlyLifeUrban_binary = case_when(
            LocEarlyLife == "I lived in a village (rural)" ~ "No",
            LocEarlyLife == "I lived in a city (urban)" ~ "Yes",
            LocEarlyLife == "I lived both in a city and a village" ~ "Yes"
        )
    )
mb <- readRDS('data/phyloseq_sampledata.RDS')

# Continuous variables (regression models)
make_input_folder_cont(df, mb, "Age", "age")
make_input_folder_cont(df, mb, "BMI", "bmi")
make_input_folder_cont(df, mb, "TotalCalories", "kcal")
make_input_folder_cont(df, mb, "Carbohydrates", "carbohydrates")
make_input_folder_cont(df, mb, "Proteins", "proteins")
make_input_folder_cont(df, mb, "Fibre", "fibre")
make_input_folder_cont(df, mb, "Fat", "fat")
make_input_folder_cont(df, mb, "SodiumInt", "sodium")
make_input_folder_cont(df, mb, "CRP", "crp")
make_input_folder_cont(df, mb, "LDL", "ldl")
make_input_folder_cont(df, mb, "GFR", "gfr")
make_input_folder_cont(df, mb, "BristolScale", "bristol")
make_input_folder_cont(df, mb, "PhysAct", "physicalactivity")

# Binary variables (classification models)
make_input_folder_bin(df, mb, "Women", "women")
make_input_folder_bin(df, mb, "FecalSample_Prob", "probiotics")
make_input_folder_bin(df, mb, "AntiHT", "bpmed")
make_input_folder_bin(df, mb, "Occupation_manual", "occupationmanual")
make_input_folder_bin(df, mb, "Active_binary", "active")
make_input_folder_bin(df, mb, "EarlyLifeUrban_binary", "earlylifeurban")
