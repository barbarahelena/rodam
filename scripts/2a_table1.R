#### RODAM characteristics

## Libraries
library(tableone)
library(dplyr)

## Output folder
resultsfolder <- "results/table1"
dir.create(resultsfolder, showWarnings = FALSE)

## Load dataset
df_new <- rio:: import("data/clinicaldata.RDS")

## Table 1
Table1 <- CreateTableOne(data = df_new, 
                         vars = c("Age", "Sex", "BMI", "CurrSmoking",
                                  "SBP", "DBP", "TC", "LDL", "HbA1c", "GFR", "CRP",
                                  "HT", "AntiHT", "DM", "DMMed", "LipidLowering",
                                  "FecalSample_AB", "FecalSample_Prob", "BristolScale",
                                  "TotalCalories", "Carbohydrates", "Proteins",  
                                  "Fat", "Fibre", "SodiumInt", "Alcohol", "AlcoholBin",
                                  "PhysAct", "Occupation_binary"), 
                         strata = c("Site"), addOverall = TRUE)
Table1 <- print(Table1, contDigits = 1, missing = TRUE, nonnormal = c("BristolScale", "PhysAct"))
Table1 <- as.data.frame(Table1)
table1 <- Table1 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(table1, "results/table1/table1.csv")
