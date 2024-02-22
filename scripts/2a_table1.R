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
                                  "AlcoholIntake", "TotalCalories", "Fibre",
                                  "SBP", "DBP", "TC", "LDL", "HbA1c",
                                  "HT", "AntiHT", "DM", "DMMed", "LipidLowering",
                                  "FecalSample_AB", "FecalSample_Prob", "BristolScale"), 
                         strata = c("Site"))
Table1 <- print(Table1, contDigits = 1, missing = TRUE, nonnormal = "BristolScale")
Table1 <- as.data.frame(Table1)
table1 <- Table1 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(table1, "results/table1/table1.csv")
