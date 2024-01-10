# Data cleaning
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(phyloseq)
library(tidyverse)
library(rio)
library(haven)

## Open phyloseq object
tab <- readRDS('data/phyloseq/rarefied/phyloseq_rarefied.RDS')
tax <- readRDS('data/phyloseq/rarefied/taxtable_rarefied.RDS')

## Inspect phyloseq object
ntaxa(tab)
nsamples(tab)
depth <- colSums(tab@otu_table)
all(depth == depth[1]) # TRUE -> rarefied

## Fix IDs
idlist <- str_remove(sample_names(tab), '_T1')
sample_names(tab) <- idlist
sample_names(tab)

### Get top 300 ASV by abundance 
ss <- taxa_sums(tab)
ss <- ss[order(ss, decreasing = T)]
ss <- ss[1:300]
top300 <- names(ss)
tax300 <- tax[rownames(tax) %in% top300, ]

sum(!is.na(tax300$Species)) / nrow(tax300) * 100 # 25.7 % of top 300 ASVs have species level
sum(!is.na(tax300$Genus)) / nrow(tax300) * 100 # 79 % of top 300 ASVs  have genus level
sum(!is.na(tax300$Family)) / nrow(tax300) * 100 # 97.7 % of top 300 ASVs have family level
sum(!is.na(tax300$Phylum)) / nrow(tax300) * 100 # 100 % of top 300 ASVs have phylum level

## Open HELIUS clinical data
df <- haven::read_sav("data/RODAM_dataset.sav") %>% 
    filter(!duplicated(RodamID))
reninaldo <- haven::read_sav("data/RODAM_reninaldo.sav") %>% 
    filter(!duplicated(RodamID))
df_compl <- full_join(df, reninaldo, by = "RodamID")

## Clean RODAM dataframe
df_new <- df_compl %>% 
    dplyr::select(ID=RodamID, Age=R2_Age, Sex=R2_Sex, Site=R2_Site,
                  Smoking=R1R2_Smoking, Alcohol=R2_alcohol_units_day, 
                  # Physical exam
                  BMI=R2_BMI, SBP=R2_BPsys_mean, 
                  DBP=R2_BPdia_mean, BPcat=R2_BPcat_ESC,
                  # Cardiometabolic
                  HT=R2_HTN_MedBP, Hyperchol=R1R2_CholDiagn, Stroke=R1R2_Stroke,
                  DM=R2_DM_Dichot, ACR_KDIGO=R2_ACR_KDIGO,
                  FramRisk=R2_FramScore,
                  # Lab
                  Na=R2_Na, UrineNa=R2_Na_Urine, Ka=R2_K, UrineK=R2_K_Urine,
                  GFR=R2_CKDEPI_eGFR_2021, ACR=R2_ACR,
                  TC = R2_Chol, LDL=R2_LDLChol, HDL=R2_HDLChol, 
                  Trig=R2_TG, HbA1c=R2_HbA1c, CRP=R2_CRP, 
                  Renin=R1_Renin, Aldo=R1_Aldosterone,
                  # Diet
                  TotalCalories=R2_TotalEnergy, Proteins=R2_PROT_day, Fat=R2_Fat_day,
                  Fibre=R2_Fiber_diet, Carbohydrates=R2_CHO, AlcoholIntake=R2_Alcohol_day,
                  SodiumInt=R2_Natrium_diet,
                  # Medication
                  DMMed=R2_DiabetesMeds,AB=R2_Antibacterials, AFung=R2_Antifungals,
                  LipidLowering=R2_Antilipidemics, 
                  # Antihypertensive medication
                  AntiHT=R2_AntihypertensivesC02, Diuretics=R2_AntihypertensivesC03,
                  CalciumAnt=R2_AntihypertensivesC08, BetaBlocker=R2_AntihypertensivesC07,
                  RAASi=R2_AntihypertensivesC09,
                  # Fecal sample 
                  FecalSample_AB=R2_Stool_AB, FecalSample_Prob=R2_Stool_Probiotics, 
                  FecalSample_Bristol=R2_Stool_BSS
    ) %>% 
    naniar::replace_with_na(., to_na = list("routing missing", "missing",                                      "nvt", "Routing missing", "Missing")) %>% 
    mutate(across(c("Sex", "Smoking", "Site", "FecalSample_AB", "FecalSample_Prob",
                    "FecalSample_Bristol",
                    "HT", "DM", "Hyperchol", "Stroke", "ACR_KDIGO",
                    "AntiHT", "Diuretics", "BPcat",
                    "CalciumAnt", "BetaBlocker", "RAASi", 
                    "LipidLowering","DMMed","AB", "AFung"),
                  ~as_factor(.x, levels = c("labels"))),
           CurrSmoking = fct_recode(Smoking, "Yes" = "Yes",
                                            "No" = "No, I have never smoked",
                                            "No" = "No, but I used to smoke"),
           CurrSmoking = fct_rev(CurrSmoking),
           FecalSample_AB = fct_recode(FecalSample_AB, "Yes" = "Yes",
                                    "No" = "No",
                                    "No" = "I do not know"),
           FecalSample_AB = fct_rev(FecalSample_AB),
           FecalSample_Prob = fct_rev(FecalSample_Prob),
           BristolScale = str_remove(FecalSample_Bristol, "Type "),
           BristolScale = as.numeric(as.character(BristolScale)),
           ARR = Aldo/Renin
    ) %>% 
    filter(Site %in% c("Urban Ghana", "Rural Ghana")) %>% 
    mutate(across(where(is.numeric), as.numeric) # all other vars to numeric, do this last
    ) %>% 
    # remove unused levels
    droplevels(.) %>% 
    filter(AB != "Yes") %>% 
    filter(BristolScale != 7)
    # filter(FecalSample_AB == "No") 205 of 1050 on antibiotics.. so skipped this step
dim(df_new)

## Select IDs present in both dataframes
sample_names_to_keep <- df_new$ID
tabprune <- prune_samples(sample_names_to_keep, tab)
df_new <- df_new[df_new$ID %in% sample_names(tab), ]
sampledata <- sample_data(df_new)
sample_names(sampledata) <- df_new$ID
rodamcomplete <- merge_phyloseq(tabprune, sampledata)

## Save files
saveRDS(df_new, file = "data/clinicaldata.RDS")
saveRDS(rodamcomplete, file = "data/phyloseq_sampledata.RDS")
saveRDS(tax, file = "data/taxtable.RDS")

