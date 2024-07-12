# Data cleaning
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(phyloseq)
library(tidyverse)
library(rio)
library(haven)

## Open phyloseq object
tab <- readRDS('data/rarefied/phyloseq_rarefied.RDS')
tax <- readRDS('data/rarefied/taxtable_rarefied.RDS')
ids <- haven::read_sav("data/HELIUS_participants_Total_ProsRod.sav")

## Inspect phyloseq object
ntaxa(tab)
nsamples(tab)
depth <- colSums(tab@otu_table)
all(depth == depth[1]) # TRUE -> rarefied

## Fix IDs
idlist <- str_remove(sample_names(tab), '_T1')
idlist <- case_when(str_detect(idlist, "X") ~ ids$RodamID[match(str_remove(idlist, "X"), ids$Heliusnr)],
          .default=idlist)
sample_names(tab) <- idlist

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
df2 <- haven::read_sav("data/ROD2_18_BV_Microbiome_HTN_20240226_RODAMbaseline_ProsRODAMfollowup_April2023_additionalVars.sav") %>% filter(!duplicated(RodamID))
reninaldo <- haven::read_sav("data/RODAM_reninaldo.sav") %>% 
    filter(!duplicated(RodamID))
df_compl <- full_join(df, reninaldo, by = "RodamID") %>% 
    full_join(., df2, by = "RodamID")

## Clean RODAM dataframe
na_strings <- c("","routing missing", "missing", 
                "nvt", "Routing missing", "Missing")
df_compl <- df_compl %>% 
    dplyr::select(ID=RodamID, Age=R2_Age, Sex=R2_Sex, Site=R2_Site,
                  Smoking=R1R2_Smoking, Alcohol=R2_alcohol_units_day,
                  Occupation_binary = R1R2_Occupation3, 
                  Occupation = R1R2_Occupation6,
                  Education = R1R2_Education,
                  Profession = R1R2_Profession,
                  # Living environment
                  LocGhana1 = R1_LiveGHA, LocGhana2 = R2_LiveGHA,
                  CityGhana1 = R1_LiveCityGHAX, CityGhana2 = R2_LiveCityGHAX,
                  VillGhana1 = R1_LiveVillageGHAX, VillGhana2 = R2_LiveVillageGHAX,
                  TimeLoc1 = R1_LiveTime, TimeLoc2 = R2_LiveTime,
                  LocEarlyLife = R2_EarlyLife, MoveEarlyLife = R2_EarlyLifeCityVil,
                  # Physical exam
                  BMI=R2_BMI, SBP=R2_BPsys_mean, 
                  DBP=R2_BPdia_mean, BPcat=R2_BPcat_ESC,
                  # Cardiometabolic
                  HT=R2_HTN_MedBP, Hyperchol=R1R2_CholDiagn, Stroke=R1R2_Stroke,
                  DM=R2_DM_Dichot, ACR_KDIGO=R2_ACR_KDIGO,
                  FramRisk=R2_FramScore, Active=R2_Active, PhysAct=R2_Ptotallevels,
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
    naniar::replace_with_na_all(condition = ~.x %in% na_strings)

df_new <- df_compl %>% 
    mutate(across(c("VillGhana1", "VillGhana2", "CityGhana1", "CityGhana2"), str_to_lower),
           VillGhana = case_when(
               !is.na(VillGhana1) ~ VillGhana1,
               !is.na(VillGhana2) ~ VillGhana2
           ),
           CityGhana = case_when(
               !is.na(CityGhana1) ~ CityGhana1,
               !is.na(CityGhana2) ~ CityGhana2
           ),
            across(c("Sex", "Smoking", "Site", "FecalSample_AB", "FecalSample_Prob",
                    "FecalSample_Bristol",
                    "HT", "DM", "Hyperchol", "Stroke", "ACR_KDIGO",
                    "AntiHT", "Diuretics", "BPcat",
                    "CalciumAnt", "BetaBlocker", "RAASi", 
                    "LipidLowering","DMMed","AB", "AFung", "LocGhana1", "LocGhana2",
                    "Occupation", "Occupation_binary", "Profession", "Education",
                    "LocEarlyLife", "MoveEarlyLife", "Active"),
                  ~as_factor(.x, levels = c("labels"))),
           LocGhana = case_when(
               !is.na(LocGhana1) ~ LocGhana1,
               !is.na(LocGhana2) ~ LocGhana2
           ),
           TimeLoc = case_when(
               !is.na(TimeLoc1) ~ TimeLoc1,
               !is.na(TimeLoc2) ~ TimeLoc2
           ),
           MoveEarlyLife = case_when(
               MoveEarlyLife == "I first lived in a village and moved to a city later on" 
                    ~ paste("VillageToCity"),
               MoveEarlyLife == "I first lived in a city and moved to a village later on" 
                    ~ paste("CityToVillage"),
               MoveEarlyLife == "I moved between village and city several times" 
                    ~ paste("MultipleMoves")
           ),
           # across(c("LocGhana", "MoveEarlyLife", "VillGhana", "CityGhana"), as.factor),
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
           AlcoholBin = case_when(Alcohol > 0 ~ "Alcohol use",
                                          Alcohol == 0 ~ "No alcohol use"),
           AlcoholBin = fct_rev(AlcoholBin),
           PhysActDay = PhysAct / 7,
           ARR = Aldo/Renin,
           AntiHT = case_when(
               AntiHT == "Yes" | Diuretics == "Yes" | CalciumAnt == "Yes" |
                   BetaBlocker == "Yes" | RAASi == "Yes" ~ paste0("Yes"),
               .default = "No"
           ),
           AntiHT = as.factor(AntiHT),
           AgeCat = case_when(Age < 50 ~ "Young", Age >= 50 ~ "Old"),
           AgeCat = as.factor(AgeCat),
           BMIcat = case_when(
               BMI < 20 ~ "<20",
               BMI >20 & BMI < 25 ~ "20-25",
               BMI > 25 & BMI < 30 ~ "25-30",
               BMI > 30 ~ ">30"
           )
    ) %>% 
    filter(Site %in% c("Urban Ghana", "Rural Ghana", "Amsterdam")) %>%
    filter(ID %in% sample_names(tab)) %>% 
    mutate(across(where(is.numeric), as.numeric) # all other vars to numeric, do this last
    ) %>% 
    # remove unused levels
    droplevels(.) %>% 
    select(-TimeLoc1, -TimeLoc2, -LocGhana1, -LocGhana2, -VillGhana1, -VillGhana2, 
           -CityGhana1, -CityGhana2) %>% 
    filter(AB != "Yes")

## Clean dietary data - two IDs have way too high intake
dietids <- c("G0421", "G2545")
df_new <- df_new %>% mutate(
    TotalCalories = case_when(ID %in% dietids ~ NA, .default = TotalCalories),
    Proteins = case_when(ID %in% dietids ~ NA, .default = Proteins),
    Fat = case_when(ID %in% dietids ~ NA, .default = Fat),
    Fibre = case_when(ID %in% dietids ~ NA, .default = Fibre),
    Carbohydrates = case_when(ID %in% dietids ~ NA, .default = Carbohydrates),
    SodiumInt = case_when(ID %in% dietids ~ NA, .default = SodiumInt),
    AlcoholIntake = log(AlcoholIntake+0.01)
)
    
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

