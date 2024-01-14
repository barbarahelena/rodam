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
                                  "FecalSample_AB", "FecalSample_Prob", "BristolScale",
                                  "Renin", "Aldo", "ARR"), 
                         strata = c("Site"))
Table1 <- print(Table1, contDigits = 1, missing = TRUE, nonnormal = "BristolScale")
Table1 <- as.data.frame(Table1)
table1 <- Table1 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(table1, "results/table1/table1.csv")

## PCA diet
df_diet <- df_new %>% select(ID, Site, TotalCalories, Fibre, Proteins, 
                             Carbohydrates, Fat, SodiumInt, AlcoholIntake) %>% 
    filter(!is.na(TotalCalories)) %>% 
    filter(TotalCalories < 7500)
df_diet2 <- df_diet %>% select(-ID, -Site)
matdiet <- as.matrix(df_diet2)
pc <- mixOmics::pca(matdiet)
pc2 <- as.data.frame(pc$variates$X)
pc2 <- pc2 %>% mutate(ID = df_diet$ID, Site = df_diet$Site)
expvar_diet <- pc$prop_expl_var$X[1:2]
loadings <- as.data.frame(pc$loadings$X)
loadings$Variables <- rownames(loadings)
loadings <- loadings %>% filter(Variables %in% c("TotalCalories", "SodiumInt")) %>% 
    mutate(Variables = fct_recode(Variables, "Sodium" = "SodiumInt"))
(pcadiet <- pc2 %>% 
        ggplot(aes(PC1, PC2)) +
        geom_point(aes(color = Site), size = 1, alpha = 0.7) +
        xlab(paste0('PC1 (', round(expvar_diet[1]*100, digits = 1),'%)')) +
        ylab(paste0('PC2 (', round(expvar_diet[2]*100, digits = 1),'%)')) +
        theme_Publication() +
        stat_ellipse(geom = "polygon", aes(color = Site, fill = Site), alpha = 0.2, type = "norm")+
        scale_color_manual(values = pal_futurama()(4)[c(3,4)]) +
        scale_fill_manual(values = pal_futurama()(4)[c(3,4)], guide = "none") +
        labs(color = "", title = "PCA diet")+
        geom_segment(data = loadings, aes(x = 0, y = 0, xend = (PC1*2500), yend = (PC2*2500)), 
                     arrow = arrow(length = unit(1/2, "picas")),
                     color = "black", linewidth = 0.9) +
        annotate("text", x = (loadings$PC1*2700), y = (loadings$PC2*2700),
                 label = loadings$Variables)
        )
ggsave(pcadiet, filename = "results/table1/PCA_diet_loading.pdf", device = "pdf", width = 5, height = 5)
ggsave(pcadiet, filename = "results/table1/PCA_diet_loading.svg", device = "svg", width = 5, height = 5)

df <- left_join(df_new, pc2, by = c("ID", "Site")) %>% select(everything(.), DietPC1=PC1, DietPC2=PC2)
saveRDS(df, "data/clinicaldata_pcdiet.RDS")

pl1 <- ggplot(df_diet, aes(x=Site, y=TotalCalories))+
    geom_violin(aes(fill=Site))+
    scale_fill_manual(values = pal_futurama()(4)[3:4], guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'Kcal', title = "Total calories")+
    ggpubr::stat_compare_means()
pl1

pl2 <- ggplot(df_diet, aes(x=Site, y=Fibre))+
    geom_violin(aes(fill=Site))+
    scale_fill_manual(values = pal_futurama()(4)[3:4], guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'grams', title = "Fibre")+
    ggpubr::stat_compare_means()
pl2

pl3 <- ggplot(df_diet, aes(x=Site, y=Proteins))+
    geom_violin(aes(fill=Site))+
    scale_fill_manual(values = pal_futurama()(4)[3:4], guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'grams', title = "Proteins")+
    ggpubr::stat_compare_means()
pl3

pl4 <- ggplot(df_diet, aes(x=Site, y=Carbohydrates))+
    geom_violin(aes(fill=Site))+
    scale_fill_manual(values = pal_futurama()(4)[3:4], guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'grams', title = "Carbohydrates")+
    ggpubr::stat_compare_means()
pl4

pl5 <- ggplot(df_diet, aes(x=Site, y=Fat))+
    geom_violin(aes(fill=Site))+
    scale_fill_manual(values = pal_futurama()(4)[3:4], guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'grams', title = "Fat")+
    ggpubr::stat_compare_means()
pl5

pl6 <- ggplot(df_diet, aes(x=Site, y=SodiumInt))+
    geom_violin(aes(fill=Site))+
    scale_fill_manual(values = pal_futurama()(4)[3:4], guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    theme(legend.position = 'none')+
    labs(x='', y = 'mg', title = "Sodium")+
    ggpubr::stat_compare_means()
pl6

pl7 <- ggplot(df_diet, aes(x=Site, y=(AlcoholIntake+0.01)))+
    geom_violin(aes(fill=Site))+
    scale_fill_manual(values = pal_futurama()(4)[3:4], guide = FALSE)+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+
    theme_Publication()+
    scale_y_log10()+
    theme(legend.position = 'none')+
    labs(x='', y = 'log10(grams)', title = "Alcohol")+
    ggpubr::stat_compare_means()
pl7

ggarrange(pl1, pl2, pl3, pl4, pl5, pl6, pl7, nrow = 2, ncol = 4,
          labels = c(LETTERS[1:7]))
ggsave("results/table1/dietaryintake.pdf", width = 12, height = 8)
ggsave("results/table1/dietaryintake.svg", width = 12, height = 8)

