## Mediation models
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(rio)
library(phyloseq)
library(dplyr)
library(tidyr)
library(forcats)
library(patchwork)

## Output folder
resultsfolder <- "results/glm_clr"
dir.create(resultsfolder, showWarnings = FALSE)

## Functions
afronden2 <- function(x) return(as.numeric(format(round(x, 2),2)))
afronden3 <- function(x) return(as.numeric(format(round(x, 3),3)))

## Opening RODAM file
rodam_mb <- readRDS("data/phyloseq_sampledata.RDS")
rodamclin <- readRDS("data/clinicaldata_pcdiet.RDS")
rodamotu <- as.data.frame(t(as(rodam_mb@otu_table, "matrix")))
tax <- readRDS("data/tax_table.RDS")
best_pred <- rio::import('rural_urban/output_XGB_class_rural_urban_2024_01_14__00-47-28/feature_importance.txt')

## Preparation for models
best_pred <- best_pred %>% arrange(-RelFeatImp) %>% slice(1:20)
dfmb <- rodamotu %>% select(best_pred$FeatName)
dfmb2 <- as.data.frame(clr(dfmb+1))
dfmb2$ID <- rownames(dfmb2)
df <- left_join(rodamclin, dfmb2, by='ID')
dim(df)
colnames(df)[(ncol(df)-19):ncol(df)] <- make.unique(tax$Tax[match(colnames(df)[(ncol(df)-19):ncol(df)], tax$ASV)])
dfmb2 <- df[(ncol(df)-19):ncol(df)]

## Mediation of diet
## Path A: location and health factors
res_health <- c()
for (i in c(7:9,28)) {
    df$health <- scale(df[[i]])
    m0 <- lm(health ~ Site, data = df)
    dietfactor <- colnames(df)[i]
    m0 <- tidy(m0, conf.int=T)[2,]
    resRow <- cbind(dietfactor, m0$estimate, m0$conf.low, m0$conf.high, m0$p.value)
    colnames(resRow) <- c("Healthfactor", "m0-est", "m0-l95", "m0-u95", "m0-p")
    res_health <- rbind(res_health, resRow)
    df$health <- NULL
}
res_health <- as.data.frame(res_health)
res_health2 <- res_health %>% 
    mutate_at(c(2:5), as.character) %>% 
    mutate_at(c(2:5), as.numeric) %>% 
    mutate_at(c(2:4), afronden2) %>% 
    mutate_at(c(5), afronden3)
openxlsx::write.xlsx(res_health2, "results/glm_clr/lm_health.xlsx")

res_health3 <- res_health2 %>% 
    pivot_longer(c(2:5), names_to=c("model", "cat"), 
                 names_prefix="m", 
                 names_sep='-',
                 values_to="value") %>% 
    pivot_wider(names_from = cat, values_from = value) %>% 
    mutate(model = factor(model, levels = c("0"), 
                          labels = c("Unadjusted")),
           Healthfactor = factor(Healthfactor, levels = colnames(df)[c(7:9,28)]),
           Healthfactor = fct_rev(Healthfactor))

(plhealth <- ggplot(res_health3, aes(x=Healthfactor, y=est)) +
        geom_hline(yintercept = 0, color = "grey40") +
        geom_point(position=position_dodge(-0.7)) +
        scale_shape_manual(values = 19)+
        # scale_y_continuous(breaks = c(-2:8))+
        geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.7)) +
        theme_Publication()+
        theme(legend.position = "bottom")+
        labs(title = "Dietary factors and location",
             x = "", y = "Difference (SD) for urban location", shape = "", color = "") +
        scale_color_cosmic() +
        coord_flip())
ggsave("results/glm_clr/lm_health.pdf", width = 4, height = 3)
ggsave("results/glm_clr/lm_health.svg", width = 4, height = 3)

## Path B: diet and ASVs
res_asvs <- c()
for (i in c(7)) {
    df$diet <- scale(df[[i]])
    dietfactor <- colnames(df)[i]
    for(j in c(55:74)){
        df$asv <- df[[j]]
        taxasv <- colnames(df)[j]
        m0 <- lm(asv ~ diet, data = df)
        m0 <- tidy(m0, conf.int=T)[2,]
        resRow <- cbind(dietfactor, taxasv, m0$estimate, m0$conf.low, m0$conf.high, m0$p.value)
        colnames(resRow) <- c("Diet", "ASV", "m0-est", "m0-l95", "m0-u95", "m0-p")
        res_asvs <- rbind(res_asvs, resRow)
    }
    df$diet <- NULL
}
res_asvs <- as.data.frame(res_asvs)
res_asvs2 <- res_asvs %>% 
    mutate_at(c(3:6), as.character) %>% 
    mutate_at(c(3:6), as.numeric) %>% 
    mutate_at(c(3:5), afronden2) %>% 
    mutate(
        `m0-q` = p.adjust(`m0-p`, 'fdr')
    ) %>% 
    mutate_at(c(6), afronden3) 

openxlsx::write.xlsx(res_asvs2, "results/glm_clr/lm_asvsdiet.xlsx")

res_asvs3 <- res_asvs2 %>% 
    pivot_longer(c(3:7), names_to=c("model", "cat"), 
                 names_prefix="m", 
                 names_sep='-',
                 values_to="value") %>% 
    pivot_wider(names_from = cat, values_from = value) %>% 
    mutate(model = factor(model, levels = c("0"), 
                          labels = c("Unadjusted")),
           ASV = factor(ASV, levels = colnames(df)[55:74]),
           ASV = fct_rev(ASV),
           Diet = factor(Diet, levels = colnames(df)[7]),
           sigq = case_when(q < 0.05 ~ paste0("q<0.05"), q >= 0.05 ~ paste0("not sig")),
           sigq = as.factor(sigq))

(plasvdiet <- ggplot(res_asvs3, aes(x=ASV, y=est, color = Diet, shape = sigq)) +
        geom_hline(yintercept = 0, color = "grey40") +
        geom_point(position=position_dodge(-0.7)) +
        scale_shape_manual(values = c(21,19))+
        # scale_y_continuous(breaks = c(-2:8))+
        geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.7)) +
        theme_Publication()+
        theme(legend.position = "bottom")+
        labs(title = "Dietary factors and ASVs",
             x = "", y = "Difference (log10 counts) for SD increase in dietary factor", shape = "", color = "") +
        scale_color_cosmic(guide = "none") +
        coord_flip() +
        facet_wrap(~Diet, nrow = 2))
ggsave("results/glm_clr/lm_asvdiet.pdf", width = 11, height = 9)
ggsave("results/glm_clr/lm_asvdiet.svg", width = 11, height = 9)

## Step 4: location and ASVs
restotal <- c()
for (i in c(55:74)) {
    df$asv <- df[[i]]
    m0 <- lm(asv ~ Site, data = df)
    m1 <- lm(asv ~ Site + BMI, data = df)
    m2 <- lm(asv ~ Site + Fibre, data = df)
    m3 <- lm(asv ~ Site + Carbohydrates, data = df)
    m4 <- lm(asv ~ Site + SodiumInt, data = df)
    
    taxasv <- colnames(df)[i]
    m0 <- tidy(m0, conf.int=T)[2,]
    m1 <- tidy(m1, conf.int=T)[2,]
    m2 <- tidy(m2, conf.int = T)[2,]
    m3 <- tidy(m3, conf.int = T)[2,]
    m4 <- tidy(m4, conf.int = T)[2,]
    
    resRow <- cbind(taxasv, m0$estimate, m0$conf.low, m0$conf.high, m0$p.value,
                    m1$estimate, m1$conf.low, m1$conf.high, m1$p.value,
                    m2$estimate, m2$conf.low, m2$conf.high, m2$p.value,
                    m3$estimate, m3$conf.low, m3$conf.high, m3$p.value,
                    m4$estimate, m4$conf.low, m4$conf.high, m4$p.value)
    colnames(resRow) <- c("ASV", 
                          "m0-est", "m0-l95", "m0-u95", "m0-p", 
                          "m1-est", "m1-l95", "m1-u95", "m1-p",
                          "m2-est", "m2-l95", "m2-u95", "m2-p",
                          "m3-est", "m3-l95", "m3-u95", "m3-p",
                          "m4-est", "m4-l95", "m4-u95", "m4-p")
    restotal <- rbind(restotal, resRow)
    df$asv <- NULL
}

restotal <- as.data.frame(restotal)
restotal2 <- restotal %>% 
    mutate_at(c(2:21), as.character) %>% 
    mutate_at(c(2:21), as.numeric) %>% 
    mutate_at(c(2:4, 6:8, 10:12, 14:16, 18:20), afronden2) %>% 
    mutate(
        `m0-q` = p.adjust(`m0-p`, 'fdr'),
        `m1-q` = p.adjust(`m1-p`, 'fdr'),
        `m2-q` = p.adjust(`m2-p`, 'fdr'),
        `m3-q` = p.adjust(`m3-p`, 'fdr'),
        `m4-q` = p.adjust(`m4-p`, 'fdr')
    ) %>% 
    mutate_at(c(5,9,13,17,21), afronden3)
openxlsx::write.xlsx(restotal2, "results/glm_clr/lm_compl_diet.xlsx")

restotal3 <- restotal2 %>% 
    pivot_longer(c(2:26), names_to=c("model", "cat"), 
                 names_prefix="m", 
                 names_sep='-',
                 values_to="value") %>% 
    pivot_wider(names_from = cat, values_from = value) %>% 
    mutate(model = factor(model, levels = c("0", "1", "2", "3", "4"), 
                          labels = c("Unadjusted", "BMI", 
                                     "Fibre", "Carbohdyrates", "Sodium")),
           ASV = factor(ASV, levels = colnames(df)[(ncol(df)-19):ncol(df)]),
           ASV = fct_rev(ASV), 
           sigq = case_when(q < 0.05 ~ paste0("q<0.05"), q >= 0.05 ~ paste0("not sig")),
           sigq = as.factor(sigq))

(pltotaal <- ggplot(restotal3, aes(x=ASV, y=est, color=model, shape = sigq)) +
        geom_hline(yintercept = 0, color = "grey40") +
        geom_point(position=position_dodge(-0.7)) +
        scale_shape_manual(values = c(21,19))+
        # scale_y_continuous(breaks = c(-2:8))+
        geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.7)) +
        theme_Publication()+
        theme(legend.position = "bottom")+
        labs(title = "ASVs with differential abundance between geographical locations",
             x = "", y = "Difference (log10-transformed counts) for urban location", shape = "", color = "") +
        scale_color_cosmic() +
        coord_flip())
ggsave("results/glm_clr/lm_compl_diet.pdf", width = 9, height = 10)
ggsave("results/glm_clr/lm_compl_diet.svg", width = 9, height = 10)

## Step 4: diet and ASVs
restotal <- c()
for (i in c(55:60)) {
    df$asv <- df[[i]]
    taxasv <- colnames(df)[i]
    for(j in c(7)){
        df$diet <- scale(df[[j]])
        m0 <- lm(asv ~ diet, data = df)
        m1 <- lm(asv ~ diet + Site, data = df)
        dietname <- colnames(df)[j]
        m0 <- tidy(m0, conf.int=T)[2,]
        m1 <- tidy(m1, conf.int=T)[2,]
        resRow <- cbind(taxasv, dietname, m0$estimate, m0$conf.low, m0$conf.high, m0$p.value,
                        m1$estimate, m1$conf.low, m1$conf.high, m1$p.value)
        colnames(resRow) <- c("ASV", "Diet",
                              "m0-est", "m0-l95", "m0-u95", "m0-p", 
                              "m1-est", "m1-l95", "m1-u95", "m1-p")
        restotal <- rbind(restotal, resRow)
        df$diet <- NULL
    }
    df$asv <- NULL
}

restotal <- as.data.frame(restotal)
restotal2 <- restotal %>% 
    mutate_at(c(3:10), as.character) %>% 
    mutate_at(c(3:10), as.numeric) %>% 
    mutate_at(c(3:5, 7:9), afronden2) %>% 
    mutate(
        `m0-q` = p.adjust(`m0-p`, 'fdr'),
        `m1-q` = p.adjust(`m1-p`, 'fdr')
    ) %>% 
    mutate_at(c(6,10), afronden3)
openxlsx::write.xlsx(restotal2, "results/glm_clr/lm_compl_diet_s.xlsx")

restotal3 <- restotal2 %>% 
    pivot_longer(c(3:12), names_to=c("model", "cat"), 
                 names_prefix="m", 
                 names_sep='-',
                 values_to="value") %>% 
    pivot_wider(names_from = cat, values_from = value) %>% 
    mutate(model = factor(model, levels = c("0", "1"), 
                          labels = c("Unadjusted", "Site")),
           ASV = factor(ASV, levels = colnames(df)[(ncol(df)-19):ncol(df)]),
           ASV = fct_rev(ASV), 
           Diet = factor(Diet, levels = colnames(df)[c(32,34,35,37)]),
           sigq = case_when(q < 0.05 ~ paste0("q<0.05"), q >= 0.05 ~ paste0("not sig")),
           sigq = as.factor(sigq))

(pltotaal <- ggplot(restotal3, aes(x=ASV, y=est, color=model, shape = sigq)) +
        geom_hline(yintercept = 0, color = "grey40") +
        geom_point(position=position_dodge(-0.7)) +
        scale_shape_manual(values = c(21,19))+
        # scale_y_continuous(breaks = c(-2:8))+
        geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.7)) +
        theme_Publication()+
        theme(legend.position = "bottom")+
        labs(title = "ASVs with differential abundance between geographical locations",
             x = "", y = "Difference (log10-transformed counts) for SD increase in macronutrient group", shape = "", color = "") +
        scale_color_cosmic() +
        coord_flip()+
        facet_wrap(~Diet))
ggsave("results/glm_clr/lm_compl_diet_2.pdf", width = 9, height = 8)
ggsave("results/glm_clr/lm_compl_diet_2.svg", width = 9, height = 8)
