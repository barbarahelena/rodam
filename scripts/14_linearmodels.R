## Linear models
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(rio)
library(phyloseq)
library(dplyr)
library(tidyr)
library(forcats)
library(patchwork)

## Output folder
resultsfolder <- "results/glm"
dir.create(resultsfolder, showWarnings = FALSE)

## Functions
theme_Empty <- function() {
    theme_minimal() +
        theme(axis.title = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
              axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(), 
              axis.ticks.y = element_blank(),
              plot.margin=unit(c(4,0,0,0), "mm")
        )
}

distribution <- function(df, col){
    var <- df[[col]]
    pl <- ggplot(df) + 
        geom_density(aes(var), fill = pal_cosmic()(6)[6], color = pal_cosmic()(6)[6]) +      
        theme_Empty()
}

afronden2 <- function(x) return(as.numeric(format(round(x, 2),2)))
afronden3 <- function(x) return(as.numeric(format(round(x, 3),3)))

## Opening RODAM file
rodam_mb <- readRDS("data/phyloseq_sampledata.RDS")
rodamclin <- readRDS("data/clinicaldata.RDS")
dietids <- c("G0421", "G2545")
rodamclin <- rodamclin %>% mutate(
    TotalCalories = case_when(ID %in% dietids ~ NA, .default = TotalCalories),
    Proteins = case_when(ID %in% dietids ~ NA, .default = Proteins),
    Fat = case_when(ID %in% dietids ~ NA, .default = Fat),
    Fibre = case_when(ID %in% dietids ~ NA, .default = Fibre),
    Carbohydrates = case_when(ID %in% dietids ~ NA, .default = Carbohydrates),
    SodiumInt = case_when(ID %in% dietids ~ NA, .default = SodiumInt)
)
rodamotu <- as.data.frame(t(as(rodam_mb@otu_table, "matrix")))
tax <- readRDS("data/tax_table.RDS")
best_pred <- rio::import('rural_urban/output_XGB_class_rural_urban_2024_01_04__00-26-33/feature_importance.txt')

## Preparation for models
best_pred <- best_pred %>% arrange(-RelFeatImp) %>% slice(1:20)
dfmb <- rodamotu %>% select(best_pred$FeatName)
dfmb <- dfmb %>% mutate(across(everything(.), ~log10(.x + 1)))
dfmb$ID <- rownames(dfmb)
df <- left_join(rodamclin, dfmb, by='ID')
dim(df)
colnames(df)[(ncol(df)-19):ncol(df)] <- make.unique(tax$Tax[match(colnames(df)[(ncol(df)-19):ncol(df)], tax$ASV)])
dfmb2 <- df[(ncol(df)-19):ncol(df)]

## Regression models
res <- c()
for (i in c((ncol(df)-19):ncol(df))) {
    df$asv <- df[[i]]
    m0 <- lm(asv ~ Site, data = df)
    m1 <- lm(asv ~ Site + Age + Sex + BMI + CurrSmoking + BristolScale, data = df)
    m2 <- lm(asv ~ Site + Age + Sex + BMI + CurrSmoking + AntiHT + DMMed + AFung + LipidLowering, data = df)
    m3 <- lm(asv ~ Site + Age + Sex + BMI + CurrSmoking + AntiHT + DMMed + AFung + LipidLowering + TotalCalories + Fibre + Fat + Proteins + Carbohydrates + SodiumInt + log10(AlcoholIntake+0.01), data = df)
    
    taxasv <- colnames(df)[i]
    m0 <- tidy(m0, conf.int=T)[2,]
    m1 <- tidy(m1, conf.int=T)[2,]
    m2 <- tidy(m2, conf.int = T)[2,]
    m3 <- tidy(m3, conf.int = T)[2,]
    
    resRow <- cbind(taxasv, m0$estimate, m0$conf.low, m0$conf.high, m0$p.value,
                    m1$estimate, m1$conf.low, m1$conf.high, m1$p.value,
                    m2$estimate, m2$conf.low, m2$conf.high, m2$p.value,
                    m3$estimate, m3$conf.low, m3$conf.high, m3$p.value)
    colnames(resRow) <- c("ASV", 
                          "m0-est", "m0-l95", "m0-u95", "m0-p", 
                          "m1-est", "m1-l95", "m1-u95", "m1-p",
                          "m2-est", "m2-l95", "m2-u95", "m2-p",
                          "m3-est", "m3-l95", "m3-u95", "m3-p")
    res <- rbind(res, resRow)
    df$asv <- NULL
}

res <- as.data.frame(res)
res2 <- res %>% 
    mutate_at(c(2:17), as.character) %>% 
    mutate_at(c(2:17), as.numeric) %>% 
    mutate_at(c(2:4, 6:8, 10:12, 14:16), afronden2) %>% 
    mutate(
        `m0-q` = p.adjust(`m0-p`, 'fdr'),
        `m1-q` = p.adjust(`m1-p`, 'fdr'),
        `m2-q` = p.adjust(`m2-p`, 'fdr'),
        `m3-q` = p.adjust(`m3-p`, 'fdr')
    ) %>% 
    mutate_at(c(5,9,13,17), afronden3)
openxlsx::write.xlsx(res2, "results/glm/lm_ruralurban.xlsx")

res3 <- res2 %>% 
    pivot_longer(c(2:21), names_to=c("model", "cat"), 
                 names_prefix="m", 
                 names_sep='-',
                 values_to="value") %>% 
    pivot_wider(names_from = cat, values_from = value) %>% 
    mutate(model = factor(model, levels = c("0", "1", "2", "3"), 
                          labels = c("Unadjusted", "Age, Sex, BMI, Smoking, Bristol scale", 
                                     "+Medication", "+Diet")),
           ASV = factor(ASV, levels = colnames(df)[(ncol(df)-19):ncol(df)]),
           ASV = fct_rev(ASV),
           sigq = case_when(q < 0.05 ~ paste0("q<0.05"), q >= 0.05 ~ paste0("not sig")),
           sigq = as.factor(sigq))

(plurb <- ggplot(res3, aes(x=ASV, y=est, color=model)) +
        geom_hline(yintercept = 0, color = "grey40") +
        geom_point(position=position_dodge(-0.7)) +
        scale_shape_manual(values = 19)+
        # scale_y_continuous(breaks = c(-2:8))+
        geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.7)) +
        theme_Publication()+
        theme(legend.position = "bottom")+
        labs(title = "ASVs with differential abundance between geographical locations",
             x = "", y = "Difference (log10-transformed counts) for urban location", shape = "", color = "") +
        scale_color_cosmic() +
        coord_flip())

plots <- lapply(c(colnames(dfmb2)), distribution, df = dfmb2)
distr <- ggarrange(plotlist = plots, ncol = 1)
plots <- plurb + distr + plot_layout(ncol = 2, widths = c(5,1))
ggsave("results/glm/lm_models.pdf", width = 9, height = 10)
ggsave("results/glm/lm_models.svg", width = 9, height = 10)
