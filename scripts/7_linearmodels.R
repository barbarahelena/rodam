## Linear models
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(rio)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(phyloseq)
library(dplyr)
library(tidyr)
library(forcats)
library(patchwork)
library(compositions)
library(broom)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90, vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(0.7)),
                axis.text.x = element_text(angle = 0), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
    
} 

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
rodamclin <- readRDS("data/clinicaldata_pcdiet.RDS")
rodamotu <- as.data.frame(t(as(rodam_mb@otu_table, "matrix")))
tax <- readRDS("data/taxtable.RDS")

## Rural - Urban predictors
best_pred <- rio::import('rural_urban/output_XGB_class_rural_urban_2024_02_05__21-24-45/feature_importance.txt')

## Preparation for models
best_pred <- best_pred %>% arrange(-RelFeatImp) %>% slice(1:20)
dfmb <- rodamotu %>% dplyr::select(best_pred$FeatName)
dfmb2 <- as.data.frame(log10(dfmb+1))
dfmb2$ID <- rownames(dfmb2)
df <- left_join(rodamclin, dfmb2, by='ID') %>% filter(Site %in% c("Rural Ghana", "Urban Ghana"))
dim(df)
colnames(df)[(ncol(df)-19):ncol(df)] <- make.unique(tax$Tax[match(colnames(df)[(ncol(df)-19):ncol(df)], tax$ASV)])
dfmb2 <- df[(ncol(df)-19):ncol(df)]

## Regression models
res <- c()
for (i in c(58:77)) {
    df$asv <- df[[i]]
    m0 <- lm(asv ~ Site, data = df)
    m1 <- lm(asv ~ Site + Age + Sex + BMI + HT, data = df)
    m3 <- lm(asv ~ Site + Age + Sex + BMI + HT + DietPC1 + DietPC2, data = df)
    
    taxasv <- colnames(df)[i]
    m0 <- tidy(m0, conf.int=T)[2,]
    m1 <- tidy(m1, conf.int=T)[2,]
    m3 <- tidy(m3, conf.int = T)[2,]
    
    resRow <- cbind(taxasv, m0$estimate, m0$conf.low, m0$conf.high, m0$p.value,
                    m1$estimate, m1$conf.low, m1$conf.high, m1$p.value,
                    m3$estimate, m3$conf.low, m3$conf.high, m3$p.value)
    colnames(resRow) <- c("ASV", 
                          "m0-est", "m0-l95", "m0-u95", "m0-p", 
                          "m1-est", "m1-l95", "m1-u95", "m1-p",
                          "m3-est", "m3-l95", "m3-u95", "m3-p")
    res <- rbind(res, resRow)
    df$asv <- NULL
}

res <- as.data.frame(res)
res2 <- res %>% 
    mutate_at(c(2:13), as.character) %>% 
    mutate_at(c(2:13), as.numeric) %>% 
    mutate_at(c(2:4, 6:8, 10:12), afronden2) %>% 
    mutate(
        `m0-q` = p.adjust(`m0-p`, 'fdr'),
        `m1-q` = p.adjust(`m1-p`, 'fdr'),
        `m3-q` = p.adjust(`m3-p`, 'fdr')
    ) %>% 
    mutate_at(c(5,9,13), afronden3)
openxlsx::write.xlsx(res2, "results/glm/lm_ruralurban.xlsx")

res3 <- res2 %>% 
    pivot_longer(c(2:16), names_to=c("model", "cat"), 
                 names_prefix="m", 
                 names_sep='-',
                 values_to="value") %>% 
    pivot_wider(names_from = cat, values_from = value) %>% 
    mutate(model = factor(model, levels = c("0", "1", "3"), 
                          labels = c("Unadjusted", "Age, Sex, BMI, Hypertension", "+Diet")),
           ASV = factor(ASV, levels = colnames(df)[(ncol(df)-19):ncol(df)]),
           ASV = fct_rev(ASV),
           sigq = case_when(q < 0.05 ~ paste0("q<0.05"), q >= 0.05 ~ paste0("not sig")),
           sigq = as.factor(sigq))

(plurb <- ggplot(res3, aes(x=ASV, y=est, color=model, shape = sigq)) +
        geom_hline(yintercept = 0, color = "grey40") +
        geom_errorbar(aes(ymin=l95,ymax=u95, color = model), width=.3, position=position_dodge(-0.7)) +
        geom_point(position=position_dodge(-0.7)) +
        scale_shape_manual(values = c(21,19))+
        # scale_y_continuous(breaks = c(-2:8))+
        theme_Publication()+
        theme(legend.position = "bottom")+
        labs(title = "Rural - urban Ghana",
             x = "", y = "Difference (log10-transformed counts) for urban location", shape = "", color = "") +
        scale_color_manual(values = pal_cosmic()(7)[c(1,6,7)]) +
        coord_flip())

plots <- lapply(c(colnames(dfmb2)), distribution, df = dfmb2)
distr <- ggarrange(plotlist = plots, ncol = 1)
plots1 <- plurb + distr + plot_layout(ncol = 2, widths = c(5,1))
# ggsave("results/glm/lm_models.pdf", width = 9, height = 10)
# ggsave("results/glm/lm_models.svg", width = 9, height = 10)

## Preparation for models
best_pred <- best_pred %>% arrange(-RelFeatImp) %>% slice(1:20)
dfmb <- rodamotu %>% dplyr::select(best_pred$FeatName)
dfmb2 <- as.data.frame(log10(dfmb+1))
dfmb2$ID <- rownames(dfmb2)
df <- left_join(rodamclin, dfmb2, by='ID') %>% filter(Site %in% c("Rural Ghana", "Amsterdam"))
dim(df)
colnames(df)[(ncol(df)-19):ncol(df)] <- make.unique(tax$Tax[match(colnames(df)[(ncol(df)-19):ncol(df)], tax$ASV)])
dfmb2 <- df[(ncol(df)-19):ncol(df)]

## Regression models
res <- c()
for (i in c(58:77)) {
    df$asv <- df[[i]]
    m0 <- lm(asv ~ Site, data = df)
    m1 <- lm(asv ~ Site + Age + Sex + BMI + HT, data = df)
    m3 <- lm(asv ~ Site + Age + Sex + BMI + HT + DietPC1 + DietPC2, data = df)
    
    taxasv <- colnames(df)[i]
    m0 <- tidy(m0, conf.int=T)[2,]
    m1 <- tidy(m1, conf.int=T)[2,]
    m3 <- tidy(m3, conf.int = T)[2,]
    
    resRow <- cbind(taxasv, m0$estimate, m0$conf.low, m0$conf.high, m0$p.value,
                    m1$estimate, m1$conf.low, m1$conf.high, m1$p.value,
                    m3$estimate, m3$conf.low, m3$conf.high, m3$p.value)
    colnames(resRow) <- c("ASV", 
                          "m0-est", "m0-l95", "m0-u95", "m0-p", 
                          "m1-est", "m1-l95", "m1-u95", "m1-p",
                          "m3-est", "m3-l95", "m3-u95", "m3-p")
    res <- rbind(res, resRow)
    df$asv <- NULL
}

res <- as.data.frame(res)
res2 <- res %>% 
    mutate_at(c(2:13), as.character) %>% 
    mutate_at(c(2:13), as.numeric) %>% 
    mutate_at(c(2:4, 6:8, 10:12), afronden2) %>% 
    mutate(
        `m0-q` = p.adjust(`m0-p`, 'fdr'),
        `m1-q` = p.adjust(`m1-p`, 'fdr'),
        `m3-q` = p.adjust(`m3-p`, 'fdr')
    ) %>% 
    mutate_at(c(5,9,13), afronden3)
openxlsx::write.xlsx(res2, "results/glm/lm_ruralurban_ruralAms.xlsx")

res3 <- res2 %>% 
    pivot_longer(c(2:16), names_to=c("model", "cat"), 
                 names_prefix="m", 
                 names_sep='-',
                 values_to="value") %>% 
    pivot_wider(names_from = cat, values_from = value) %>% 
    mutate(model = factor(model, levels = c("0", "1","3"), 
                          labels = c("Unadjusted", "Age, Sex, BMI, Hypertension", "+Diet")),
           ASV = factor(ASV, levels = colnames(df)[(ncol(df)-19):ncol(df)]),
           ASV = fct_rev(ASV),
           sigq = case_when(q < 0.05 ~ paste0("q<0.05"), q >= 0.05 ~ paste0("not sig")),
           sigq = as.factor(sigq))

(plurb_Ams <- ggplot(res3, aes(x=ASV, y=est, color=model, shape = sigq)) +
        geom_hline(yintercept = 0, color = "grey40") +
        geom_point(position=position_dodge(-0.7)) +
        scale_shape_manual(values = c(21,19))+
        # scale_y_continuous(breaks = c(-2:8))+
        geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.7)) +
        theme_Publication()+
        theme(legend.position = "bottom")+
        labs(title = "Rural Ghana - Amsterdam",
             x = "", y = "Difference (log10-transformed counts) for Amsterdam location", shape = "", color = "") +
        scale_color_manual(values = pal_cosmic("hallmarks_light")(8)[c(1,6,7)]) +
        coord_flip())

plots <- lapply(c(colnames(dfmb2)), distribution, df = dfmb2)
distr <- ggarrange(plotlist = plots, ncol = 1)
plots2 <- plurb_Ams + distr + plot_layout(ncol = 2, widths = c(5,1))
# ggsave("results/glm/lm_models_rural_Amsterdam.pdf", width = 9, height = 10)
# ggsave("results/glm/lm_models.svg", width = 9, height = 10)


# Urban - Amsterdam model
best_pred <- rio::import('urban_ams/output_XGB_class_urban_ams_2024_02_21__14-01-19/feature_importance.txt')

## Preparation for models
best_pred <- best_pred %>% arrange(-RelFeatImp) %>% slice(1:20)
dfmb <- rodamotu %>% dplyr::select(best_pred$FeatName)
dfmb2 <- as.data.frame(log10(dfmb+1))
dfmb2$ID <- rownames(dfmb2)
df <- left_join(rodamclin, dfmb2, by='ID') %>% 
        filter(Site %in% c("Urban Ghana", "Amsterdam")) %>% 
        droplevels(.)
dim(df)
colnames(df)[(ncol(df)-19):ncol(df)] <- make.unique(tax$Tax[match(colnames(df)[(ncol(df)-19):ncol(df)], tax$ASV)])
dfmb2 <- df[(ncol(df)-19):ncol(df)]

## Regression models
res <- c()
for (i in c(58:77)) {
    df$asv <- df[[i]]
    m0 <- lm(asv ~ Site, data = df)
    m1 <- lm(asv ~ Site + Age + Sex + BMI + HT, data = df)
    m3 <- lm(asv ~ Site + Age + Sex + BMI + HT + DietPC1 + DietPC2, data = df)
    
    taxasv <- colnames(df)[i]
    m0 <- tidy(m0, conf.int=T)[2,]
    m1 <- tidy(m1, conf.int=T)[2,]
    m3 <- tidy(m3, conf.int = T)[2,]
    
    resRow <- cbind(taxasv, m0$estimate, m0$conf.low, m0$conf.high, m0$p.value,
                    m1$estimate, m1$conf.low, m1$conf.high, m1$p.value,
                    m3$estimate, m3$conf.low, m3$conf.high, m3$p.value)
    colnames(resRow) <- c("ASV", 
                          "m0-est", "m0-l95", "m0-u95", "m0-p", 
                          "m1-est", "m1-l95", "m1-u95", "m1-p",
                          "m3-est", "m3-l95", "m3-u95", "m3-p")
    res <- rbind(res, resRow)
    df$asv <- NULL
}

res <- as.data.frame(res)
res2 <- res %>% 
    mutate_at(c(2:13), as.character) %>% 
    mutate_at(c(2:13), as.numeric) %>% 
    mutate_at(c(2:4, 6:8, 10:12), afronden2) %>% 
    mutate(
        `m0-q` = p.adjust(`m0-p`, 'fdr'),
        `m1-q` = p.adjust(`m1-p`, 'fdr'),
        `m3-q` = p.adjust(`m3-p`, 'fdr')
    ) %>% 
    mutate_at(c(5,9,13), afronden3)
openxlsx::write.xlsx(res2, "results/glm/lm_urbanams.xlsx")

res3 <- res2 %>% 
    pivot_longer(c(2:16), names_to=c("model", "cat"), 
                 names_prefix="m", 
                 names_sep='-',
                 values_to="value") %>% 
    pivot_wider(names_from = cat, values_from = value) %>% 
    mutate(model = factor(model, levels = c("0", "1", "3"), 
                          labels = c("Unadjusted", "Age, Sex, BMI, Hypertension", "+Diet")),
           ASV = factor(ASV, levels = colnames(df)[(ncol(df)-19):ncol(df)]),
           ASV = fct_rev(ASV),
           sigq = case_when(q < 0.05 ~ paste0("q<0.05"), q >= 0.05 ~ paste0("not sig")),
           sigq = factor(sigq, levels = c("q<0.05", "not sig")))

res3$ASV <- fct_recode(res3$ASV, "Clostridium sensu stricto 1 cel/disp/saud" = 
                           "Clostridium sensu stricto 1 celatum/disporicum/saudiense")

(plams <- ggplot(res3, aes(x=ASV, y=est, color=model, shape = sigq)) +
        geom_hline(yintercept = 0, color = "grey40") +
        geom_point(position=position_dodge(-0.7)) +
        scale_shape_manual(values = c(19,21))+ #this one only has 1 level
        # scale_y_continuous(breaks = c(-2:8))+
        geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.7)) +
        theme_Publication()+
        theme(legend.position = "bottom")+
        labs(title = "Urban Ghana - Amsterdam",
             x = "", y = "Difference (log10-transformed counts) for Amsterdam location", 
             shape = "", color = "") +
        scale_color_manual(values = pal_cosmic("hallmarks_light")(8)[c(1,6,7)]) +
        coord_flip())

plots <- lapply(c(colnames(dfmb2)), distribution, df = dfmb2)
distr <- ggarrange(plotlist = plots, ncol = 1)
plots3 <- plams + distr + plot_layout(ncol = 2, widths = c(5,1))
# ggsave("results/glm/lm_models_urbanams.pdf", width = 9, height = 10)
# ggsave("results/glm/lm_models_urbanams.svg", width = 9, height = 10)


## Preparation for models
best_pred <- best_pred %>% arrange(-RelFeatImp) %>% slice(1:20)
dfmb <- rodamotu %>% dplyr::select(best_pred$FeatName)
dfmb2 <- as.data.frame(log10(dfmb+1))
dfmb2$ID <- rownames(dfmb2)
df <- left_join(rodamclin, dfmb2, by='ID') %>% filter(Site %in% c("Rural Ghana", "Amsterdam")) %>% 
    droplevels(.)
dim(df)
colnames(df)[(ncol(df)-19):ncol(df)] <- make.unique(tax$Tax[match(colnames(df)[(ncol(df)-19):ncol(df)], tax$ASV)])
dfmb2 <- df[(ncol(df)-19):ncol(df)]

## Regression models
res <- c()
for (i in c(58:77)) {
    df$asv <- df[[i]]
    m0 <- lm(asv ~ Site, data = df)
    m1 <- lm(asv ~ Site + Age + Sex + BMI + HT, data = df)
    m3 <- lm(asv ~ Site + Age + Sex + BMI + HT + DietPC1 + DietPC2, data = df)
    
    taxasv <- colnames(df)[i]
    m0 <- tidy(m0, conf.int=T)[2,]
    m1 <- tidy(m1, conf.int=T)[2,]
    m3 <- tidy(m3, conf.int = T)[2,]
    
    resRow <- cbind(taxasv, m0$estimate, m0$conf.low, m0$conf.high, m0$p.value,
                    m1$estimate, m1$conf.low, m1$conf.high, m1$p.value,
                    m3$estimate, m3$conf.low, m3$conf.high, m3$p.value)
    colnames(resRow) <- c("ASV", 
                          "m0-est", "m0-l95", "m0-u95", "m0-p", 
                          "m1-est", "m1-l95", "m1-u95", "m1-p",
                          "m3-est", "m3-l95", "m3-u95", "m3-p")
    res <- rbind(res, resRow)
    df$asv <- NULL
}

res <- as.data.frame(res)
res2 <- res %>% 
    mutate_at(c(2:13), as.character) %>% 
    mutate_at(c(2:13), as.numeric) %>% 
    mutate_at(c(2:4, 6:8, 10:12), afronden2) %>% 
    mutate(
        `m0-q` = p.adjust(`m0-p`, 'fdr'),
        `m1-q` = p.adjust(`m1-p`, 'fdr'),
        `m3-q` = p.adjust(`m3-p`, 'fdr')
    ) %>% 
    mutate_at(c(5,9,13), afronden3)
openxlsx::write.xlsx(res2, "results/glm/lm_urbanams_ruralAms.xlsx")

res3 <- res2 %>% 
    pivot_longer(c(2:16), names_to=c("model", "cat"), 
                 names_prefix="m", 
                 names_sep='-',
                 values_to="value") %>% 
    pivot_wider(names_from = cat, values_from = value) %>% 
    mutate(model = factor(model, levels = c("0", "1", "3"), 
                          labels = c("Unadjusted", "Age, Sex, BMI, Hypertension", "+Diet")),
           ASV = factor(ASV, levels = colnames(df)[(ncol(df)-19):ncol(df)]),
           ASV = fct_rev(ASV),
           sigq = case_when(q < 0.05 ~ paste0("q<0.05"), q >= 0.05 ~ paste0("not sig")),
           sigq = as.factor(sigq))
res3$ASV <- fct_recode(res3$ASV, "Clostridium sensu stricto 1 cel/disp/saud" = 
                           "Clostridium sensu stricto 1 celatum/disporicum/saudiense")

(plams_rural <- ggplot(res3, aes(x=ASV, y=est, color=model, shape = sigq)) +
        geom_hline(yintercept = 0, color = "grey40") +
        geom_point(position=position_dodge(-0.7)) +
        scale_shape_manual(values = c(19,21))+
        geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.7)) +
        theme_Publication()+
        theme(legend.position = "bottom")+
        labs(title = "Rural Ghana vs Amsterdam",
             x = "", y = "Difference (log10-transformed counts) for Amsterdam location", 
             shape = "", color = "") +
        scale_color_manual(values = pal_cosmic("hallmarks_light")(8)[c(1,6,7)]) +
        coord_flip())

plots <- lapply(c(colnames(dfmb2)), distribution, df = dfmb2)
distr <- ggarrange(plotlist = plots, ncol = 1)
plots4 <- plams_rural + distr + plot_layout(ncol = 2, widths = c(5,1))
# ggsave("results/glm/lm_models_urbanams_rural.pdf", width = 9, height = 10)
# ggsave("results/glm/lm_models_urbanams_rural.svg", width = 9, height = 10)

gplot1 <- ggarrange(plots1, plots2, ncol = 2, labels = c("A", "B"))
gplot1 <- annotate_figure(gplot1, top = text_grob("Best predictors rural vs urban Ghana"))
gplot2 <- ggarrange(plots3, plots4, ncol = 2, labels = c("C", "D"))
gplot2 <- annotate_figure(gplot2, top = text_grob("Best predictors urban Ghana vs Amsterdam"))
ggarrange(gplot1, gplot2, nrow = 2)
ggsave("results/glm/allmodels.pdf", width = 15, height = 10)

ggarrange(plots1, plots3, nrow = 1, ncol = 2, labels = c("A", "B"))
ggsave("results/glm/twomodels.pdf", width = 15, height = 10)
ggsave("results/glm/twomodels.svg", width = 15, height = 10)
