## Process results machine learning

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(stringr)
library(forcats)

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

source("scripts/functions.R")

path_true <- 'rural_urban_dich_genus/output_XGB_class_rural_urban_dich_genus_2024_02_28__13-03-39'
plot_feature_importance_class(path_true, 20)
plot_feature_importance_color_microbiome(path_true, 20)

path_true <- 'urban_ams_dich_genus/output_XGB_class_urban_ams_dich_genus_2024_02_28__13-03-41'
plot_feature_importance_class(path_true, 20)
plot_feature_importance_color_microbiome(path_true, 20)

## Data
mb <- readRDS("data/phyloseq_sampledata.RDS")
tax <- readRDS("data/taxtable_genus.RDS")
tax <- as.data.frame(tax)
tax$ASV <- rownames(tax)
tax$Tax <- str_c(tax$Genus, " spp.")
feat1 <- rio::import('rural_urban_dich_genus/output_XGB_class_rural_urban_dich_genus_2024_02_28__13-03-39/feature_importance.txt') %>% slice(1:10)
feat2 <- rio::import('urban_ams_dich_genus/output_XGB_class_urban_ams_dich_genus_2024_02_28__13-03-41/feature_importance.txt') %>% slice(1:10)
feats <- rbind(feat1, feat2) %>% filter(!duplicated(FeatName))
dim(feats)

clin <- readRDS("data/clinicaldata.RDS")

## Data Model Rural - urban
mat <- t(as(mb@otu_table, "matrix"))
mat <- mat[,feat1$FeatName]
mat <- ifelse(mat > 0, 1, 0)
colnames(mat) <- tax$Tax[match(colnames(mat), tax$ASV)]
dfmb <- as.data.frame(mat)
dfmb$ID <- rownames(mat)

df <- left_join(clin, dfmb, by = "ID") %>% 
    mutate(across(c(56:65), as.factor)) %>% 
    mutate(across(c(56:65), ~fct_recode(.x, "Present"="1", "Absent"="0")))

dflong <- df %>% select(ID, Site, 56:65) %>% pivot_longer(., 3:12, names_to = "ASV")
dfsum <- dflong %>% group_by(Site, ASV) %>% summarise(percentage = mean(value == "Present")*100)

(pl1 <- ggplot(data = dfsum, aes(x = ASV, group = Site, fill = Site, y = percentage)) +
        geom_bar(stat = "identity", position = "dodge") + 
        theme_Publication() +
        scale_fill_manual(values = pal_cosmic()(4)[c(2:4)]) +
        labs(y = "Presence %", x = "", fill = "", title = "Predictors rural - urban Ghana") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave(pl1, filename = "results/vanishingmicrobes_genus_1.pdf", width = 10, height = 5, )

## Data model urban-ams
mat <- t(as(mb@otu_table, "matrix"))
mat <- mat[,feat2$FeatName]
mat <- ifelse(mat > 0, 1, 0)
colnames(mat) <- tax$Tax[match(colnames(mat), tax$ASV)]
dfmb <- as.data.frame(mat)
dfmb$ID <- rownames(mat)

df <- left_join(clin, dfmb, by = "ID") %>% 
    mutate(across(c(56:65), as.factor)) %>% 
    mutate(across(c(56:65), ~fct_recode(.x, "Present"="1", "Absent"="0")))

dflong <- df %>% select(ID, Site, 56:65) %>% pivot_longer(., 3:12, names_to = "ASV")
dfsum <- dflong %>% group_by(Site, ASV) %>% summarise(percentage = mean(value == "Present")*100)

(pl2 <- ggplot(data = dfsum, aes(x = ASV, group = Site, fill = Site, y = percentage)) +
        geom_bar(stat = "identity", position = "dodge") + 
        theme_Publication() +
        scale_fill_manual(values = pal_cosmic()(4)[c(2:4)]) +
        labs(y = "Presence %", x = "", fill = "", title = "Predictors urban Ghana - Amsterdam") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave(pl2, filename = "results/vanishingmicrobes_genus_2.pdf", width = 10, height = 5)

ggarrange(pl1, pl2, nrow = 2, 
          common.legend = TRUE, 
          legend = "bottom",
          labels = c("A", "B"))
ggsave("results/vanishingmicrobes_genus.pdf", width = 10, height = 10)
